/********************************************************************************** 
*  field.cc
* ---------------------------------------------------------------------------------
* contains functions we use to evolve the field (secodn order accurate)
* we use the algorithms developed by Lousto and Price in arXiv:gr-qc/0503001
* and arXiv:gr-qc/9705071
*
* Authors: Peter Zimmerman and Roland Haas
* File Created on Nov 15, 2010
*********************************************************************************/   
#include <cmath>
#include <omp.h>


#include "field.h"
#include "header.h"
#include "helper.h"
#include "initial_data.h"
#include "mathhelper.h"
#include "particle.h"
#include "r_from_rstar.h"
#include "simdata.h"

#ifdef HAVE_IS_NAN
#    define have_isnan 1
#endif

#define do_copy 0

int low_pos = 0;

// the next array contains pointers that point to the areas currently 
// used to store the field for past, future and present 
// use void shuffle_fields(void) to resort them for the next timestep 
//
// these grids are staggered, that is the grids for eg. past and present or 
// offset by Delta_rstar/2. Schematically it looks like this
//
// future :  x x x x x x x x x x
// present: x x x x x x x x x x
// past   :  x x x x x x x x x x
// 
///               
double *psi[N_SLICES];
double **fields[] = {psi};

Field::Field(SimulationParameter_t fs)
{
    fsimParams = fs;

    particle = new Particle(fsimParams);  // Note: Particle gets the undecimated Delta_t from the commandline

    // ease typing for parameter keys
    q       = fsimParams.q;
    m       = fsimParams.m;
    Q       = fsimParams.Q;
    Delta_t = fsimParams.Delta_t;
    el      = fsimParams.el;
    el_max  = fsimParams.el_max; 
    t_min   = fsimParams.t_min;
    t_max   = fsimParams.t_max;
    t_switch = fsimParams.t_switch;

    // set ints
    N_cells             = 0;
    cell_left           = 0;
    cell_right          = 0;
    cell_inner_boundary = 0;
    cell_outer_boundary = 0;
    observer_cell       = 0;

    // set doubles 
    Delta_t    /= fsimParams.chop;
    Delta_rstar = 2.0*Delta_t;

    // excision vars
    do_excision = 1; // default to no excision
    x_t_excision = 5000; // defalt to 1000

}

Field::~Field()
{
    free(particle);
}
//
// void Field::initializeField(void)
//---------------------------------------------------------------------------------   
// PURPOSE:
//  sets up the field data on the inial timeslices using the data on the current
//  slice and the current first time derivative
//  assumes that the initial timeslice is even (t=0)
// INPUTS:
//  none
// RESULTS:
//  fill psi_past and psi_present with useful data
// GLOBALS:
//  psi
//---------------------------------------------------------------------------------   
void Field::initializeField(void)
{
    if (fsimParams.drdt_init > 0) 
        panic("Only radial infall supported.");
    /* 
     determine size of domain required
     
         t ^
           |                      3*t_max/2
     t_max + --------- \--------------\--
           |            \              \
           |             \              \  (t_{1/2}, r^*_{max})
           |              \             /  
           |               \           /
     t=0   | - - - - - - - - - - - - -/- - - -> r^*
                       ^    ^         ^
                 r^*_low   r^*_0     r^*_high 
    */

    // interpolate to get the time and phase space data on the initial slice
    //
    pos_tmin = particle->interpolate_position(fsimParams.t_min, &low_pos); 
    pos_tmax = particle->interpolate_position(fsimParams.t_max, &low_pos);     

    double r_low  = fsimParams.r_init;
    double r_high = fsimParams.r_init;

    // set rstar_low at top-left corner of future cone
    double rstar_low  = particle->rstar_from_r(r_low) - t_max;
    // set rstar_high to be r_init and 4 cells of right buffer
    double rstar_high = particle->rstar_from_r(r_high) + 4*Delta_rstar;
    
    //
    // field memory allocation
    //

    // number of cells needed to cover the future light cone on the left 
    // and the half out-going light cone on the right 
    //
    N_cells = (int)ceil(((rstar_high+t_max/2.) - (rstar_low-t_max)) / Delta_rstar);

    //
    // variables used in the inlines for getting rstar from cell
    // and vice-versa
    //
    // left edge of the computational domain
    rstar_min =  rstar_low;
    // right edge of the computational domain
    rstar_max = t_max/2. + rstar_high + Delta_rstar; 

    // get memory for fields on the slices
    for (int i = 0; i < DIM(psi); i++)
        psi[i] = new double[N_cells];

    // get memory for radial coordinates 
    for (int i = 0; i < DIM(r_from_cell); i++) 
        r_from_cell[i] = new double[N_cells];

    // fill in radial coordinates for all cells
    for (int c = 0; c <= N_cells; c++)
    {
        r_from_cell[EVEN][c] = r_from_rstar(rstar_from_cell(EVEN, c), Q, particle->get_rEqrstar());
        r_from_cell[ODD][c]  = r_from_rstar(rstar_from_cell(ODD,  c), Q, particle->get_rEqrstar());
    }

    //
    // make sure that observer falls on a cell of the un'decimated' grid
    //
    observer_cell = fsimParams.chop * (int)floor( (fsimParams.observer-rstar_min) /
            (fsimParams.chop*Delta_rstar) );

    //
    // grid boundary generation 
    //

    // implement the inner boundary at the RN event horizon 
    //
    for (cell_inner_boundary = N_cells-1 ; cell_inner_boundary > 0 ; cell_inner_boundary--)
    {
        // reduce size of domain to region where (numerically) r(rstar) > r_+
        if ( r_from_cell[ODD][cell_inner_boundary] <= (M+sqrt(M*M-Q*Q)) ) // should be r_+
            break;
    }
    if (cell_inner_boundary == 0)
        warn("Innermost radius is not numerically equal to r_plus.");

    // implement outer boundary
    // 
#if have_isnan
    if (!isnan(fsimParams.rstar_outer_boundary))
        cell_outer_boundary = cell_from_rstar(ODD, fsimParams.rstar_outer_boundary);
    else
        cell_outer_boundary = N_cells;      
#endif

    //
    //  Set the initial data
    //
    if (q==0)
    {
        dbg_msg("q=0 passed, assuming Gausian initial data", "");

        if (0 != fsimParams.drdt_init)
            panic("Please specify dr/dt = 0 for Gaussian data");

        double psi0_cell_left = psi0_gaussian(rstar_from_cell(EVEN, cell_left)); 

        if (fabs(psi0_cell_left) > 1e-12)
            panic("Gaussian peak is too wide, psi0 = %g does not vanish at (u0,v0)", psi0_cell_left);

        if ((cell_left - observer_cell) % fsimParams.chop != 0) // observer_cell is always on top of a coarse cell
            panic("Initial particle location is not a coarse cell. Please change it slightly.");    
        
        //
        // Gaussian initial data for q=0
        //

        // find the cells which bound the active numerical domain
        cell_right = cell_left = cell_from_rstar(EVEN, particle->rstar_from_r(fsimParams.r_init));

        for (int c = cell_left; c <= cell_outer_boundary; c++)
        {
            // set field initial data to gaussians
            psi[PRESENT][c] = psi0_gaussian( (c-cell_left)*Delta_rstar + rstar_from_cell(EVEN,c) );
            psi[PAST][c] = psi[FUTURE][c] = psi0_gaussian( (c-cell_left)*Delta_rstar + Delta_t + rstar_from_cell(ODD,c) ); // note the Delta_t shift in the argument for psi0

        }
    }
    else 
    {
        //
        // set initial 'active domain' which is position dependent
        //
        cell_left  = cell_from_r(EVEN, pos_tmin.r);  // left flanking cell 
        cell_right = cell_left + 4;                  // right flanking cell with a 1 cell buffer
       
        //
        // Unphysical (psi=0) initial data for sourced case
        //

        for (int c = 0; c <= N_cells ; c++)
        {
            // set field initial data to zero
            //
            if (fsimParams.drdt_init == 0 || t_switch != 0)
            {
                psi[ANCIENT][c] = psi_initial_data(particle->ancient.r, r_from_cell[EVEN][c]);  
                psi[PAST][c]    = psi_initial_data(particle->past.r,    r_from_cell[ODD ][c]); 
                psi[PRESENT][c] = psi_initial_data(particle->present.r, r_from_cell[EVEN][c]);
                psi[FUTURE][c]  = psi_initial_data(particle->future.r,  r_from_cell[ODD ][c]);  // here future position is still present since we haven't stepped yet
            }
            else
                psi[ANCIENT][c] = psi[PRESENT][c] = psi[PAST][c] = psi[FUTURE][c] = 0.0;
       }
    }

}

// void Field::evolveField()
//------------------------------------------------------
// PURPOSE : Evolves the field 
// INPUTS  : none
// RESULTS : fills the field array with values
// GLOBALS : psi, fields
//------------------------------------------------------
//
//----------------------------------------------------------------------
// PLEASE NOTE: due to the staggered grid (and Delta_rstar = 2 Delta_t)
// we have that:
// For a timestep where x_t is even:
//  \psi(t+h, r^*)   = psi[FUTURE][c]
//  \psi(t-h, r^*)   = psi[PAST][c] 
//  \psi(t, r^* + h) = psi[PRESENT][c+1] 
//  \psi(t, r^* - h) = psi[PRESENT][c]
// For a timestep where x_t is odd:
//  \psi(t+h, r^*)   = psi[FUTURE][c+1]
//  \psi(t-h, r^*)   = psi[PAST][C+1] 
//  \psi(t, r^* + h) = psi[PRESENT][c+1] 
//  \psi(t, r^* - h) = psi[PRESENT][c]
//
//
// sketch of grid structure:
//
//            0       1       2       3       4       5       6       7
//        0       1       2       3       4       5       6       7
// t= 7       x       x       x       x       x       x       x       x  
//                                  /   \.                             
// t= 6   x       x       x       x       x       x       x       x    
//                              /           \.                         
// t= 5       x       x       x       x       x       x       x       x  
//                          /                   \.                     
// t= 4   x       x       x       x       x       x       x       x    
//                      /                           \.                 
// t= 3       x       x       x       x       x       x       x       x  
//                  /                                   \.             
// t= 2   x       x       x       x       x       x       x       x    
//              /                                           \.         
// t= 1       x       x       x       x       x       x       x       x
//          /                                                  \.
// t= 0   x       x       x       x       x       x       x       x    
//                                                                  
// t=-1       x       x       x       x       x       x       x       x
//       0h  1h  2h  3h  4h  5h  6h  7h  8h  9h 10h
//
//----------------------------------------------------------------------
void Field::evolveField(void)
{
    // 
    // evolve field from time t=0 to t=t_max
    //
    for (x_t = 0; x_t*Delta_t <= t_max; x_t++)
    {
        int isodd = x_t&1;
        int cell_start, cell_stop;

        //================//
        // odd timestep   //
        //================// 
        if (isodd) 
        {
            if (cell_left < cell_inner_boundary) 
            {
                cell_start = cell_inner_boundary;

                psi[FUTURE][cell_start] = psi[PRESENT][cell_start];
            }
            else
                cell_start = cell_left;

            if (cell_right > cell_outer_boundary)
                cell_stop = cell_outer_boundary; // -2 in field_second_order
            else
                cell_stop = cell_right;

            // loop over all the cells on the current time slice
            //
#pragma omp parallel for firstprivate(cell_start, cell_stop)  
            for (int c = cell_start; c < cell_stop; c++)
            {
                double r  = r_from_cell[!isodd][c+isodd];     // radial coordinate for the middle of the cell
                double V0 = ( 1-Delta_t*Delta_t/2. * ((el*(el+1))/SQR(r) * (1-2*M/r+SQR(Q/r))) );

                // update the future value of the field
                psi[FUTURE][c+1] = -psi[PAST][c+1] +  V0 * (psi[PRESENT][c] + psi[PRESENT][c+1]);
            }/*-- End of OpenMP parallel region --*/
            
            // add source term
            //
            int c         = cell_from_rstar(isodd, particle->present.rstar);  // cell containing the particle 
            double r      = r_from_cell[!isodd][c+isodd];                     // radial coordinate at middle of sourced cell 
            double rstar0 = rstar_from_cell(!isodd,c+isodd);                  // radial coordinate at middle of sourced cell 

            double Vb = (el*(el+1))/SQR(r) * (1-2*M/r+SQR(Q/r));  // potential at middle/bottom of cell

            double A1,A2,A3,A4; // areas wrt u,v ie. A1(u,v) = 4 * A1(t, r^*) etc

            double t1, t2;           // times when the particle enter and leave the source cell
            double rstar1, rstar2;   // places where the particle enters and leave the source cell 

            // set (t1,rstar1) and (t2,rstar2) to the coordinates for which the p
            // particle enters and leave the cell
            //
            particle->findEnterleave(&t1, &rstar1, &t2, &rstar2, rstar0);           
            // use these coordinats to construct the areas for integration
            //
            particle->findAreas(rstar0, t1,rstar1, t2,rstar2, &A1, &A2, &A3, &A4);

            double source;
            integrateSource(t1,t2, rstar0, source);  // set the value of source 

            // update psi
            //
            psi[FUTURE][c+1] = -psi[PAST   ][c+1] * (1 + Vb/4*(A2-A3)) 
                              + psi[PRESENT][c+1] * (1 - Vb/4*(A4+A3)) 
                              + psi[PRESENT][c  ] * (1 - Vb/4*(A1+A3)) 
                              - 0.25*(1-Vb/4*A3)*source;

            // copy data on the excision boundary
#if do_copy 
            if (do_excision && (x_t>x_t_excision))
                psi[FUTURE][cell_left]=psi[PRESENT][cell_left];
#endif 
            // 
            // adjust domain to cover dom. of dep.
            //
            if (2*x_t*Delta_t < t_max) // grow the domain (this is really overkill for infall unless rdot(t=0) > 0 ) 
                cell_right++;
            else 
                (void)0;
        }
        //================//
        // even timestep  //
        //================// 
        else 
        {
            // we implement an inner boundary condition at the event horizon
            //
            if (cell_left < cell_inner_boundary) 
                cell_start = cell_inner_boundary; // no boundary value to impose here
            else
                cell_start = cell_left;

            if (cell_right >= cell_outer_boundary)
            {
                cell_stop = cell_outer_boundary;

                psi[FUTURE][cell_outer_boundary] = psi[PRESENT][cell_outer_boundary] 
                    + (psi[PRESENT][cell_outer_boundary] - psi[PAST][cell_outer_boundary-1]) * (1-2*Delta_t/fsimParams.rstar_outer_boundary);
            }
            else
                cell_stop = cell_right;

            // loop over all the cells on the current time slice
            //
#pragma omp parallel for firstprivate(cell_start, cell_stop)  
            for (int c = cell_start ; c < cell_stop ; c++)
            {
                double r  = r_from_cell[!isodd][c+isodd];
                double V0 = ( 1-Delta_t*Delta_t/2. * ((el*(el+1))/SQR(r) * (1-2*M/r+SQR(Q/r))) );

                // update the future value of the field
                psi[FUTURE][c+0] = -psi[PAST][c+0] + V0 * (psi[PRESENT][c] + psi[PRESENT][c+1]);
            }/*-- End of OpenMP parallel region --*/

            // add source term
            int c         = cell_from_rstar(isodd, particle->present.rstar);  // cell containing the particle 
            double r      = r_from_cell[!isodd][c+isodd];                     // radial coordinate at middle of sourced cell
            double rstar0 = rstar_from_cell(!isodd,c+isodd);                  // radial coordinate at middle of sourced cell

            double Vb = (el*(el+1))/SQR(r) * (1-2*M/r+SQR(Q/r));  // potential at middle/bottom of cell
            
            double A1,A2,A3,A4; // wrt u,v ie. A1(u,v) = 4 * A1(t, r^*) etc
            double t1,t2,rstar1,rstar2;

            particle->findEnterleave(&t1, &rstar1, &t2, &rstar2, rstar0);
            particle->findAreas(rstar0, t1,rstar1, t2,rstar2, &A1, &A2, &A3, &A4);

            double source;
            integrateSource(t1,t2, rstar0, source);

            // update psi for the sourced cell
            //
            psi[FUTURE][c+0] =  -psi[PAST   ][c+0] * (1 + Vb/4*(A2-A3)) 
                               + psi[PRESENT][c+1] * (1 - Vb/4*(A4+A3)) 
                               + psi[PRESENT][c  ] * (1 - Vb/4*(A1+A3))  
                               - 0.25*(1-Vb/4*A3) * source;

            //  compute flux at the center of the second to rightmost cell
            //  
            if (x_t%fsimParams.chop == 0 && fsimParams.flux_fh != NULL)
            {
                double t_buff = 40*Delta_t; // 20 cell buffer to reduce junk
                if ((x_t*Delta_t) > (t_max/2.+4.*Delta_rstar+t_buff))
                {
                    int flux_cell = cell_right-2;
                    double dpsi_dt = 1./(2.*Delta_t) * (psi[FUTURE][flux_cell] - psi[PAST][flux_cell]) ;
                    double flux    = 1./(4.*M_PI) * SQR(dpsi_dt) / (el*(el+1.));
                    double u = x_t*Delta_t - (rstar_from_cell(EVEN,flux_cell)+ 0.5*Delta_rstar);
                    fprintf(fsimParams.flux_fh, "%22.17e %22.17e %22.17e %22.17e\n", x_t*Delta_t,u, psi[PRESENT][flux_cell], flux);
                }
            }

            // 
            // Boundary conditions
            //
            if (do_excision && (x_t == x_t_excision))
            {
                int num_exc = 30;
                if (cell_from_rstar(EVEN, particle->present.rstar) <= (cell_left+num_exc))
                    panic("Excised domain includes the particle cell at time t=%f",x_t*Delta_t);
                else
                {
                    //int n = num_exc;
                    for (int n = cell_left+num_exc; n>0; --n) // excise cells 
                        psi[PRESENT][n] = psi[PAST][n] = psi[ANCIENT][n] = 0;
                    cell_left += num_exc;
                }

                //fprintf(stderr, "%22.17e %22.17e %22.17e\n", x_t*Delta_t, psi[PRESENT][cell_left], psi[PRESENT][cell_left+1]);
            }
#if do_copy
            else if (do_excision && (x_t > x_t_excision))
            {
                // copy the data along the ray v=constant for inner boundary
                psi[FUTURE][cell_left-1] = psi[PRESENT][cell_left];
                fprintf(stderr, "%22.17e %22.17e %22.17e\n", x_t*Delta_t, psi[PRESENT][cell_left], psi[PRESENT][cell_left+1]);
                // adjust domain to cover dom. of dep.
                cell_left--;
            }
#endif
            else
            {
                // adjust domain to cover dom. of dep.
                cell_left--;
                //fprintf(stderr, "%22.17e %22.17e %22.17e\n", x_t*Delta_t, psi[PRESENT][cell_left], psi[PRESENT][cell_left+1]);
            }

            if (2*x_t*Delta_t > t_max) 
                cell_right--;

            // output the data
            //
            if (x_t % fsimParams.chop == 0)
            {
                // dump particle data
                particle->dump();
                
                // dump mulitpoles at particle
                dump(x_t);

                // dump at observer 
                dumpAtObserver(x_t);
            }
        }

        // step the particle
        particle->step();  // interpolate to get the position of the particle on the next time slice
        // 
        // shuffle fields
        //
        for (int i = 0 ; i < DIM(fields) ; i++)
        {
            double *tmp        = fields[i][ANCIENT];
            fields[i][ANCIENT] = fields[i][PAST];
            fields[i][PAST]    = fields[i][PRESENT];
            fields[i][PRESENT] = fields[i][FUTURE];
            fields[i][FUTURE]  = tmp;
        }

        // status line
        if (!fsimParams.quiet && fmod(x_t*Delta_t,t_max/100) < Delta_t)
            fprintf(stderr, "%5d\r", int(x_t*Delta_t));
    }                      

}

// void Field::integrateSource(double t_1, double t_2, double x0, double & source)
//------------------------------------------------------
// PURPOSE: 
 //  integrates the RHS of the wave equation up to fourth order in h
// INPUTS:
//  t_1     : point where the particle enters the cell (full time)
//  t_2     : point where the particle leaves the cell
//  x0     : coordinate at the middle of the cell
//  source : source term
// RESULTS:
//  sets source :
//  integrates the source term using a Gauss-Legendre Method
//  and treats the boundary terms using Haas' re-tooling of the
//  way Price and Lousto split the source; i.e., writing it
//  as S(t,r) = f0(t)*G(t)*\delta(r-r0(t)) + f(r)*H(t)\delta^\prime(r-r0(t))
// GLOBALS:
//  psi
// NOTES : 
// ----------------------------------------------------------------
//
void Field::integrateSource(const double t_1, const double t_2, const double x0, double & source)
{
    double E = getParticle()->get_E();
    // maple code to get the acceleration

    // integrate using the DGL to calculate intermediate points
    // we use a four-point Gaussian quadrature rule (which should be enough since it has an error O(h^8))
    // coeffs from Abramowitz and Stegun p.916 (see p.10 11.07.2006)
    static const double x[4] = {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
    static const double w[4] = { 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454};
    double t0 = (t_1+t_2)/2; // (a+b)/2: offset
    double dt = (t_2-t_1)/2; // (b-a)/2: scale factor

    gsl_pos pos;
    

    source = 0;

    int lr; // +1 or -1 depending of enter and leave positions 

    double prefact_el = 4.*M_PI*q*sqrt((2.*el+1.)/4./M_PI);

    double rdotdot = 0;

    // 
    // Time Integral of the G(t,r0(t)) function
    //
    for (int i=0; i<DIM(x); i++)
    {
        pos = getParticle()->interpolate_position(t0+dt*x[i], &low_pos);

        double r0      = pos.r;  // worldline radial coordinate
        double f0      = 1. - 2.*M/r0 + SQR(Q/r0);  // metric function evaluated on the worldline r0(t)
        //double rdot2   = SQR(f0/(E-q*Q/r0)) * (SQR(E-q*Q/r0) - SQR(m)*f0); // (dr0/dt)^2
        
        double rdot2 = SQR(m*f0/(E-q*Q/r0)) * SQR(pos.rdot);

        
        if (t0+dt*x[i] > t_switch)
        {
            #include "accel.cc"
            rdotdot = d2r_dt2; // maple for d^2r0/dt^2
        }

        // update the integral value (only G since H/f is ind. of r)
        source += w[i] * (rdotdot/f0 + 2.0*rdot2/(SQR(f0*r0))*(SQR(Q)/r0 - M));
    }
           
        // apply scaling coming from the integral measure
        source *= dt; 

        // 
        // Boundary Terms  (The H's)
        //
        pos = getParticle()->interpolate_position(t_1, &low_pos); // position of the particle at time t_1 (= T_b in notebook)

        lr = (pos.rstar > x0 ? 1 : -1); // where is the particle wrt the center line ?

        // add first boundary term
        source += (lr + m*pos.rdot/(E-q*Q/pos.r)); // contribution from crossing at (t_1,r1)
        
        pos = getParticle()->interpolate_position(t_2, &low_pos); // position of the particle at time t_2 (= T_t in notebook) 

        lr = (pos.rstar > x0 ? 1 : -1);

        // add second boundary term
        source += (lr - m*pos.rdot/(E-q*Q/pos.r)); 

        // add the pre-factor and the 2 coming from the diffeo du*dv -> dt*dr^*
        source *= 2.*prefact_el;  


        //const double t_on_slice = x_t*Delta_t;
        //source *= transitionSource(t_on_slice);
        //fprintf(stderr, "%f %f %.15e %.15e %.15e %.15e \n", t_1, t_2, t_on_slice, pos.r, pos.rdot, source);
}                          

void Field::writeHeader(int argc, char **argv)
{

    Header header;
    header.initialize(argc, argv);
    
    // add to header
    header.add("Delta_t",    fsimParams.Delta_t);
    header.add("t_min",      fsimParams.t_min);
    header.add("t_max",      fsimParams.t_max);
    header.add("chop",       fsimParams.chop);
    header.add("r(t=0)",     fsimParams.r_init);
    header.add("dr/dt(t=0)", fsimParams.drdt_init);
    header.add("obsever",    fsimParams.observer);
    header.add("obs_cell",   observer_cell);
    header.add("el",         fsimParams.el);
    header.add("q",          fsimParams.q);
    header.add("m",          fsimParams.m);
    header.add("E",          particle->get_E());
    
    if (fsimParams.field_fh != NULL)
        header.write(fsimParams.field_fh);    
    if (fsimParams.multipole_fh != NULL)
        header.write(fsimParams.multipole_fh);    
    if (fsimParams.flux_fh != NULL)
        header.write(fsimParams.flux_fh);    
}

void Field::dumpAtObserver(int t_counter)
{

    if (fsimParams.field_fh == NULL)// || t_counter%fsimParams.chop != 0)
        return;
    if (t_counter&1)
        panic("Can only dump on even slices.");

    fprintf(fsimParams.field_fh, "%22.17e %22.17e\n", t_counter*Delta_t, psi[PRESENT][observer_cell]);
}

