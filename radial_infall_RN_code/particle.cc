// particle.cc
//
// File Created: Sat Sept 18, 2010
// File Copied for RN modifications: Thu Jan 6, 2011 
//
#include <math.h>
#include <float.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

#include "geometry.h"
#include "header.h"
#include "hermite.h"  
#include "helper.h"
#include "mathhelper.h"
#include "particle.h"
#include "r_from_rstar.h"

#define N_SLICES 4
#define GSL_DIMENSION 2
#define RKV_IND_RSTAR 0
#define RKV_IND_VRSTAR 1


// default constructor
Particle::Particle(SimulationParameter_t ps)
{
    psimParams = ps;
    // set fields to invalid data so that initialize has something to do
    header_written    = false;
    intersect_low_pos = -1;
    N_positions       = -1;
    internal_pos      = NULL;
    solver            = NULL;
    E                 = NAN;
    Q                 = psimParams.Q;
    q                 = psimParams.q;
    m                 = psimParams.m;
    r_switch          = get_r_switch(Q);
}

// default destructor
Particle::~Particle()
{
    if (internal_pos)
        free(internal_pos);

    if (solver)
        gsl_root_fdfsolver_free(solver);
}              

//
// geodesic_rhs
// ---------------------------------------------
// PURPOSE: 
//  calculates the right-hand side of the second-order evolution equation in rstar_dot
// INPUTS:
//  t	: current time (unused)
//  x	: a vector storing current rstar and v values
//  v	: the derivatives are stored here
// PARAMS  : (unused)
// GLOBALS:
//  M	: black hole mass
// ---------------------------------------------                
int geodesic_rhs(double t, const double x[], double v[], void *params)
{
    const class Particle *pt = (class Particle *)params; 


    const double r = r_from_rstar(x[RKV_IND_RSTAR], pt->Q, pt->r_switch);
    const double f = 1.0 - 2.0*M/r + SQR(pt->Q/r);

    v[RKV_IND_RSTAR]  = x[RKV_IND_VRSTAR];                                                                                 // dr^*/dt 
    v[RKV_IND_VRSTAR] = - SQR(pt->m)*f * (r*M*pt->E + pt->Q*((M-r)*pt->q - pt->E*pt->Q)) / CUBE(pt->E * r - pt->q*pt->Q);  // d^2r^*/dt^2

    return GSL_SUCCESS;
}
     
static const double step_particle_epsilon = 4*DBL_EPSILON;

void Particle::evolveParticle()
{
    int status;  // GSL status

    // set the data in traj_data
    header_written = false;

    int chop = psimParams.chop;
    psimParams.Delta_t /= chop;

    double t_min     = psimParams.t_min; 
    double t_max     = psimParams.t_max;  
    double t_switch  = psimParams.t_switch;

    double Delta_t   = psimParams.Delta_t; 
    double r_init    = psimParams.r_init;
    double drdt_init = psimParams.drdt_init;  // velocity in coordinate time
    
    double f_init = 1. - 2.*M/r_init + SQR(Q/r_init);  // f(r_init) in the RN metric 

    // energy of the particle 
    E = m * sqrt(f_init)/sqrt(1.-SQR(drdt_init/f_init)) + q*Q/r_init;  // needs to be in terms of coordinate time here

    // initial time counter
    x_t = 0;
     
    solver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_steffenson);  // http://www.gnu.org/software/gsl/manual/html_node/Initializing-the-Solver.html

    intersect_low_pos = 0;

    // set the gsl-error handler to something useful
    gsl_error_handler_t *old_handler;
    gsl_error_header = "initialize_particle";
    old_handler = gsl_set_error_handler(gsl_panic);

    //
    // set up the gsl ODE solver (Runge-Kutta-Fehlberg 4-5 algorithm)
    //
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step *stepper      = gsl_odeiv_step_alloc(T, GSL_DIMENSION);
    gsl_odeiv_control *control   = gsl_odeiv_control_y_new(step_particle_epsilon, step_particle_epsilon);
    gsl_odeiv_evolve *evolve     = gsl_odeiv_evolve_alloc(GSL_DIMENSION);
    gsl_odeiv_system sys         = {geodesic_rhs, NULL, GSL_DIMENSION, this};

    // get some memory for the position vector
    N_positions = (int)MAX((t_max-t_min)/Delta_t/(double)chop, N_SLICES); // just a guess, really
    if (NULL == (internal_pos = (struct internal_pos *)malloc(sizeof(gsl_pos) * N_positions)))
        panic("Out of memory to store the particle's position.");        
                                                                                                    
    // set up the initial conditions for rstar and vrstar
    double rstar_init  = rstar_from_r(r_init); // get initial r^* value from lookup tables
    double vrstar_init = drdt_init/f_init;     // get dr^*/dt, which is function of coordinate time, from dr/dt at t=0

    double x[GSL_DIMENSION];  // vector to store rstar(t) and vrstar(t)

    // set the initial values for the trajectory arrays
    x[RKV_IND_RSTAR ] = rstar_init;  // initial position
    x[RKV_IND_VRSTAR] = vrstar_init; // initial velocity    

    // time counters etc. 
    int stp;      // index for internal_pos, which saves the data
    double t;	  // current simulation time
    double t_end; // when to stop the simulation of the particle
    double h;	  // step size
    
    /* -------------------------*\
     * step into the past first *
    \*--------------------------*/

    stp   = 0;
    t     = 0.0;
    t_end = t_min - 3.0*Delta_t; 
    h     = -Delta_t; // h<0 since we're going into the past

    // 
    // evolve the system backwards 3 steps first
    //
    if (t_switch == 0)
    {
        do
        {
            status = gsl_odeiv_evolve_apply(evolve, control, stepper, &sys, &t, t_end, &h, x);

            if (stp == N_positions) // we need more space
            {
                N_positions += (int)MAX((t_max-t_min)/Delta_t/(double)chop, N_SLICES);
                if (NULL == (internal_pos = 
                            (struct internal_pos *)realloc(internal_pos, sizeof(gsl_pos) * N_positions)))
                    panic("Out of memory to store the particle's position.");
            }

            // store the quantities for later use
            internal_pos[stp].t      = t;
            internal_pos[stp].rstar  = x[RKV_IND_RSTAR ];
            internal_pos[stp].vrstar = x[RKV_IND_VRSTAR];

            // increment step
            stp++;
        } while (GSL_SUCCESS == status && t != t_end);   

        if (status != GSL_SUCCESS)
            panic("Evolution of the particle failed.");        

        // invert the order of the data points so that early points are first
        for (int j=0 ; j<stp/2 ; j++)
        {
            struct internal_pos swp;
            swp = internal_pos[j];
            internal_pos[j] = internal_pos[stp-1-j];
            internal_pos[stp-1-j] = swp;
        }                       
    }
    else  // t_switch is on
    {      
        // set the stp=0 element of internal_pos
        internal_pos[stp].t      = t_end;       // set the time to t_min - 3*Delta_t
        internal_pos[stp].rstar  = rstar_init;  // particle remains static for all t<t_switch
        internal_pos[stp].vrstar = 0.0;         // "
        stp++;
    }

    // save the position and velocity prior to the switch 
    internal_pos[stp].t      = t_switch - 10.0*Delta_t;  // ten time steps before switching on the particle, this helps the interpolation
    internal_pos[stp].rstar  = rstar_init;
    internal_pos[stp].vrstar = 0.0;        

    // increment stp
    stp++;
    
    // save the initial position & velocity data of the input phase space point
    internal_pos[stp].t      = t_switch;
    internal_pos[stp].rstar  = rstar_init;
    internal_pos[stp].vrstar = vrstar_init;

    /*---------------------------*\
     * step into the future next *
    \*---------------------------*/

    gsl_odeiv_evolve_reset(evolve);

    // fill in the initial position and velocity and step the simulation
    x[RKV_IND_RSTAR]  = rstar_init;
    x[RKV_IND_VRSTAR] = vrstar_init;

    t = t_switch; // reset time to switch time (which defaults to zero) 
    t_end = t_max + 2*Delta_t;  // end time
    h = Delta_t; // h>0 as we are moving futurewise

    // evolve the worldline in time
    do {
        status = gsl_odeiv_evolve_apply(evolve, control, stepper, &sys, &t, t_end, &h, x);

        if (stp == N_positions) // we need more space
        {
            N_positions += (int)MAX( (t_max-t_min)/Delta_t/chop, N_SLICES);

            if (NULL == (internal_pos = (struct internal_pos *)realloc(internal_pos, sizeof(gsl_pos) * N_positions)))
                panic("Out of memory to store the particle's position.");
        }

        // store the quantities for later use
        internal_pos[stp].t = t;
        internal_pos[stp].rstar  = x[RKV_IND_RSTAR ];
        internal_pos[stp].vrstar = x[RKV_IND_VRSTAR];
        
        stp++;

    } while (GSL_SUCCESS == status && t != t_end); // gsl guarantees that t==t_end finally

    if (status != GSL_SUCCESS)
        panic("Evolution of the particle failed.");

    // adjust the size of the position vector
    N_positions = stp;
    if (NULL == (internal_pos = (struct internal_pos *)realloc(internal_pos, sizeof(gsl_pos) * N_positions)))
        panic("Could not shrink the final position vector, this is odd.");                                   

    // now fill in the slices with something that is not utter nonsense (but still wrong)
    int low=0; // dummy 

    ancient = interpolate_position(-2*Delta_t, &low);
    past    = interpolate_position(-1*Delta_t, &low);
    present = interpolate_position( 0*Delta_t, &low);
    future  = interpolate_position( 1*Delta_t, &low);

    // tear down gsl ODE solver
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(stepper);

    // restore error handler
    gsl_set_error_handler(old_handler);
}                        

//
// get_bracket
// -----------------------------------------------------------------------------
// PURPOSE: finds the two entries in internal_pos that bracketing entries satisfying 
//  the condition cmp(entry) == 0 and returns a pointer to the first one
// INPUTS:
//  cmp : comparison function, is passed the current element and the param
//		: has to return 0 if the element is the one searched for, <0 if
//		: it is too low, and >0 if it is too high
//  param	: passed to cmp
//  low_pos	: used to speed up the search, cache this value between calls
// RESULTS:
//  updates low_pos to point to the lower one of the two entries
// GLOBALS:
//  internal_pos : all the particle's positions
//  N_positions	: DIM(internal_pos)
// -----------------------------------------------------------------------------
// 
void Particle::get_bracket(bracket_cmp cmp, void *param, int *low_pos)
{
    // hunt variables
    int high, low;
    int jump;

    // sanitize low_pos given by he caller
    if (*low_pos < 0)
        *low_pos = 0;
    else if (*low_pos > N_positions-2)
        *low_pos = N_positions-2;

    // hunt for an index with time value larger than the current 
    // simulation time.
    if (cmp(&internal_pos[*low_pos],param) <= 0) // search for high
    {
        jump = 1; // initial jump size
        low = *low_pos;
        high = low + 1;
        while (high < N_positions && cmp(&internal_pos[high],param) < 0)
        {
            jump *= 2;
            low = high;
            high += jump;
        }
    }
    else    // search for low_pos
    {
        jump = -1; // initial jump size
        high = *low_pos;
        low = high - 1;
        while (low >= 0 && cmp(&internal_pos[low], param) > 0)
        {
            jump *= 2;
            high = low;
            low += jump;
        }
    }

    // limit range to valid values
    if (high >= N_positions)
    {
        low = N_positions - 2;
        high = N_positions - 1;
    }
    else if (low < 0)
    {
        low = 0;
        high = 1;
    }

    while(high - low > 1) // zero in on the correct range
    {
        int med = (low + high)/2;

        if (cmp(&internal_pos[med], param) > 0)
            high = med;
        else
            low = med;
    }
    *low_pos = low;

    return;
}
                            
int time_cmp(struct Particle::internal_pos *ip, void *param)
{

    return (ip->t - *((double *)param) < 0 ? -1 : 1);
}                                                                                                               

    
//
// interpolate_position
//---------------------------------------------------------------------------------
// PURPOSE: Hermite interpolation with hunting
// INPUTS:
//  t		: time where we want to know the position
//  out		: storage for the result
//  low_pos	: used to speed up the search, cache this value between calls
// RESULTS:
//  sets out to the particle's position and four-velocity u^r at time t
// GLOBALS:
//  M		          : mass of the black hole
//---------------------------------------------------------------------------------
gsl_pos Particle::interpolate_position(const double t, int *low_pos)
{
    double rstarDot; // drstar/dt

    gsl_pos out;

    // interpolation variables
    double x[2]; // function values used for the Hermite interpolation
    double v[2]; // first derivatives used for the Hermite interpolation

    // aquire to entries which straddle the correct time
    get_bracket(time_cmp, (void *)&t, low_pos);

    // aquire the first derivatives at the upper and lower points
    x[0]= internal_pos[*low_pos+0].rstar;
    v[0]= internal_pos[*low_pos+0].vrstar;

    
    x[1] = internal_pos[*low_pos+1].rstar;
    v[1] = internal_pos[*low_pos+1].vrstar;

    // interpolate
    out.rstar = hermite_interpolate( 
            internal_pos[*low_pos+0].t, 
            internal_pos[*low_pos+1].t, 
            x[0], x[1],
            v[0], v[1], 
            t );

    rstarDot = hermite_interpolate_df( 
            internal_pos[*low_pos+0].t,
            internal_pos[*low_pos+1].t, 
            x[0], x[1], 
            v[0], v[1], 
            t );

    // calculate the derived quantities
    out.t = t;
    out.r = r_from_rstar(out.rstar, Q, r_switch);
    //
    // NOTE: we set pos.rdot to u^r = dr/d\tau
    //
    out.rdot  = (E - q*Q/out.r)/m * rstarDot; 

    return out;
}                        


//
// void step(void)
// ------------------------------------------------------------------------
// this function evolves the position of the particle by one timestep h
// it has to be called after each timestep
//
// INPUTS:
//  none
// RESULTS:
//  updates pos_future to reflect the particle's new position
// GLOBALS:
//  r_init            : initial position
//  drdt_init         : initial velocity in coordinate time
//  M		          : the black hole mass
//  Delta_t	          : temporal stepsize
//  x_t		          : the current time counter
//  internal_position : positions of the particle
// ------------------------------------------------------------------------
void Particle::step(void)
{
    // hunt variables
    static int low_pos = 0; // caches the index in internal_position

    // shuffle old values
    ancient = past;    
    past    = present;
    present = future;

    // step forward
    ++x_t;
    future = interpolate_position((x_t+1)*psimParams.Delta_t, &low_pos);
}                        


Header header;

// void dump_particle()
// ----------------------------------------------------------------------
// PURPOSE: writes the current position and velocity of the particle
// INPUTS:  filehandle (C) to write do
// RESULTS: writes out the contents of pos_present
// GLOBALS:
//  fh	    : output file name
//  header_written : bool for header
//  gsl_pos present : contains present data
//
void Particle::dump(void)
{  

    if (NULL == psimParams.particle_fh || x_t%psimParams.chop != 0 )
        return;
    // write header
    if (!header_written)
    {
        header.write(psimParams.particle_fh);
        header_written = true;
    }

    double f0_times_dtdtau = (E - q*Q/present.r)/m; 

    double drdt = m*(1.0-2.0*M/present.r+SQR(Q/present.r))/f0_times_dtdtau * present.rdot;

    fprintf(psimParams.particle_fh, "%f %21.18e %21.18e %21.18e %21.18e %21.18e\n", present.t, present.r, present.rstar, present.rdot, drdt, f0_times_dtdtau);
}                            


void Particle::writeHeader(int argc, char **argv)
{
    // add all data to the header
    header.initialize(argc, argv);
    header.add("Delta_t", psimParams.Delta_t);
    header.add("t_min",   psimParams.t_min);
    header.add("t_max",   psimParams.t_max);
    header.add("chop",    psimParams.chop);
    header.add("r0",      psimParams.r_init);
    header.add("drdt0",   psimParams.drdt_init);  
    header.add("E",       get_E());
}
