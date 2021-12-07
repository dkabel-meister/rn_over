// ----------------------------------------------------------------------
// extractor.cc
// ----------------------------------------------------------------------
// Extracts the field at the particle using the junction conditions
//
// File Created on Tuesday November 20, 2010
// Authors: Peter Zimmerman and Roland Haas
// ----------------------------------------------------------------------
#include "field.h"
#include "particle.h"
#include "helper.h"
#include "mathhelper.h"

void Field::getJumps(double *jumps, gsl_pos *pos)
{
    double prefact = 4*M_PI*q*sqrt((2*el+1)/(4*M_PI));
    double &jump_psi=jumps[0]; 

    jump_psi = prefact;
}              

// void Field::extrapolateRight(const double del, const double* const vals, 
//      double *field0, double *field_rstar)
// --------------------------------------------------------------------------------
// PURPOSE:
//  computes 4th order extrapolating polynomial on the right side of the particle.
// INPUTS:
//  del         : r^* - r^_0
//  vals        : field values to the right of the particle (i.e. psi[cell0+1]
//  field0      : right-extrapolated field at the particle
//  field_rstar : right-extrapoled rstar-derivative of field
// RESULTS:
//  fills field0 and field_rstar
// GLOBALS:
//  none
//---------------------------------------------------------------------------------   
void Field::extrapolateRight(const double del, const double* const vals, double *field0, double *field_rstar)
{
      *field0 = ((vals[1-1]/384.0-vals[2-1]/96.0+vals[5-1]/384.0+vals[3-1]/64.0-vals[4-1]/96.0)*del*del*
	del*del+((-3.0/16.0*vals[2-1]+vals[3-1]/4.0-7.0/48.0*vals[4-1]+5.0/96.0*vals[1-1]+vals[5-1]/32.0)*del*
	del*del+((-13.0/12.0*vals[2-1]+19.0/16.0*vals[3-1]-7.0/12.0*vals[4-1]+11.0/96.0*vals[5-1]+35.0/96.0
	*vals[1-1])*del*del+((-2.0*vals[2-1]+3.0/2.0*vals[3-1]-2.0/3.0*vals[4-1]+25.0/24.0*vals[1-1]+vals[5-1]/8.0)*
	del+vals[1-1]*Delta_t)*Delta_t)*Delta_t)*Delta_t)/(Delta_t*Delta_t*Delta_t*Delta_t);

      *field_rstar = ((-vals[1-1]/96.0-vals[5-1]/96.0+vals[2-1]/24.0-vals[3-1]/16.0+vals[4-1]/24.0)*del*
	del*del+((-5.0/32.0*vals[1-1]+9.0/16.0*vals[2-1]-3.0/4.0*vals[3-1]+7.0/16.0*vals[4-1]-3.0/32.0*vals[5-1]
	)*del*del+((-35.0/48.0*vals[1-1]+13.0/6.0*vals[2-1]-19.0/8.0*vals[3-1]+7.0/6.0*vals[4-1]-11.0/48.0*
	vals[5-1])*del+(2.0*vals[2-1]-3.0/2.0*vals[3-1]+2.0/3.0*vals[4-1]-25.0/24.0*vals[1-1]-vals[5-1]/8.0)*Delta_t)
	*Delta_t)*Delta_t)/(Delta_t*Delta_t*Delta_t*Delta_t);
}

// void Field::extrapolateLeft(const double del, const double* const vals, 
//      double *field0, double *field_rstar)
// --------------------------------------------------------------------------------
// PURPOSE:
//  computes 4th order extrapolating polynomial on the right side of the particle.
// INPUTS:
//  del         : r^* - r^_0
//  vals        : field values to the right of the particle (i.e. psi[cell0+1]
//  field0      : left-extrapolated field at the particle
//  field_rstar : left-extrapoled rstar-derivative of field
// RESULTS:
//  fills field0 and field_rstar
// GLOBALS:
//  none
//---------------------------------------------------------------------------------   
void Field::extrapolateLeft(const double del, const double* const vals, double *field0, double *field_rstar)
{
      *field0 = ((vals[1-1]/384.0-vals[2-1]/96.0+vals[5-1]/384.0+vals[3-1]/64.0-vals[4-1]/96.0)*del*del*
	del*del+((-vals[1-1]/32.0+7.0/48.0*vals[2-1]-vals[3-1]/4.0+3.0/16.0*vals[4-1]-5.0/96.0*vals[5-1])*del*
	del*del+((11.0/96.0*vals[1-1]-7.0/12.0*vals[2-1]+19.0/16.0*vals[3-1]-13.0/12.0*vals[4-1]+35.0/96.0*
	vals[5-1])*del*del+((2.0/3.0*vals[2-1]-3.0/2.0*vals[3-1]+2.0*vals[4-1]-vals[1-1]/8.0-25.0/24.0*vals[5-1])*del
	+Delta_t*vals[5-1])*Delta_t)*Delta_t)*Delta_t)/(Delta_t*Delta_t*Delta_t*Delta_t);

      *field_rstar = ((-vals[1-1]/96.0-vals[5-1]/96.0+vals[2-1]/24.0-vals[3-1]/16.0+vals[4-1]/24.0)*del*
	del*del+((3.0/32.0*vals[1-1]-7.0/16.0*vals[2-1]+3.0/4.0*vals[3-1]-9.0/16.0*vals[4-1]+5.0/32.0*vals[5-1])
	*del*del+((-11.0/48.0*vals[1-1]+7.0/6.0*vals[2-1]-19.0/8.0*vals[3-1]+13.0/6.0*vals[4-1]-35.0/48.0*
	vals[5-1])*del+(-2.0/3.0*vals[2-1]+3.0/2.0*vals[3-1]-2.0*vals[4-1]+vals[1-1]/8.0+25.0/24.0*vals[5-1])*Delta_t
	)*Delta_t)*Delta_t)/(Delta_t*Delta_t*Delta_t*Delta_t);
}

// void Field::interpolate(double *field_vals, int isodd, 
//      gsl_pos *pos, double *field0, double *field_rstar)
// --------------------------------------------------------------------------------
// PURPOSE:
//  computes 4th order extrapolating polynomials on each side of the particle given the
//  particle data and field data on each side
// INPUTS:
//  field_vals  : values of psi
//  isodd       : timestep parity specifier
//  pos         : position data
//  field0      : field
//  field_rstar : field derivative
// RESULTS:
//  fills field0 and field_rstar for both sides
// GLOBALS:
//  none
//---------------------------------------------------------------------------------   
void Field::interpolate(const double* const field_vals, const int isodd, 
        gsl_pos *pos, double *field0, double *field_rstar)
{
    int cell0;   // cell containing the particle
    double del;  // r^* - r^*_0

    cell0 = cell_from_rstar(isodd, pos->rstar); 

    del = rstar_from_cell(isodd, cell0) - pos->rstar;

    extrapolateLeft(del, &field_vals[cell0-4], &field0[0], &field_rstar[0]);

    del = rstar_from_cell(isodd, cell0+1) - pos->rstar;
    
    extrapolateRight(del, &field_vals[cell0+1], &field0[1], &field_rstar[1]);
}

// void Field::dump(const int t)
// --------------------------------------------------------------------------------
// PURPOSE:
//   Extract data at the worldline and write it to a file on timestep x_t
// INPUTS:
//   t : integer time
// RESULTS:
//  writes field data to file 
// GLOBALS:
//  field class
//---------------------------------------------------------------------------------    

void Field::dump(const int t)
{   
    // const double Right = 1;
    // const double Left =  0;

    if (NULL == fsimParams.multipole_fh || t%fsimParams.chop != 0)
        return;
    if (t&1)
        panic("Can only dump on even slices.");   

    double psi0[2];
    double psi_t[2], psi_rstar[2];
    double psi_future[2], psi_past[2]; 

    double jumps[1];
    

    // jumps to check my data
    getJumps(&jumps[0], &particle->present);
    
    // time derivative along world line
    interpolate(psi[FUTURE], !(t&1), &particle->future, psi_future, psi_rstar); // NB psi_rstar is a dummy
    interpolate(psi[PAST]  , !(t&1), &particle->past  , psi_past,   psi_rstar); // NB psi_rstar is a dummy
    
    // field and rstar derivative
    interpolate(psi[PRESENT], t&1, &particle->present, psi0, psi_rstar);

    //
    // regularize the field
    //
    double regularField[NUM_REGS];

    regularize(psi0[1], regularField); // compute (F_{tr}^R)_\ell using right extrapolated field points

    // partial time derivative (NB we need psi_rstar for this)
    for (int lr = 0 ; lr <= 1 ; lr++)
        psi_t[lr] = (psi_future[lr]-psi_past[lr])/(2*Delta_t)-particle->present.rdot/particle->get_E()*psi_rstar[lr];
    
    //fprintf(fsimParams.multipole_fh, "%.18e ",x_t*Delta_t);
    fprintf(fsimParams.multipole_fh, "%.18e ", particle->present.r );


    for (int reg_idx = 0; reg_idx < NUM_REGS; reg_idx++)
    {
        fprintf(fsimParams.multipole_fh, "%.18e ", regularField[reg_idx]); 
    }
    fprintf(fsimParams.multipole_fh, "\n");
}
