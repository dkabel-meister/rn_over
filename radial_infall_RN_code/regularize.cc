// regularize.cc
// ----------------------------------------------
// Computes and subtracts the singular component 
// of the retarded field via the regularization
// parameters.
//
// Author: Peter Zimmerman
// File Created on December 15th 2010
// File Copied for RN modifications: Fri Jan 14, 2011 
//

#include "field.h"
#include "mathhelper.h"
#include "helper.h"

     
// field::regularize(double retPsi, double &regField)
// ------------------------------------------------------
// PURPOSE:
//  computes the regularization parameters 
//  and sets the regular field 
// INPUTS: 
//  int el_max       :  max el value (used for tail)
//  double retPsi    : retarded psi on the current slice
//  double &regField : regular field being computed, this is r0^2*(F_{tr})_el
// RESULTS : 
//  returns r^2 rescaled field( F_{tr})_el
// GLOBALS:
//  M
// MODIFIED GLOBALS: none
// -------------------------------------------------------
void Field::regularize(double retPsi, double *regField) 
{
    // use sign for r>r0
    signed int sign_Delta = 1; 

    // get particle data
    const double r0 = particle->present.r;
    const double E  = getParticle()->get_E();

    double t = x_t/Delta_t;

//include for D coefficient 
#include "D_reg.cc"
    
    // A coefficient
    const double A = -sign_Delta;

    // B coefficient
    double B = 0.0;
    
    // C coefficient is zero

    // D coeff
    double D = 0.0;

  if ((x_t == 0 && fsimParams.drdt_init == 0) || (t < t_switch))
  {
      const double f0  = 1.0 - 2.0*M/r0 + SQR(Q/r0);
      const double df0_dr = 2.*M/SQR(r0)-2.*SQR(Q/r0)/r0;
      B = r0*df0_dr/(4.0*sqrt(f0))-sqrt(f0)/2.0; 
      D = 1.0/(f0*sqrt(f0)) * ( 9.0/16.0*(SQR(M)-SQR(Q))/SQR(r0) - 15.0/16.0*M*(SQR(M)-SQR(Q))/(r0*SQR(r0)) + 3.0/8.0*SQR(Q)*(SQR(M)-SQR(Q))/(SQR(r0)*SQR(r0)) );
  }
  else 
  {
        B = q*Q/(m*r0)- E/(2.0*m);
        D = Dreg; // from maple
  }
    double r2_retField_el = -q;  // el=0 value of r^2*(F_tr)_el

    if (el > 0)
        r2_retField_el = -sqrt((2.0*el+1.0)/4.0/M_PI) * retPsi;
    
    //
    // regularize 
    //
    regField[0] = r2_retField_el;                                 // retarded field
    regField[1] = regField[0]- q*(el+1.0/2.0)*A;             // remove A term 
    regField[2] = regField[1] - q*B;                         // remove A and B terms
    regField[3] = regField[2] -  q*D/((el-1/2.)*(el+3/2.));  // remove A, B, and D terms

    // add the tail term, which is used to speed up the 
    // convergence of the mode sum
    const double tail_coeff = 4.*(el_max+1.)/( (2.*el_max+1.)*(2*el_max+3.) );
    regField[4] = regField[2] + tail_coeff * q*D;
}


