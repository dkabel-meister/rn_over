// initial_data.cc
// contains functions describing the initial data of the system
//
// File created on Feb 01, 2006
// Last Change: Mon Sept 6 2010
// Authors: Roland Haas and Peter Zimmerman

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

#include "initial_data.h"
#include "mathhelper.h"
#include "helper.h"


/******************************************************
 * Funtions to compute the gaussian initial data      *
 ******************************************************/ 

// parameters describing the Gaussian peak (have to make sure it goes to zero
// fast enough near the initial location of the particle)
/*  WORKS WITH
 * ./generator --dt .25 --tmax 300 --el 3 --r0 100 --rdot0 0.0 --chop 4 --field_fh field.test --observer 100 --q 0 
 */
double gaussian_peak = 200;
double gaussian_width = 10;
// double psi0_gaussian(double r)
// ---------------------------------------------------------------
// PURPOSE: 
//  Generates gaussian initial data for the field psi
// INPUTS:
//  rstar : Tortoise coordinate where we want to evaluate the field
// RESULTS:
//  value of the field at the given position
// GLOBALS:
//  M : the black hole mass
//-----------------------------------------------------------------
double psi0_gaussian(const double rstar)
{
    return exp(-(rstar-gaussian_peak)*(rstar-gaussian_peak)/(2.*gaussian_width*gaussian_width));
}


//
// Funtions to compute the static charge initial data 
//

InitialData::InitialData() { ; }

// double legendre_Pl(const int, const double)
// ----------------------------------------------------
// PURPOSE: Calculates the legendre polynomials using 
//  recursion, but is not limited to arguments in the range
//  |x| <= 1, as in the GSL implementation
// INPUTS: 
//  l : multipole index
//  x : argument
// RESULTS: 
//  return P_{\ell}(x)
// GLOBALS: none
// -----------------------------------------------------
double InitialData::legendre_Pl(const int l, const double x)
{
    if (l==0)
    {
        return 1.0;
    }
    else if (l==1)
    {
        return x;
    }
    else
    {
        /* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */
        double p_ellm2 = 1.0;    /* P_0(x) = P_{l-2} */
        double p_ellm1 = x;      /* P_1(x) = P_{l-1} */
        double p_ell = p_ellm1;

        static int ell;

        for (ell=2; ell <= l; ell++)
        {
            p_ell = (x*(2*ell-1)*p_ellm1 - (ell-1)*p_ellm2) / ell;
            p_ellm2 = p_ellm1;
            p_ellm1 = p_ell;
        }

        return p_ell;
    }
}

// double legendre_Q0(const double)
// ----------------------------------------------------
// PURPOSE: Calculates the l=0 legendre function of 
//  second kind  rcursion, 
//  but is not limited to arguments in the range
//  |x| <= 1, as in the GSL implementation
// INPUTS: 
//  x : argument
// RESULTS: 
//  return Q_{0}(x)
// GLOBALS: none
// -----------------------------------------------------
double InitialData::legendre_Q0(const double x)
{
    if (x <= -1.0 || x == 1.0)
        panic("Invalid argument to legendreQ.\n");

    if (x<1.0)
        return 0.5*log((1.0+x)/(1.0-x));
    else
        return 0.5 * log((x+1.0)/(x-1.0));
}

// static double legendre_Q1(const double)
// ----------------------------------------------------
// PURPOSE: Calculates the l=1 legendre function of 
//  second kind  rcursion, 
//  but is not limited to arguments in the range
//  |x| <= 1, as in the GSL implementation
// INPUTS: 
//  x : argument
// RESULTS: 
//  return Q_{1}(x)
// GLOBALS: none
// -----------------------------------------------------
double InitialData::legendre_Q1(const double x)
{
    if (x <= -1.0 || x == 1.0)
        panic("Invalid argument to legendreQ.\n");

    if (x<1.0)
        return  0.5 * x * (log((1.0+x)/(1.0-x))) - 1.0;
    else
        return  0.5 * x * log((x+1.0)/(x-1.0)) - 1.0;
}


// static double legendre_Ql(const int, const double)
// ----------------------------------------------------
// PURPOSE: Calculates legendre polynomials of the
//  second kind using routines from GSL
// INPUTS: 
//  l : multipole index
//  x : argument
// RESULTS: 
//  return Q_{\ell}(x)
// GLOBALS: none
// -----------------------------------------------------   
double InitialData::legendre_Ql(const int l, const double x)
{
        if (l==0)
            return InitialData::legendre_Q0(x);
        else if (l==1)
            return InitialData::legendre_Q1(x);
        else
        {
            gsl_sf_result res;

            /* use GSL */
            int status = gsl_sf_legendre_Ql_e(l, x, &res);

            if (status == GSL_SUCCESS)
                return res.val;
            else
                // TODO: use gsl error handler
                panic("GSL failed to compute the legendre Q function.\n");
        }

}
// double phi_initial_data(const int, const double, const double)
// -----------------------------------------------------------------
// PURPOSE: Calculates the static initial data using the mode 
//  decomposition in Eric's notes.
// INPUTS: 
//  el : multipole index
//  r0 : source point on the initial slice
//  r  : field point on the initial slice
// RESULTS: 
//  returns psi := phi1_l (m=0)
// GLOBALS: 
//  M : mass of the black hole
// -----------------------------------------------------
double Field::psi_initial_data(const double r0, const double r)
{
    double phi = 0;

    if (el < 0)
        panic("Can not compute initial data if el<0.\n");

    // compute \Phi_{<} and \Phi_{>} for \ell=0
    if (el == 0)
    {
        if (r<r0)
            phi =  0.0;
        else
            phi = sqrt(4.0*M_PI*q);
    }
    else  // \ell > 0
    {
        const double lambda = M;
        const double mu = sqrt(SQR(M)-SQR(Q));

        const double x0 = (r0 - lambda)/mu;
        const double x = (r - lambda)/mu;

        const double pre = 4.*M_PI*q * el/(el+1.) * sqrt((2.*el+1)/(4.*M_PI)) * 1./(x0+lambda/mu);

        if (r<r0)
        {
            phi = pre*(x0*InitialData::legendre_Ql(el,x0)-InitialData::legendre_Ql(el-1,x0)) *
            ( (el*x + (el+1.)*lambda/mu)*InitialData::legendre_Pl(el,x) + InitialData::legendre_Pl(el-1,x) ); // Eq 4.60a   
        }
        else if (r>r0)
        {
            phi = pre*(x0*InitialData::legendre_Pl(el,x0)-InitialData::legendre_Pl(el-1,x0)) *
            ( (el*x + (el+1.)*lambda/mu)*InitialData::legendre_Ql(el,x) + InitialData::legendre_Ql(el-1,x) );  // Eq 4.60b
        }
        else
            panic("Can not compute initial data at the particle r=r0");
    }
            return phi;

}





