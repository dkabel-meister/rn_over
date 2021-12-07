// r_from_rstar contains a solver for rstar = r + 2M ln(r/2M - 1) for r
//
// File created on Jan 24, 2006
// Authors: Roland Haas and Peter Zimmerman
// 

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "r_from_rstar.h"
#include "helper.h"
#include "mathhelper.h"

//
// static double rstar_finder_Schw(const double r, void *params)
// -------------------------------------------------------------
// helper function to solve r0^* = r + 2*M*log(r/(2*M)-1) for r
// we have two variants used for r^* > 4M and r^* < 4M
// for r^* > 4M we use the equation above, below
// in fact we calulate:
// (r/2M - 1)*exp((r-r0^*)/2M) - 1  (1)
// which share to correct roots with the equation above
// inputs:
// r	  : guess for r
// params : (double *)r0^*
// results:
// the value of eq. (1)
// globals:
// M	  : black hole mass
// 
static double rstar_exp_finder_Schw(const double r, void *params)
{
    const double rstar0 = *(double *)params;

    return (r/(2*M) - 1)*exp((r-rstar0)/(2*M)) - 1;
}
// derivative of rstar_finder
static double rstar_exp_finder_Schw_df(const double r, void *params)
{
    const double rstar0 = *(double *)params;

    return r/(4*SQR(M))*exp((r-rstar0)/(2*M));
}   
// both at once
static void rstar_exp_finder_Schw_fdf(const double r, void *params, double *y, double *dy)
{
    const double rstar0 = *(double *)params;
    double Exp;

    Exp = exp((r-rstar0)/(2*M));
    *y  = (r/(2*M) - 1)*Exp - 1;
    *dy = r/(4*SQR(M))*Exp;
}
// same thing using a ln
static double rstar_ln_finder_Schw(const double r, void *params)
{
    const double rstar0 = *(double *)params;

    return r + 2*M*log(r/(2*M)-1) - rstar0;
}
// derivative of rstar_finder
static double rstar_ln_finder_Schw_df(const double r, void *params)
{
    double rstar0 = *(double *)params;

    return r/(r - 2*M);
    rstar0 += 0; // make gcc happy
}   
// both at once
static void rstar_ln_finder_Schw_fdf(const double r, void *params, double *y, double *dy)
{
    const double rstar0 = *(double *)params;

    *y  = r + 2*M*log(r/(2*M)-1) - rstar0;
    *dy = r/(r - 2*M);
}


//
// double r_from_rstar_Schw(double rstar)
// ------------------------------------------------------------------------
// converts rstar to r by numerically solving rstar = r + 2*M*ln(r/2M - 1)
// inputs:
// rstar : the tortoise coordinate r^*
// results:
// returns the corresponding Schwarzschild coordinate r
// globals:
// M				: black hole mass
// r_from_r_star_epsilon	: result is accepted if successive iterations 
// 				: are closer than this
// r_from_rstar_max_iter	: max. iterations allowable to find a solution
// 
// we stop the iteration if we do not find a solution after this many iterations: */
static const unsigned int r_from_rstar_max_iter = 100;
static const double r_from_rstar_epsilon = DBL_EPSILON;  // DBL_ESPILON ~ 2e-16

static double r_from_rstar_Schw(double rstar)
{
    double low_r, high_r; // initial guesses
    double r, r0;	  // curernt and previous guesses
    unsigned int iter;	  // current iteration number
    // gsl specific stuff
    int gsl_success;
    const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
    static gsl_root_fdfsolver *s = NULL;
    gsl_function_fdf fdf;
    gsl_error_handler_t *old_handler;

    // set the error handler to something useful
    gsl_error_header = "r_from_rstar";
    old_handler = gsl_set_error_handler(gsl_panic);

    // initialize the solver
    if(NULL == s)
        s = gsl_root_fdfsolver_alloc (T);

    // make our initial guess */
    // Note: we have that r < 4M <=> r^* < 4M
    if(rstar < 4*M) // the ln dominates */
    {
        // we use: r^* = r + 2M ln(r/2M-1)
        //  <=>    r = 2M (exp((r^* - r)/2M)+1)
        //  and  therefore since r>0:
        //  exp(r^*/2M) >  exp((r^*-r)/2M)
        //  and thus:
        //  r <  2M (exp(r^*/2M)+1)
        low_r  = 2*M; // that's the horizon, it cannot possibly be less
        high_r = 2*M*(exp(rstar/(2*M))+1);
        // initialize the solver 
        fdf.f = rstar_exp_finder_Schw;
        fdf.df = rstar_exp_finder_Schw_df;
        fdf.fdf = rstar_exp_finder_Schw_fdf;
    }
    else // r dominates */
    {
        // we use: r^* = r + 2M ln(r/2M-1)
        // since r > 4M the ln(...) is non-negative and thus
        // r^* > r
        // on the other hand we have:
        // ln(r/2M - 1) = ln(1+(r-4M)/2M) < (r-4M)/2M and thus
        // r > (r^*+4M)/2
        high_r = rstar;
        low_r  = (rstar+4*M)/2;
        // initialize the solver 
        fdf.f = rstar_ln_finder_Schw;
        fdf.df = rstar_ln_finder_Schw_df;
        fdf.fdf = rstar_ln_finder_Schw_fdf;
    }

    // for large and negative rstar, then r=2M this might happen
    if(GSL_SUCCESS == gsl_root_test_interval(low_r, high_r, r_from_rstar_epsilon, r_from_rstar_epsilon))
        return (high_r + low_r)/2;

    // initialize the solver
    fdf.params = &rstar;

    r = (low_r+high_r)/2;
    gsl_root_fdfsolver_set(s, &fdf, r);

    // iterate until we have found the solution
    iter = 0;
    do
    {
        iter++;
        r0 = r;

        if(GSL_SUCCESS != (gsl_success = gsl_root_fdfsolver_iterate(s)))
            break;

        r = gsl_root_fdfsolver_root(s);
        gsl_success = gsl_root_test_delta (r, r0, r_from_rstar_epsilon, r_from_rstar_epsilon);
    } while (GSL_CONTINUE == gsl_success && iter++ < r_from_rstar_max_iter);


    if(GSL_SUCCESS != gsl_success) // really should never occur
        panic("Could not find a acceptable solution for r after %d iteration\nBest guess: %g, interval:[%.19g,%.19g]", iter, r, r, r0);

    // restore error handler
    gsl_set_error_handler(old_handler);
    gsl_error_header = "default_header";

    return r;
}

//
// Helper function for Reissner Nordstrom 
//
double Rstar::rPlus(const double Q_bh) { return M+sqrt(M*M-Q_bh*Q_bh); }
double Rstar::rMinus(const double Q_bh) { return M-sqrt(M*M-Q_bh*Q_bh); }   
double Rstar::alphaPlus(double Q_bh) { return SQR(rPlus(Q_bh))/(rPlus(Q_bh)-rMinus(Q_bh)); }
double Rstar::alphaMinus(double Q_bh) { return SQR(rMinus(Q_bh))/(rPlus(Q_bh)-rMinus(Q_bh));}
double Rstar::logArg(double r, double Q_bh) { return log(r - rPlus(Q_bh)) * alphaPlus(Q_bh) - log(r - rMinus(Q_bh)) * alphaMinus(Q_bh); }


double Rstar::rstar_exp_finder_RN_f(const double r, const double rstar0, const double Q_bh)
{
    double diff=rPlus(Q_bh)-rMinus(Q_bh);
    double alpha_plus  = SQR(rPlus(Q_bh))/diff; 
    double alpha_minus = SQR(rMinus(Q_bh))/diff; 
    return exp(rstar0-r) - pow(r-rPlus(Q_bh), alpha_plus)/pow(r-rMinus(Q_bh), alpha_minus);
}
double Rstar::rstar_exp_finder_RN_df(const double r, const double rstar0, const double Q_bh)
{
    double diff = rPlus(Q_bh)-rMinus(Q_bh);
    double alpha_plus = SQR(rPlus(Q_bh))/diff; 
    double alpha_minus = SQR(rMinus(Q_bh))/diff; 
    return  - exp(rstar0-r) 
        - alpha_plus  * pow(r-rPlus(Q_bh),alpha_plus) / ( (r-rPlus(Q_bh)) * pow(r-rMinus(Q_bh),alpha_minus) ) 
        + alpha_minus * pow(r-rPlus(Q_bh),alpha_plus) / ( (r-rMinus(Q_bh))* pow(r-rMinus(Q_bh),alpha_minus) ); 
}
double Rstar::rstar_ln_finder_RN_f(const double r, const double rstar0, const double Q_bh)
{ 
    double div = rPlus(Q_bh)-rMinus(Q_bh);
    return r + SQR(rPlus(Q_bh))/div * log(r-rPlus(Q_bh)) - SQR(rMinus(Q_bh))/div * log(r-rMinus(Q_bh)) - rstar0;
}
double Rstar::rstar_ln_finder_RN_df(const double r, const double Q_bh)
{ 
    double div = rPlus(Q_bh)-rMinus(Q_bh); 
    return 1.0 + SQR(rPlus(Q_bh))/div/(r-rPlus(Q_bh)) - SQR(rMinus(Q_bh))/div/(r-rMinus(Q_bh));
}
double Rstar::rstar_from_r_RN(const double r, const double Q_bh)
{
    double div = rPlus(Q_bh)-rMinus(Q_bh);
    return r + SQR(rPlus(Q_bh))/div * log(r-rPlus(Q_bh)) - SQR(rMinus(Q_bh))/div * log(r-rMinus(Q_bh));
}                                                                              

// double get_r_switch(const double Q)
// ------------------------------------------------------------------------
// PURPOSE : use bisection to obtain the point where r = rstar
// inputs:
//  Q : black hole charge
// RESULTS:
//  returns the coordinate r=rstar
// globals:
// M				: black hole mass
// r_from_r_star_epsilon	: result is accepted if successive iterations 
// 				: are closer than this
// r_from_rstar_max_iter	: max. iterations allowable to find a solution
// -------------------------------------------------------------------------
double get_r_switch(const double Q)
{
        // 
        // bisect for the point where r=rstar
        //
        double bound_low  = 2.4; // as Q->1 the log arg asymptotes to 2.42...
        double bound_high = 4.1; // Q = 0 value plus bisection buffer
        double bound_mid;
        double log_arg=0;
        double dx=0;

        if (Rstar::logArg(bound_low,Q)*Rstar::logArg(bound_high,Q)>0)
            panic("Roots must be bracketed to do bisection\n");

        double r_switch = Rstar::logArg(bound_low,Q) < 0.0 ? dx=bound_high-bound_low : dx=bound_low-bound_high;

        do
        {
            log_arg = Rstar::logArg(bound_mid = r_switch + (dx *= 0.5), Q);

            if (log_arg <= 0.0) 
                r_switch = bound_mid;

        } while(fabs(dx)>4*r_from_rstar_epsilon || log_arg==0);
        
        return r_switch;
}
 //
// double r_from_rstar(const double rstar, const double Q_bh, const double r_switch)
// ----------------------------------------------------------------------------------
//  PURPOSE :converts rstar to r by numerically solving 
//      r^* = r + r_+^2/(r_+-r_-) \ln(r-r_+) - r_-^2/(r_+-r_-)\ln(r-r_-)
// INPUTS:
//  rstar : the tortoise coordinate r^*
//  Q_bh  : black hole charge
//  r_switch : place where r=r^*
// RESULTS:
//  returns the corresponding Schwarzschild coordinate r
// GLOBALS:
//  M				: black hole mass
//  r_from_r_star_epsilon	: result is accepted if successive iterations 
// 				: are closer than this
//  r_from_rstar_max_iter	: max. iterations allowable to find a solution
// 
// we stop the iteration if we do not find a solution after this many iterations 
// --------------------------------------------------------------------------------------
double r_from_rstar(const double rstar, const double Q_bh, const double r_switch)
{
    double r; 
    double r0;
    unsigned int iter;	  // current iteration number

    r = r_from_rstar_Schw(rstar);  // seed with Schwarzchild
    // 
    // Root solve for r using Newton's method 
    //
    if (Q_bh == 0)
    {
        return r;
    }
    else
    {
        if (rstar == r_switch)
        {
            r = rstar;
        }
        else if (rstar < r_switch)  
        {
            iter = 0;
            do
            {
                iter++;
                r0 = r;
                r = r0 - Rstar::rstar_exp_finder_RN_f(r0,rstar, Q_bh) / Rstar::rstar_exp_finder_RN_df(r0, rstar, Q_bh);
            } while (iter < r_from_rstar_max_iter && fabs(r-r0) > 1024*r_from_rstar_epsilon);

            //if (iter >= r_from_rstar_max_iter && fabs(r-r0) > 1024*r_from_rstar_epsilon)  
              //  panic("Could not find a acceptable solution for r for rstar=%g after %d iteration\nBest guess: %g, interval:[%.19g,%.19g]. Error=%.15e", iter, rstar, r, r, r0, fabs(r-r0));
        }
        else
        {
            // iterate until we have found the solution
            iter = 0;
            do
            {
                iter++;
                r0 = r;
                r = r0 - Rstar::rstar_ln_finder_RN_f(r0,rstar, Q_bh) / Rstar::rstar_ln_finder_RN_df(r0, Q_bh);
            } while (iter < r_from_rstar_max_iter && fabs(r-r0) > 1024*r_from_rstar_epsilon);

             //if (iter >= r_from_rstar_max_iter && fabs(r-r0) > 1024*r_from_rstar_epsilon)  
               // panic("Could not find a acceptable solution for r for rstar=%g after %d iteration\nBest guess: %g, interval:[%.19g,%.19g]", iter, rstar, r, r, r0);
        }
    }

    return r;
}                                     
