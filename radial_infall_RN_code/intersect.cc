// intersect.cc
// contains functions calculating intersections between the world line and a straight line
//
// File created on Feb 24, 2006
// Author: Roland Haass Peter Zimmerman

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "geometry.h"
#include "helper.h"
#include "hermite.h"
#include "intersect.h"
#include "mathhelper.h"
#include "particle.h"

// local helper functions 
static double part_position_intersect(double t, void *vp);
static double part_position_derivative(double t, void *vp);
static void part_position_fdf(double t, void *vp, double *f, double *df);

//
// part_position_intersect
// --------------------------------------------------------------------------------
// PURPOSE: 
//  function to calculate r_p^*(t) - r^*_{edge}(t)
//  where r^*_{edge} is on the straight line forming an edge of the cell
//  the particle passes through. Which edge it is is determined
//  by rstar0 and slope in *params (see struct part_position_params)
// INPUTS:
//  t      : time
//  params : points to a struct positionParams describing the edge
// RESULTS:
//  returns r_p^*(t) - r^*_{edge}(t)
// GLOBALS:
//  M : black hole mass
// NOTE: GSL requires the pointer to void in the function's arguments which is why
//  there is the cast below.
// ----------------------------------------------------------------------------------
static double part_position_intersect(double t, void *vp)
{
    struct positionParams *params = (struct positionParams *)vp;
    // interpolate to get r* 
    double rstar = hermite_interpolate(
			params->low.t, params->high.t, 
			params->low.rstar, params->high.rstar, 
			params->low.vrstar, params->high.vrstar, t );   

    return rstar - (params->rstar0 + params->slope*(t-params->t0));
}

//
// static double part_position_derivative(double t, void *vp)
// ----------------------------------------------------------------------
// PURPOSE: 
//  function to calculate dr_p^*/dt - dr^*_{edge}/dt
// INPUTS:
//  t      : time
//  params : points to a struct part_position_params describing the edge
// RESULTS:
//  returns dr_p^*/dt - slope
// GLOBALS:
//  M : black hole mass
// ----------------------------------------------------------------------

static double part_position_derivative(double t, void *vp)
{
    struct positionParams *params = (struct positionParams *)vp;

    double vrstar = hermite_interpolate_df(
			params->low.t, params->high.t, 
			params->low.rstar, params->high.rstar, 
			params->low.vrstar, params->high.vrstar, t );    

    return vrstar - params->slope;
}

// static void part_position_fdf(double t, void *void_params, double *f, double *df)
// -----------------------------------------------------------------------------------
// PURPOSE: 
//  function to calculate dr_p^*/dt - dr^*_{edge}/dt
// INPUTS:
//  t      : time
//  params : points to a struct part_position_params describing the edge
// RESULTS:
//  returns dr_p^*/dt - slope
// GLOBALS:
//  M : black hole mass
// NOTE: GSL requires the pointer to void in the function's arguments
// -----------------------------------------------------------------------------------

static void part_position_fdf(double t, void *vp, double *f, double *df)
{
    struct positionParams *params = (struct positionParams *)vp;
    
    double  rstar = hermite_interpolate(
			params->low.t, params->high.t, 
			params->low.rstar, params->high.rstar, 
			params->low.vrstar, params->high.vrstar, t );

    double vrstar = hermite_interpolate_df(
			params->low.t, params->high.t, 
			params->low.rstar, params->high.rstar, 
			params->low.vrstar, params->high.vrstar, t );    

    *f  = rstar - (params->rstar0 + params->slope*(t-params->t0));
    *df = vrstar - params->slope;
}


int lightcone_cmp(struct Particle::internal_pos *ip, void *vp)
{
    struct positionParams *par = (struct positionParams *)vp;
    double rstar;

    rstar = ip->rstar;

    return (((-par->slope)*(rstar - (par->rstar0 + par->slope*(ip->t-par->t0)))) < 0 ? -1 : 1);
}      

static const double find_enterleave_epsabs = 16*DBL_EPSILON;
static const double find_enterleave_epsrel = 128*DBL_EPSILON;
static const int find_enterleave_max_iter = 10;

//
// void Particle::intersect(double t0, double t1, double rstar0, double slope, double *t, double *rstar)
// -------------------------------------------------------------------------------------------------------
// PURPUSE: 
//  Uses the Steffenson algorithm for the gsl library to intersect the particle's trajectory from 
//  0 to Delta_t with the straight line rstar0+slope*t
// INPUTS:
//  t0	    : min time (either past or present)
//  t1	    : max time (either present or future)
//  rstar0  : vertex of the cell at t=t0 (left/right or bottom)
//  slope   : slope of the line forming the edge
//  t,rstar : will hold the coordinates of the intersection
// RESULTS:
//  updates t, rstar to hold the coordinates of the intersection
// -------------------------------------------------------------------------------------------------------
void Particle::intersect(const double t0, const double t1, const double rstar0, const double slope, double *t, double *rstar)
{ 
    double Delta_t = psimParams.Delta_t;

    int iter;				// number of iterations done so far
    double t_old=NAN;			// previous iterations 't' value

    struct positionParams params;	// extra arguments for 'f'

    gsl_function_fdf f;			// holds the function we want to set to zero

    // set the gsl error handler to something useful
    gsl_error_handler_t *old_handler;
    gsl_error_header = "find_enterleave_intersect";
    old_handler      = gsl_set_error_handler(gsl_panic);

    // initialize the difference function
    f.f   = part_position_intersect;
    f.df  = part_position_derivative;
    f.fdf = part_position_fdf;

    // set params
    params.t0     = t0;
    params.rstar0 = rstar0;
    params.slope  = slope;

    // interpolation variables
    double x[2];	// function values used for the Hermite interpolation
    double v[2];	// first derivatives used for the Hermite interpolation

    // set intersect_low_pos to point to the first entry in params such that
    // lightcone_cmp(entry)==0
    get_bracket(lightcone_cmp, (void *)&params, &intersect_low_pos);  

    // acquire the first derivatives at the upper and lower points
    x[0] = internal_pos[intersect_low_pos+0].rstar;
    v[0] = internal_pos[intersect_low_pos+0].vrstar;
    x[1] = internal_pos[intersect_low_pos+1].rstar;
    v[1] = internal_pos[intersect_low_pos+1].vrstar;

    params.low.t      = internal_pos[intersect_low_pos].t;
    params.low.rstar  = x[0];
    params.low.vrstar = v[0];

    params.high.t      = internal_pos[intersect_low_pos+1].t;
    params.high.rstar  = x[1];
    params.high.vrstar = v[1];

    f.params = &params;

    // initialize the solver to find t
    gsl_pos pos0[2];

    int low_pos = intersect_low_pos;

    pos0[0] = interpolate_position(t0, &low_pos);
    pos0[1] = interpolate_position(t1, &low_pos);

    // intersect_line is called intersect in odd and defined in geometry.h
    // this initializes the solver to find t
    Geometry::intersect_line(pos0[0].t, pos0[0].rstar, pos0[1].t, pos0[1].rstar, pos0[0].t,
            rstar0, pos0[0].t+Delta_t, rstar0+slope*Delta_t, t, rstar);  

    if (GSL_SUCCESS != gsl_root_fdfsolver_set(solver, &f, *t))
        panic("Could not initialize solver to find enter/leave time");

    // now iterate until we find t
    iter = 0;
    do
    {
        if (iter++ > find_enterleave_max_iter)
        {
            panic("Could not find the intersection of the world line and the cell edge in %d iterations, current guess %g+/-%g.",
                    find_enterleave_max_iter, *t, fabs(*t - t_old)/2);
        }

        t_old = *t; 

        if (GSL_SUCCESS != gsl_root_fdfsolver_iterate(solver))
            panic("gsl_root_fdfsolver failed");
        
        *t = gsl_root_fdfsolver_root(solver);
    
    } while (GSL_CONTINUE == gsl_root_test_delta(t_old, *t, find_enterleave_epsabs, find_enterleave_epsrel));


    // find rstar via interpolation
    *rstar = interpolate_position(*t, &intersect_low_pos).rstar;
                
    if ((*t-t0)*find_enterleave_epsrel + find_enterleave_epsabs < 0 ||  
            (t1-*t)*find_enterleave_epsrel + find_enterleave_epsabs < 0)
        panic("Intersect failed. Time is not in range %g <= %g <= %g", t0, *t, t1);

    // restore error handler
    gsl_set_error_handler(old_handler);
}
//
// void Particle::findEnterleave
// ---------------------------------------------------------------------------------------
// PURPOSE: calculates the points where the particle enters and leaves the cell
// INPUTS:
//  x0   : centre rstar coordinate of the cell
// RESULTS:
//  rstar1, t1   : point where the particle enters the cell (t1 is the full time)
//  rstar2, t2   : point where the particle leaves the cell (t2 is the full time)
// ---------------------------------------------------------------------------------------    
void Particle::findEnterleave(double *t1, double *rstar1, double *t2, double *rstar2, const double x0)
{
    // find t1 and rstar1
    if (past.rstar < x0) // entering to the left
        intersect(past.t, present.t, x0, -1., t1, rstar1);
    else // entering to the right
        intersect(past.t, present.t, x0, +1., t1, rstar1);

    // find t2 and rstar2
    if (future.rstar < x0) // leaving to the left
        intersect(present.t, future.t, x0-psimParams.Delta_t, +1., t2, rstar2);
    else // leaving to the right
        intersect(present.t, future.t, x0+psimParams.Delta_t, -1., t2, rstar2);
}        

// void Particle::findAreas(
//      const double x0, const double r1, const double t1, 
//      const double r2, const double t2,  
//      double &A1, double &A2, double &A3, double &A4)
//---------------------------------------------------------------------------------------------
// PURPOSE : this function calculates the areas that the passage of the particle cuts 
//  out of the diamond centered at (x0,t0)
// INPUTS:
//  x0       : centre coordinate of the cell
//  r1, t1   : point where the particle enters the cell (t1 is the difference to x_t*Delta_t)
//  r2, t2   : point where the particle leaves the cell (t1 is the difference to x_t*Delta_t)
// RESULTS:
//  A1   : area of the upper triangle
//  A2   : area of the left pentangle or hexangle
//  A3   : area of the lower triangle
//  A4   : area of the rigt pentangle or hexangle
// GLOBALS:
//---------------------------------------------------------------------------------------------
void Particle::findAreas( const double x0, 
        const double t1, const double rstar1, const double t2, const double rstar2,
        double *A1, double *A2, double *A3, double *A4)                
{
    const double h  = psimParams.Delta_t;    // NOTE: this is the chopped Delta_t
    const double t0 = present.t;

    // A3 and A2 are always triangles
    *A2 = 2*area_triangle(x0,t0-h, rstar1,t1, rstar1-2*(rstar1-x0),t1);
    *A3 = 2*area_triangle(x0,t0+h, rstar2,t2, rstar2-2*(rstar2-x0),t2);

    // A1 and A4 consist of upper and lower halves which might
    // be triangles or quadrangles
    if (past.rstar < x0) // A1_lower is a triangle, A4_lower is a quadrangle
    {
        // rstar1,t1 is to left
        *A1 = 2*area_triangle(rstar1,t1, present.rstar,t0, x0-h,t0);
        *A4 = 2*area_quadrangle(rstar1,t1, rstar1-2*(rstar1-x0),t1, x0+h,t0, present.rstar,t0);
    }
    else		  // A4_lower is a triangle, A1_lower is a quadrangle
    {
        // rstar1,t1 is right
        *A1 = 2*area_quadrangle(rstar1,t1, rstar1-2*(rstar1-x0),t1, x0-h,t0, present.rstar,t0);
        *A4 = 2*area_triangle(rstar1,t1, present.rstar,t0, x0+h,t0);
    }

    if (future.rstar < x0) // A1_upper is a triangle, A4_upper is a quadrangle
    {
        // rstar2,t2 is to left
        *A1 += 2*area_triangle(rstar2,t2, present.rstar,t0, x0-h,t0);
        *A4 += 2*area_quadrangle(rstar2,t2, rstar2-2*(rstar2-x0),t2, x0+h,t0, present.rstar,t0);
    }
    else		  // A4_upper is a triangle, A1_upper is a quadrangle
    {
        // rstar2,t2 is right
        *A1 += 2*area_quadrangle(rstar2,t2, rstar2-2*(rstar2-x0),t2, x0-h,t0, present.rstar,t0);
        *A4 += 2*area_triangle(rstar2,t2, present.rstar,t0, x0+h,t0);
    }

    if (fabs(*A1+*A2+*A3+*A4-4*SQR(h)) > 1e-5*SQR(h))
        panic("Areas: A1 = %g, A2 = %g, A3 = %g, A4 = %g, sum = %g, target = %g, \n", 
                *A1, *A2, *A3, *A4, *A1+*A2+*A3+*A4, 4*SQR(h));
}                       

