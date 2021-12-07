// vim: nowrap smarttab smartindent sw=4 softtabstop=4
// hermite.c
// contains the function used to perform a Hermite interpolation
// Author: Roland Haas
// File Created: Wed Sep 13, 2006
// Last Change: Sat Feb 10 11:41 2007
//

#include "hermite.h"

//
// hermite_interpolate
// calculates the hermite interpolation at position t
// inputs:
// t0	    : lower time
// y0,v0    : function value and first derivative at t0
// t1	    : upper time
// y1,v1    : function value and first derivative at t1
// t	    : we want to know y(t) here
// globals:
// none
// 
double hermite_interpolate(double t0, double t1, double x0, double x1, double v0, double v1, double t)
{
    // some helper objects
    double lambda;	// rescaled time so that t0==0, t1==1
    double xp0, xp1;	// derivatives with respect to lambda
    double dt;		// t1-t0 
    
    // rescale everything
    dt = t1 - t0;
    lambda = (t - t0)/dt;
    xp0 = v0 * dt;
    xp1 = v1 * dt;
    
    // formula found in hermite_interpolate.mws
    return x0+lambda*(xp0+lambda*((-2*xp0-3*x0+3*x1-xp1)+lambda*(2*x0+xp0-2*x1+xp1)));
}

//
// hermite_interpolate_df
// calculates the derivative of the hermite interpolation at position t
// inputs:
// t0	    : lower time
// y0,v0    : function value and first derivative at t0
// t1	    : upper time
// y1,v1    : function value and first derivative at t1
// t	    : we want to know y(t) here
// globals:
// none
// 
double hermite_interpolate_df(double t0, double t1, double x0, double x1, double v0, double v1, double t)
{
    // some helper objects
    double lambda;	// rescaled time so that t0==0, t1==1
    double xp0, xp1;	// derivatives with respect to lambda
    double dt;		// t1-t0 
    
    // rescale everything
    dt = t1 - t0;
    lambda = (t - t0)/dt;
    xp0 = v0 * dt;
    xp1 = v1 * dt;
    
    // formula found in hermite_interpolate.mws
    return (xp0+lambda*(2*(-2*xp0-3*x0+3*x1-xp1)+lambda*3*(2*x0+xp0-2*x1+xp1)))/dt;
}
