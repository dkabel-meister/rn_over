// mathhalper.cc
// contains functions to do some math which are not included in the 
// C/C++ libray
//
// File created on Feb 24, 2006
// Author: Roland Haas and Peter Zimmerman

#include <cmath>
#include <cfloat>
#include <vector>

#include "mathhelper.h"
#include "helper.h"


//
// double gaussian(double x, double x0, double w)
// calculates a Gaussian: exp(-(x-x0)^2/(2*w^2))
//
// inputs:
// x	: where to evaluate
// x0	: peak position
// w	: width
// results:
// vaue of the Gaussian
// globals:
// none
// 
double gaussian(double x, double x0, double w)
{
    return 1/(w*sqrt(2*M_PI))*exp(-SQR(x-x0)/(2*SQR(w)));
}

//
// int sq_roots(double a, double b, double c, double sols[2])
// calcualtes the root(s) of 
// a x^2 + bx + c = 0
// and returns them in sols[], the return value of the function
// itself is the number of solutions found (-1,0,1,2)
// -1 signals that a=b=c=0
//
// inputs:
// a,b,c    : coefficients of the polynomial (see above)
// sols	    : array in which to save the solutions
// results:
// numer of solutions found
// globals:
// none
//
int sq_roots(double a, double b, double c, double sols[2])
{
    double det;
    int nsols;

    if(a == 0) // degenerate equation, it is linear
    {
	if(b == 0) // even worse, it does not even depend on x
	{
	    if(c == 0) // ok, we have 0=0, everything will do
	    {
		nsols = -1;
	    }
	    else // not solvable
	    {
		nsols = 0;
	    }
	}
	else // proper linear equation
	{
	    sols[0] = -c/b;
	    nsols = 1;
	}
    }
    else    // proper quadratic equation
    {
	det = b*b-4*a*c;

	if(det < 0) // no solution
	{
	    nsols = 0;
	}
	else if(det == 0) // one solution
	{
	    sols[0] = -b/(2*a);
	    nsols = 1;
	}
	else // two solutions
	{
	    det = sqrt(det);
	    
	    if(b < 0) // avoid extinction
	    {
		sols[0] = (-b+det)/(2*a);
		sols[1] = -2*c/(b-det);
	    }
	    else
	    {
		sols[0] = -2*c/(b+det);
		sols[1] = (-b-det)/(2*a);
	    }
	    
	    nsols = 2;
	}
    }
    
    return nsols;
}

