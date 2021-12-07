// vim: sw=4 softtabstop=4 smartindent smarttab nowrap
// hermite.h
// Header file describing the function performing a hermite interpolation 
//
// Author: Roland Haas
// File created Sep 13, 2006
// Last Change: Sat Feb 10 11:41 2007
//


#ifndef _HERMITE_H_
#define _HERMITE_H_

extern double hermite_interpolate(double t0, double t1, double x0, double x1, double v0, double v1, double t);
extern double hermite_interpolate_df(double t0, double t1, double x0, double x1, double v0, double v1, double t);

#endif //_HERMITE_H_
