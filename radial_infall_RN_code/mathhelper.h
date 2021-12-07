// mathhelper.h
// Header file describing the functions dealign in various mathematical aspects
// missing from the standard library
//
// Authors: Roland Haas and Peter Zimmerman
// File created Jun 18, 2006
//


#ifndef _MATHHELPER_H_
#define _MATHHELPER_H_

#include <gsl/gsl_math.h>
#include "complex.h"

// helpfull macros...
template<typename T> static inline T SQR(T a) {return a*a;}
template<typename T> static inline T CUBE(T a) {return a*a*a;}
template<typename T> static inline T QUAD(T a) {return SQR(a)*SQR(a);}
template<typename T, typename S> static inline T MAX(T x, S y) {return (x > y ? x : y);}
template<typename T, typename S> static inline T MIN(T x, S y) {return (x < y ? x : y);}

// sign bit (incl. that of zero)
static inline int sign(int a) {return  a == 0 ? 0 : a < 0 ? -1 : 1;}

// faster version of pow(double , double) for integer powers
#define pow_2(x) gsl_pow_2(x)
#define pow_3(x) gsl_pow_3(x)
#define pow_4(x) gsl_pow_4(x)
#define pow_5(x) gsl_pow_5(x)
#define pow_6(x) gsl_pow_6(x)
#define pow_7(x) gsl_pow_7(x)
#define pow_8(x) gsl_pow_8(x)
#define pow_9(x) gsl_pow_9(x)

#undef isnan
#define isnan(x) __isnan(x)

//
// various functions to do mathematical calculations
// 
extern int sq_roots(double a, double b, double c, double sols[2]);
extern double gaussian(double x, double x0, double w);
#endif //_MATHHELPER_H_
