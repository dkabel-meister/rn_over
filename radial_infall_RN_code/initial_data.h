// initial_data.h
// Header file describing the functions describing the inital data of the 
// field problem
//
// Authors: Roland Haas & Peter Zimmerman
// File created Jun 18, 2006
//

#ifndef _INITIAL_DATA_H_
#define _INITIAL_DATA_H_

#include "field.h"

extern double gaussian_peak;
extern double gaussian_width;


// gaussian initial data
double psi0_gaussian(const double rstar);

class InitialData
{
    public:
        InitialData();

        static double legendre_Pl(const int el, const double x);
        static double legendre_Q0(const double x);
        static double legendre_Q1(const double x);
        static double legendre_Ql(const int el, const double x);
};


#endif  //_INITIAL_DATA_H_
