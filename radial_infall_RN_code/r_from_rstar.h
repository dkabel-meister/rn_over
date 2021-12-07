// r_from_rstar.h
// Header file describing the functions used to speed up the calculation of 
// r from rstar on grid points
//
// Authors: Roland Haas and Peter Zimmerman
// File created Jun 18, 2006
// Last modified: Mon Jan 10 2011
//


#ifndef _R_FROM_RSTAR_H_
#define _R_FROM_RSTAR_H_

#include <cmath>
#include "mathhelper.h"


inline double rstar_from_r_Schw(const double r) { return r+2*M*log(r/(2*M)-1); }

class Rstar 
{
    public:
        Rstar() {;};

        static double rstar_from_r_RN(const double r, const double Q);
        static double rstar_exp_finder_RN_f(const double r, const double rstar, const double Q);
        static double rstar_exp_finder_RN_df(const double r, const double rstar, const double Q);
        static double rstar_ln_finder_RN_f(const double r, const double rstar, const double Q);
        static double rstar_ln_finder_RN_df(const double r, const double Q);
        static double rPlus(const double Q); 
        static double rMinus(const double Q); 
        static double alphaPlus(const double Q);
        static double alphaMinus(const double Q);
        static double logArg(const double, const double Q);
};

double r_from_rstar(const double rstar, const double Q, const double r_switch);
double get_r_switch(const double Q);

#endif //_R_FROM_RSTAR_H_
