// simdata.h
// Header file describing the spacetime and discretisation we use
//
// Authors: Roland Haas and Peter Zimmerman
// File created Jul 24, 2006


#ifndef _SIMDATA_H_
#define _SIMDATA_H_

typedef struct 
{
        int el;
        int el_max;
        int chop;
        int quiet;
        int x_t;

        double q;
        double Q;
        double m;
        double Delta_t;
        double t_min;
        double t_max;
        double t_switch;
        double r_init;
        double drdt_init;
        double rstar_inner_boundary; 
        double rstar_outer_boundary;
        double observer;

        // files
        FILE *particle_fh;
        FILE *field_fh;
        FILE *multipole_fh;
        FILE *flux_fh;
} SimulationParameter_t;

#endif // _SIMDATA_H_
