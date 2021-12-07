// particle.h
// Header file describing the functions dealign with the simulation of 
// the particle
//
// Author: Roland Haas and Peter Zimmerman
// File Created: Sat Sept 17, 2010 (Schw)
// File Copied for RN modifications: Thu Jan 6, 2011 
///
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <cstdio>
#include <cmath>
#include <gsl/gsl_roots.h>
#include "geometry.h"
#include "simdata.h"
#include "mathhelper.h"

extern class Particle *p;

struct gsl_pos
{
    double t, r, rstar, rdot;
};

class Particle : public Geometry
{
    public:
        Particle(SimulationParameter_t ps);  // default constructor
        ~Particle(); // default destructor

        // initializes particle data
        void evolveParticle();

        // calculate the paritcle's position at a given time using a Hermite interpolation
        gsl_pos interpolate_position(const double t, int *low_pos);

        // 
        // functions defined in intersect.cc
        //
        void intersect(const double t0, const double t1, 
                const double rstar0, const double slope, double *t, double *rstar); 
        void findEnterleave(double *t1, double *rstar1, double *t2, double *rstar2, const double x0); 
        void findAreas(const double x0, double r1, double t1, double r2, double t2, 
                double *A1, double *A2, double *A3, double *A4); 

        inline double get_E(void) { return E; };
        inline double get_rEqrstar(void) { return r_switch; };

        void writeHeader(int argc, char **argv);
        double rstar_from_r(const double r);
        
        // step the particle forward in time and shuffle positions
        void step(void);

        // dump the data to file *fh
        void dump(void);

        // positions updated by step particle
        // treat as read only
        gsl_pos future, past, present, ancient;

        struct internal_pos // particle positions (Called hermite data in odd)
        {
            double t;
            double rstar;
            double vrstar;
        };                                     
        
    private:
        SimulationParameter_t psimParams;

        bool header_written;

        int x_t;           // time step counter
        int N_positions; 
        
        double E;          // particle energy 
        double m;          // particle mass   
        double Q;          // BH charge 
        double q;          // particle charge
        double r_switch;

        FILE *fh;

        // precomputed positions (these the the values of (t,r^*, dr^*/dt before interpolation)
        struct internal_pos *internal_pos;

        friend int geodesic_rhs(double t, const double x[], double v[], void *params); // NOTE: friend allow access outside Particle class 

        typedef int bracket_cmp(struct Particle::internal_pos *pos, void *params);

        friend bracket_cmp time_cmp;
        friend bracket_cmp lightcone_cmp;

        void get_bracket(bracket_cmp (*cmp), void *param, int *low_pos);

        // intersections with light rays
        gsl_root_fdfsolver *solver; // control structure for gsl

        int intersect_low_pos;

};

inline double Particle::rstar_from_r(const double r)
{
        double div = (M+sqrt(M*M-Q*Q) - (M-sqrt(M*M-Q*Q)));
        double rs = r + SQR(M+sqrt(M*M-Q*Q))/div * log(r- (M+sqrt(M*M-Q*Q))) 
            - SQR(M-sqrt(M*M-Q*Q))/div * log(r-(M-sqrt(M*M-Q*Q)));
        return rs;

}        
#endif // _PARTICLE_H_
