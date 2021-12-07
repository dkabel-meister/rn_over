#ifndef _FIELD_H_
#define _FIELD_H_

#include "particle.h"
#include "simdata.h"
#include "initial_data.h"

#define MAX_CMD_ARGS DIM(cmdline)
#define NUM_REGS 5 // retarded, A, B, D and tail

#define N_SLICES 4


#define FUTURE   3
#define PRESENT  2
#define PAST     1  
#define ANCIENT  0

#define EVEN     0
#define ODD      1
#define cell_min 0

extern int low_pos;
//
// field values
//
extern double *psi[N_SLICES];
extern double **fields[]; 

extern class Field *field;

class Field 
{
    public:
        Field(SimulationParameter_t fs);
        ~Field();

        double rstar_from_cell(int isodd, int cell);
        int cell_from_rstar(bool isodd, double rstar);
        int  cell_from_r(int isodd, double r);

        double transitionSource(double t);

        void initializeField(void); 
        double psi_initial_data(const double, const double);
        void evolveField(void);

        void integrateSource(const double t1, const double t2, const double x0, double &source);
        
        void regularize(double retPsi, double *regField);

        void writeHeader(int argc, char **argv); 
        void dumpAtObserver(int t_counter);
        void dump(const int t_counter);

        // 
        // member functions for jumps (defined in extractor.cc)
        // 
        void getJumps(double *jumps, gsl_pos *pos);
        void extrapolateRight(const double del, const double* const vals, double *field0, double *field_rstar);
        void extrapolateLeft(const double del, const double* const vals, double *field0, double *field_rstar);
        void interpolate(const double* const field, const int isodd, 
                gsl_pos *pos, double *field0, double *field_rstar);

        Particle *getParticle() { return particle; }

    private:
        Particle *particle;
        SimulationParameter_t fsimParams;

        gsl_pos pos_tmin, pos_tmax;

        bool do_excision;
        int x_t_excision;

        int el;
        int el_max;
        int x_t;
        int N_cells;
        int cell_left, cell_right;
        int cell_inner_boundary, cell_outer_boundary;
        int observer_cell;
        
        double m;
        double q;
        double Q;
        double Delta_t, Delta_rstar; 
        double t_min, t_max;
        double t_switch;
        double rstar_min, rstar_max;

        // translation table from cell numbers to r values
        double *r_from_cell[2];

};

inline double Field::rstar_from_cell(int isodd, int cell)  
{ 
    return rstar_min + cell*2*Delta_t + (isodd ? Delta_t : 0); 
}

inline int Field::cell_from_rstar(bool isodd, double rstar) 
{ 
    if (isodd) 
        return (int)floor((rstar-rstar_min-Delta_t)/(2*Delta_t)); 
    else 
        return (int)floor((rstar-rstar_min)/(2*Delta_t));
}

inline int Field::cell_from_r(int isodd, double r) 
{
    int c;
    for (c=0;r_from_cell[isodd][c] <= r;c++) ; return c-1;
}              

#endif // _FIELD_H_
