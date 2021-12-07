#ifndef _INTERSECT_H_
#define _INTERSECT_H_

#include "particle.h"


//
// structure to hold parameters used to intersect the particle's trajectory
// and the edges of the cell it passes through
// 
struct positionParams 
{ 
    double t0;	    // time from which to start (past or present)
    double rstar0;  // vertex of the edge on the past/current timeslice
    double slope;   // slope of the edge, +1 for LR and UL edge, -1 for LL and UR edges

    struct Particle::internal_pos low, high;
};

#endif //  _INTERSECT_H_
