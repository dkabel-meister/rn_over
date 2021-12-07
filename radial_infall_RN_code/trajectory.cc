// trajectory.cc : simulates the geodesic orbit of a particle around a
//  RN black hole
// File created on Sat Sept 18, 2010
// File copied for RN modifications: Fri Jan 7, 2011
// Authors: Roland Haas and Peter Zimmerman
///


// include all the header files we need
//
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <fenv.h>

#include "cmdline.h"
#include "header.h"
#include "helper.h"
#include "mathhelper.h"
#include "particle.h"
#include "r_from_rstar.h"
#include "trajectory.h"
#include "simdata.h"

static void defaultSetInputParams(SimulationParameter_t  &op);
// 
// objects that make up the simulation
//
class Particle *particle;   

int main(int argc, char **argv)
{   
    SimulationParameter_t orb_params;
    
    defaultSetInputParams(orb_params);

    struct cmdline cmdline[] = 
    {
        {TYPE_INT,                    "--chop",                   {&orb_params.chop}},
        {TYPE_DOUBLE,                 "--dt",                     {&orb_params.Delta_t}},
        {TYPE_DOUBLE,                 "--tmax",                   {&orb_params.t_max}},
        {TYPE_DOUBLE,                 "--r0",                     {&orb_params.r_init}},
        {TYPE_DOUBLE,                 "--drdt0",                  {&orb_params.drdt_init}},
        {TYPE_OPTIONAL | TYPE_INT,    "--el",                     {&orb_params.el}},
        {TYPE_OPTIONAL | TYPE_INT,    "--quiet",                  {&orb_params.quiet}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--tmin",                   {&orb_params.t_min}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--rstar_inner_boundary",   {&orb_params.rstar_inner_boundary}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--rstar_outer_boundary",   {&orb_params.rstar_outer_boundary}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--observer",               {&orb_params.observer}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--Q",                      {&orb_params.Q}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--q",                      {&orb_params.q}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--m",                      {&orb_params.m}},
        {TYPE_OPTIONAL | TYPE_FILE,   "--field_fh",               {&orb_params.field_fh}},             
        {TYPE_OPTIONAL | TYPE_FILE,   "--particle_fh",            {&orb_params.particle_fh}}             
    };                   

    // squirrel away the initial command line arguments for panic()
    //
    copy_args(argc, argv);    

    // parse the command line arguments
    //
    parse_commandline(argc, argv, cmdline, MAX_CMD_ARGS);

    // fill orbsim with data from the command line
    //
    particle = new Particle(orb_params);

    // run the particle evolution without interpolating to the slices
    //
    particle->evolveParticle();

    // write useful information to the header (we evolve first to get the energy)
    //
    particle->writeHeader(argc, argv);

    particle->dump(); // dump the initial phase-space point

    // run the simulation and fill in the particle data on the slices
    //
    while (particle->present.t < orb_params.t_max)
    {
        if (!orb_params.quiet && fabs(fmod(particle->present.t, orb_params.t_max/100)) < orb_params.Delta_t)
            fprintf(stderr,"%5.0f\r", particle->present.t);

        particle->step();
        particle->dump();
    }                             

    fputc('\n', stderr);

    return 0; // if we ever get this far, everything went ok, panic()
    // contains its own exit statement 
}

static void defaultSetInputParams(SimulationParameter_t &op)
{
    op.x_t = 0;

    // Optional arguments:
    // ---------------------
    // integers
    op.quiet = 1;  // default is shhh!
    // doubles
    op.q                    = 1.0;  // default to unit charge
    op.m                    = 1.0;  // default to unit mass
    op.Q                    = 0.0;  // default to Schw
    op.t_min                = 0.0;
    op.rstar_inner_boundary = NAN;
    op.rstar_outer_boundary = NAN;
    op.observer             = NAN;
    // FILE pointers
    op.field_fh     = NULL;
    op.flux_fh      = NULL;
    op.multipole_fh = NULL;
    op.particle_fh  = NULL;
}                               
