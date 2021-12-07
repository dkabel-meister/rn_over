#define DUMP 0

#include <unistd.h>
#include <cstdio>
#include <cmath>
#include <zlib.h>

#include "complex.h"
#include "cmdline.h"
#include "field.h"
#include "geometry.h"
#include "header.h"
#include "helper.h"
#include "mathhelper.h"
#include "particle.h"
#include "r_from_rstar.h"
#include "simdata.h"

static void defaultSetInputParams(SimulationParameter_t &fp);

class Field *field;
class Particle *particle;

int main(int argc, char **argv)
{
    SimulationParameter_t field_params;

    defaultSetInputParams(field_params);

    struct cmdline cmdline[] = 
    {
        {TYPE_INT,                    "--el",                   {&field_params.el}},
        {TYPE_INT,                    "--chop",                 {&field_params.chop}},
        {TYPE_DOUBLE,                 "--dt",                   {&field_params.Delta_t}},
        {TYPE_DOUBLE,                 "--tmax",                 {&field_params.t_max}},
        {TYPE_DOUBLE,                 "--r0",                   {&field_params.r_init}},
        {TYPE_DOUBLE,                 "--drdt0",                {&field_params.drdt_init}},
        {TYPE_OPTIONAL | TYPE_INT,    "--elmax",                {&field_params.el_max}},
        {TYPE_OPTIONAL | TYPE_INT,    "--quiet",                {&field_params.quiet}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--q",                    {&field_params.q}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--m",                    {&field_params.m}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--tmin",                 {&field_params.t_min}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--t_switch",             {&field_params.t_switch}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--rstar_inner_boundary", {&field_params.rstar_inner_boundary}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--rstar_outer_boundary", {&field_params.rstar_outer_boundary}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--observer",             {&field_params.observer}},
        {TYPE_OPTIONAL | TYPE_DOUBLE, "--Q",                    {&field_params.Q}},
        {TYPE_OPTIONAL | TYPE_FILE,   "--field_fh",             {&field_params.field_fh}},             
        {TYPE_OPTIONAL | TYPE_FILE,   "--flux_fh",              {&field_params.flux_fh}},             
        {TYPE_OPTIONAL | TYPE_FILE,   "--multipole_fh",         {&field_params.multipole_fh}},             
        {TYPE_OPTIONAL | TYPE_FILE,   "--particle_fh",          {&field_params.particle_fh}}             
    };

    copy_args(argc, argv);
    parse_commandline(argc, argv, cmdline, MAX_CMD_ARGS);

    field = new Field(field_params); // also sets pointer to Particle :) 
    
    // be nice to the CPU
    if (field_params.quiet)
        nice(10);

    if (field_params.t_min != 0)
        panic("t_min must be 0 (for now).");

    // initialize the particle
    //
    field->getParticle()->evolveParticle();

    // initialize the grid vars and the field arrays
    //
    field->initializeField();

    // write useful information to the header
    //
    field->writeHeader(argc, argv);
    field->getParticle()->writeHeader(argc, argv);

    // evolve the field 
    // This also does the interpolating of the raw particle data to the grid points 
    //
    field->evolveField();
}

//
// set default values
//
static void defaultSetInputParams(SimulationParameter_t &fp)
{
    fp.x_t = 0;

    // Set the optional arguments:
    // ---------------------
    // integers
    fp.quiet = 1;
    fp.el_max = fp.el;
    // doubles
    fp.t_min                = 0.0;
    fp.t_switch             = 0.0;
    fp.Q                    = 0.0;
    fp.q                    = 1.0;
    fp.m                    = 1.0;
    fp.rstar_inner_boundary = NAN;
    fp.rstar_outer_boundary = NAN;
    fp.observer             = NAN;
    // FILE pointers
    fp.field_fh     = NULL;
    fp.flux_fh      = NULL;
    fp.multipole_fh = NULL;
    fp.particle_fh  = NULL;
}
