// helper.h
// Header file describing the functions used prnt error messages and other 
// helpfull stuff
//
// Author: Roland Haas
// File created Jun 18, 2006
// Last Change: Mon May 26 13:45 2008
//


#ifndef _HELPER_H_
#define _HELPER_H_

#include <cstdio>

// log a debug message
#ifdef DEBUG
#   define dbg_msg(msg, ...) fprintf(stderr, "%8s %04d " msg "\n", __FILE__, __LINE__, __VA_ARGS__)
#   define dbg_msgonce(msg, ...) do { static int _dbg_msgonce ; if(!(_dbg_msgonce)) { _dbg_msgonce = 1; fprintf(stderr, "%8s %04d " msg "\n", __FILE__, __LINE__, __VA_ARGS__) ; } } while(0)
#   define dbg_showvar(var) fprintf(stderr, "%8s %04d: " #var " = %.23g\n", __FILE__, __LINE__, (double)(var))
#   define dbg_ping fprintf(stderr, "%8s %04d\n", __FILE__, __LINE__)
#   define dbg_showonce(var) do { static int _dbg_showonce ; if(!(_dbg_showonce)) { _dbg_showonce = 1; dbg_showvar(var) ; } } while(0)
#   define dbg_marker {}
#else
#   define dbg_msg(msg,...) panic("Leftover debug statement in %s line %d.\n", __FILE__, __LINE__);
#   define dbg_msgonce(var) panic("Leftover debug statement in %s line %d.\n", __FILE__, __LINE__);
#   define dbg_showvar(var) panic("Leftover debug statement in %s line %d.\n", __FILE__, __LINE__);
#   define dbg_ping panic("Leftover debug statement in %s line %d.\n", __FILE__, __LINE__);
#   define dbg_showonce(var) panic("Leftover debug statement in %s line %d.\n", __FILE__, __LINE__);
#   define dbg_marker panic("Leftover debug statement in %s line %d.\n", __FILE__, __LINE__);
#endif

// helpfull macros...
#define DIM(v) ((int)(sizeof(v)/sizeof(v[0])))

// make copies of the command line arguments
extern void copy_args(int argc, char **argv);

__attribute__ ((__noreturn__))
extern void panic(const char *message, ...);

extern const char *gsl_error_header;
__attribute__ ((__noreturn__))
extern void gsl_panic(const char *reason, const char *file, int line, int gsl_errno);

extern void warn(const char *message, ...);
#endif //_HELPER_H_
