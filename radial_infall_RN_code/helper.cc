// helper.cpp contains various helper functions used throughout the simulation
//
// File created on Jan 25, 2006
// Last Change: Thu Jan 31 14:22 2008
// Author: Roland Haas
//

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <gsl/gsl_errno.h>
#include <unistd.h>

#include "helper.h"

// copy of the command line arguments for later use
static char **cmd_args;
static int cmd_argcount;

// header printed out in front of gsl error messages, use to tell me what I did
const char *gsl_error_header = "<unknown>";

void copy_args(int argc, char **argv)
{
	cmd_args = argv;
	cmd_argcount = argc;
}

//
// void panic(const char *message, ...)
// this function should be called whenever a fatal inconsistency in the simulation
// is detected. It will quit the program after displaying the message and dumping
// all state information accessible to it.
// The arguments are printf like.
//
void panic(const char *message, ...)
{
	va_list ap;

	fputs("PANIC:", stderr);
	for(int i = 0 ; i < cmd_argcount ; i++)
	    fprintf(stderr, " %s", cmd_args[i]);
	fputc('\n', stderr);
	
	va_start(ap, message);
	vfprintf(stderr, message, ap);
	fputc('\n', stderr);
	// TODO:
	// dump data here, maybe use abort()
//	dbg_marker;
//	sleep(10);
#if 0
	{ 
		char ddd_buf[128];
		snprintf(ddd_buf, DIM(ddd_buf), "ddd %s %d", cmd_args[0], getpid());
		system(ddd_buf);
	}
#endif

	abort();
}

// void gsl_panic(const char *reason, const char *file, int line, int gsl_errno)
// my error handler for the gsl. Reproduces the gsl standard output and prepends
// it with "gsl_error_header: "
// inputs:
// reason	: gsl specific reason
// file		: gsl source file triggering the error
// line		: line in file
// errno	: one of the constants in gsl_errno
// globals:
// none
// results:
// quits the program
void gsl_panic(const char *reason, const char *file, int line, int gsl_errno)
{
	fprintf(stderr, "%s: %s:%d: ERROR: %s (%s)\n", gsl_error_header, file, line, reason, gsl_strerror(gsl_errno));
	// TODO:
	// dump data here, maybe use abort()
	panic("gsl error occured");
}

//
// void warn(const char *message, ...)
// this function should be called whenever a non-fatal inconsistency in the simulation
// is detected. It will display the message 
// The arguments are printf like.
//
void warn(const char *message, ...)
{
	va_list ap;

	fputs("WARNING:", stderr);
	va_start(ap, message);
	vfprintf(stderr, message, ap);
	fputc('\n', stderr);
}
