// cmdline.h
// Header file describing the functions used to parse the command line
//
// Author: Roland Haas
// File created Jul 24, 2006
// Last Change: Fri May 23 15:00 2008
//


#ifndef _CMDLINE_H_
#define _CMDLINE_H_

#include <stdio.h>
#include <zlib.h>


// 
// very poor implementation of a command line parser
// 
#define TYPE_DOUBLE 0
#define TYPE_INT 1
#define TYPE_FILE 2
#define TYPE_STRING 3
#define TYPE_INPUTFILE 4
#define TYPE_MAX 5
#define TYPE_OPTIONAL (1<<15) // the argument is optional, do not complain if it is not found
#define TYPE_FLAG_MASK TYPE_OPTIONAL
struct cmdline
{
    int type;
    const char *name;
    union
    {
	    void *dummy; // I can only initialize the first member of a union, 
		// and a void * can by definiton point to anything
	    double *double_data;
	    int *int_data;
	    FILE **file_data;
	    gzFile *inputfile_data;
	    char **str_data;
    };
};

int parse_commandline(int argc, char **argv, struct cmdline cmdline[], int max_cmd_args);

#endif // _CMDLINE_H_
