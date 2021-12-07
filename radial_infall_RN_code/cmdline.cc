// cmdline.c contains a very simplistic implementation of a command line parser
//
// File created on Mon Jul 24, 2006
// Last Change: Fri May 23 15:00 2008
// Author: Roland Haas
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "helper.h"
#include "cmdline.h"
#include "header.h"

//
// poor man's implementation of a command line parser
// inputs:
// argc, argv	: the command line argument count and array
// cmdline	: describes how the arguments are to be interpreted
// max_cmd_args	: DIM(cmdline)
// globals:
// none
// results:
// fills in the data pointers in cmdline[] and returns the index of the argument
// '--' if such an argument exists or 0 otherwise
// 
int parse_commandline(int argc, char **argv, struct cmdline cmdline[], int max_cmd_args)
{
    int *found_required = (int *)__builtin_alloca(sizeof(int)*max_cmd_args);
    static const char *typenames[TYPE_MAX] = {"double", "int",  "filename", "string"};

    // print help text if asked to do so
    if(argc == 2 && strcmp(argv[1], "--help") == 0)
    {
	    fprintf(stderr, "%s\n", VERSTRING);
	    fprintf(stderr, "call syntax: %s ", argv[0]);
	    for(int c = 0 ; c < max_cmd_args ; c++)
	    {
	        if(cmdline[c].type & TYPE_OPTIONAL)		
		        fprintf(stderr, "[%s {%s}] ", cmdline[c].name, typenames[cmdline[c].type & ~TYPE_FLAG_MASK]);
	        else		
		        fprintf(stderr, "%s {%s} ", cmdline[c].name, typenames[cmdline[c].type & ~TYPE_FLAG_MASK]);
	    }
	    fputc('\n', stderr);
	    fprintf(stderr, "By default all data is discarded, any unassigned filenames will point to /dev/null.\n");
	    exit(0);
    }

    if(0 == (argc & 1)) // an odd number of extra arguments is always invalid
	    panic("Malformed command line, try %s --help\n", argv[0]);

    // mark that we have not seen any required arguments so far :-)
    for(int c = 0 ; c < max_cmd_args ; c++)
	    found_required[c] = 0;

    // act on each argument
    for(int a = 1 ; a < argc ; a+=2)
    {
	    int c;

	    if(0 == strcmp(argv[a], "--"))
	        return (a+1);

	    for(c = 0 ; c < max_cmd_args ; c++)
	    {
	        if(0 == strcmp(argv[a], cmdline[c].name)) // found a known argument
	        {
		        if(a + 1 == argc) // we miss the data for this argument
		        {
		            fprintf(stderr, "Missing data for last argument '%s'\n", argv[a]);
		            exit(1);
		        }

		        // mark that we have seen this argument
		        found_required[c] = 1;

		        switch(cmdline[c].type & ~TYPE_FLAG_MASK)
		        {
		            case TYPE_DOUBLE:
			            *(cmdline[c].double_data) = atof(argv[a + 1]);
			            break;
		            case TYPE_INT:
			            *(cmdline[c].int_data) = atoi(argv[a + 1]);
			            break;
		            case TYPE_FILE:
			            if(0 == strcmp("-", argv[a + 1]))
			            {
			                *(cmdline[c].file_data) = stdout;
			            }
			            else
			            {
			                *(cmdline[c].file_data) = fopen(argv[a + 1], "w+");
			                if(NULL == *(cmdline[c].file_data))
			                {
				                panic("Could not open file '%s' for '%s'\n", argv[a+1], argv[a]);
			                }
			            }
			            break;
		            case TYPE_INPUTFILE:
			            if(0 == strcmp("-", argv[a + 1]))
			            {
			                *(cmdline[c].inputfile_data) = gzdopen(fileno(stdin), "r");
			            }
			            else
			            {
			                *(cmdline[c].inputfile_data) = gzopen(argv[a + 1], "r");
			                if(NULL == *(cmdline[c].file_data))
			                {
				                panic("Could not open file '%s' for '%s'\n", argv[a+1], argv[a]);
			                }
			            }
			            break;
		            case TYPE_STRING:
			            *(cmdline[c].str_data) = strdup(argv[a + 1]);
			            break;
		            default:
			            panic("Unknown argument type in parsing the command line");
			            break;
		        }

		        break; // break out of the inner 'for' loop since we have found out how to parse this argument
	        }
	    }
	    if(max_cmd_args == c) // we could not find a match for this argument
	        panic("Unknown command line argument '%s', try %s --help\n", argv[a], argv[0]);
    }

    // check that we have seen all the required arguments
    for(int c = 0 ; c < max_cmd_args ; c++)
    {
	    if(!(cmdline[c].type & TYPE_OPTIONAL) && 0 == found_required[c])
	        panic("Required argument %s is missing.", cmdline[c].name);
    }

    return 0;
}
