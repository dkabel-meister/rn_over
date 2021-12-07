// header.cc 
// funcion to create and check headers of data fields
//
// File created on Tue May 08, 2007
// Last Change: Mon Sep 07 18:03 2009
// Author: Roland Haas
//

#include <zlib.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>

#include "header.h"
#include "helper.h"

// constructor
Header::Header()
{
    // set all fields to invalid values so that initialize can do something...
    started = -1;
    argc = -1;
    argv = NULL;
    data = NULL;
}

// read(gzFile fh)
// read in a header from the filehandle *fh
// inputs:
// fh	    - where to read from
// results:
// none
void Header::read(gzFile fh)
{
    char buf[1024]; // buffers one line of input
    struct tm tm;   // time of call
    char *args, *name, *val;	// pointers to various subparts of buf
    int bytes, len; // used to verify that all bytes have been used in type conversion
    int tmp_int;       // holds an integer
    double tmp_dbl;    // holds a double
    double tmp_vec[4]; // holds a 4-vector

    errno = 0;

    // read and check first lines
    if(NULL == gzgets(fh, buf, DIM(buf)))
	    panic("Could not read from header: %s", strerror(errno));
    if(0 != strncmp("# Simulation started on ", buf, 24) || NULL == strptime(&buf[24], "%a %b %d %H:%M:%S %Y\n", &tm))
	    panic("Malformed header line:\n%s", buf);

    if(NULL == gzgets(fh, buf, DIM(buf)))
	    panic("Could not read from header: %s", strerror(errno));
    if(0 != strncmp("# $VER: ", buf, 8))
	    panic("Malformed header line:\n%s", buf);
    buf[strlen(buf)-1] = '\0'; // get rid of the '\n'
    if(0 != strcmp(&buf[8], VERSTRING))
	    warn("Version mismatch between code and data.");

    if(NULL == gzgets(fh, buf, DIM(buf)))
	    panic("Could not read from header: %s", strerror(errno));
    if(0 != strncmp("# cmd: ", buf, 7))
	    panic("Malformed header line:\n%s", buf);

    // save away command line used
    buf[strlen(buf)-1] = '\0'; // get rid of the '\n'
    add("cmd", &buf[7]);

    // fill in data fields
    // 
    started = mktime(&tm); 

    // fill in the tagged data
    data = NULL;
    while(NULL != gzgets(fh, buf, DIM(buf)))
    {
	    if(0 == strcmp(HEADER_END, buf))
	        return;

	    if(0 != strncmp("# val: ", buf, 7))
	        panic("Malformed header line:\n%s", buf);

	    args = strtok(&buf[7], " \n");
	    while(args != NULL)
	    {
	        name = args;
	        val = strchr(args, '=');
	        if(val == NULL)
		        panic("Malformed name=value pair: %s", args);
	        *(val++) = '\0';

	        // decide on the type of val and store it accordingly
	        len = (int)strlen(val);
	        if(1 == sscanf(val, "%d%n", &tmp_int, &bytes) && bytes == len)
		        add(strdup(name), tmp_int);
	        else if(1 == sscanf(val, "%lf%n", &tmp_dbl, &bytes) && bytes == len)
		        add(strdup(name), tmp_dbl);
	        else if(1 == sscanf(val, "(%lf,%lf,%lf,%lf)%n", &tmp_vec[0], &tmp_vec[1], &tmp_vec[2], &tmp_vec[3], &bytes) && bytes == len)
		        add(strdup(name), tmp_vec);
	        else
		        add(strdup(name), strdup(val));

	        args = strtok(NULL, " \n");
	    }
    }
}

// write(FILE *fh)
// write out a header to the filehandle *fh
// inputs:
// fh	    - where to write to
// results:
// none
static int format(FILE *fh, int &len, const char *fmt, ...)
{
    int add;
    va_list data;

    // find out if we have to start a new line
    va_start(data, fmt);
    add = vsnprintf(NULL, 0, fmt, data);
    va_end(data);
    if(len + add >= 80)
    {
	    fputs("\n# val:", fh);
	    len = 6;
    }

    // output data and update length counter
    va_start(data, fmt);
    vfprintf(fh, fmt, data);
    va_end(data);
    len += add;

    fflush(fh);

    return add;
}
void Header::write(FILE *fh)
{
    int len;

    // write down the time the program was started and its version
    fprintf(fh, "# Simulation started on %s", asctime(localtime(&started)));
    fprintf(fh, "# %s\n", VERSTRING);
    // write down the command line
    fputs("# cmd:", fh);
    for(int i = 0 ; i < argc ; i++)
    {
	    fputc(' ', fh);fputs(argv[i], fh);
    }
    // newline will be printed further below

    // add all the name value pairs added via header::add
    len = 80; // makes sure that a new line is started
    for(struct data *data_ptr = data ; data_ptr != NULL ; data_ptr = data_ptr->next)
    {
	    switch(data_ptr->type)
	    {
	        case type_dbl:
		        format(fh, len, " %s=%e", data_ptr->name, data_ptr->dbl_val);
		        break;
	        case type_int:
		        format(fh, len, " %s=%d", data_ptr->name, data_ptr->int_val);
		        break;
	        case type_chr:
		        format(fh, len, " %s=%s", data_ptr->name, data_ptr->chr_val);
		        break;
	        case type_vec:
		        format(fh, len, " %s=(%e,%e,%e,%e)", data_ptr->name, data_ptr->vec_val[0], data_ptr->vec_val[1], data_ptr->vec_val[2], data_ptr->vec_val[3]);
		        break;
	        default:
		        panic("Unknown type used in %s line %d.", __LINE__, __LINE__);
		        break;
	    }
    }
    fputs("\n" HEADER_END, fh); // end line and write end of header marker
}

// add data items to a header
void Header::initialize(int argcount, char **argvector)
{
    started = time(NULL);
    argc = argcount;
    argv = argvector;
    data = NULL;
}

void Header::add(const char *name, int val)
{
    struct data *curr;

    // allocate a new one and link it into the list
    curr = new struct data;
    curr->next = data;
    data = curr;

    // enter data
    curr->type = type_int;
    curr->name = name;
    curr->int_val = val;
}

void Header::add(const char *name, double val)
{
    struct data *curr;

    // allocate a new one and link it into the list
    curr = new struct data;
    curr->next = data;
    data = curr;

    // enter data
    curr->type = type_dbl;
    curr->name = name;
    curr->dbl_val = val;
}

void Header::add(const char *name, char *val)
{
    struct data *curr;

    // allocate a new one and link it into the list
    curr = new struct data;
    curr->next = data;
    data = curr;

    // enter data
    curr->type = type_chr;
    curr->name = name;
    curr->chr_val = strdup(val);
}

void Header::add(const char *name, double val[4])
{
    struct data *curr;

    // allocate a new one and link it into the list
    curr = new struct data;
    curr->next = data;
    data = curr;

    // enter data
    curr->type = type_vec;
    curr->name = name;
    curr->vec_val = new double[4];
    for(int i = 0 ; i < 4 ; i++)
	    curr->vec_val[i] = val[i];
}

// convert type into its name
const char *Header::get_typename(type tp)
{
    static const char *typenames[] = {"double", "int", "string", "4-vector"};
    return typenames[tp];
}

// query data items
void Header::query(const char *name, int *val)
{
    struct data *curr;
    for(curr = data ; curr!= NULL ; curr = curr->next)
    {
	    if(0 == strcmp(curr->name, name))
	        break;
    }

    if(curr == NULL)
	    panic("Cannot find data element named %s.", name);

    if(curr->type != type_int)
	    panic("Type mismatch for element %s in header; %s requested but value is %s.", name, get_typename(type_int), get_typename(curr->type));

    *val = curr->int_val;
}

void Header::query(const char *name, double *val)
{
    struct data *curr;
    for(curr = data ; curr!= NULL ; curr = curr->next)
    {
	    if(0 == strcmp(curr->name, name))
	        break;
    }

    if(curr == NULL)
	    panic("Cannot find data element named %s.", name);

    if(curr->type != type_dbl)
	    panic("Type mismatch for element %s in header; %s requested but value is %s.", name, get_typename(type_dbl), get_typename(curr->type));

    *val = curr->dbl_val;
}

void Header::query(const char *name, char **val)
{
    struct data *curr;
    for(curr = data ; curr!= NULL ; curr = curr->next)
    {
	    if(0 == strcmp(curr->name, name))
	        break;
    }

    if(curr == NULL)
	    panic("Cannot find data element named %s.", name);

    if(curr->type != type_chr)
	    panic("Type mismatch for element %s in header; %s requested but value is %s.", name, get_typename(type_chr), get_typename(curr->type));

    *val = curr->chr_val;
}

void Header::query(const char *name, double **val)
{
    struct data *curr;
    for(curr = data ; curr!= NULL ; curr = curr->next)
    {
	    if(0 == strcmp(curr->name, name))
	        break;
    }

    if(curr == NULL)
	    panic("Cannot find data element named %s.", name);

    if(curr->type != type_vec)
	    panic("Type mismatch for element %s in header; %s requested but value is %s.", name, get_typename(type_vec), get_typename(curr->type));

    *val = curr->vec_val;
}
