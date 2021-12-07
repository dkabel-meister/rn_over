// header.h
// Header file describing the functions used to write and check headers
// of the output files
//
// Author: Roland Haas
// File created Jun 18, 2006
// Last Change: Fri Apr 18 08:46 2008
//


#ifndef _HEADER_H_
#define _HEADER_H_

#include <zlib.h>
#include <stdio.h>
#include <time.h>

// Header to write to files
#define HEADER_SIZE 1024
// end of header marker
#define HEADER_END "# END OF HEADER\n"

extern char VERSTRING[];

class Header
{
    public:
        // default constructor
        Header();

        // read and write a header
        void read(gzFile fh);
        void write(FILE *fh);

        // add data items to a header
        void initialize(int argc, char **argv);
        void add(const char *name, int val);
        void add(const char *name, double val);
        void add(const char *name, char *val);
        void add(const char *name, double *val);

        // query daa items
        void query(const char *name, int *val);
        void query(const char *name, double *val);
        void query(const char *name, char **val);
        void query(const char *name, double **val);

    private:
        // time we started
        time_t started;
        // command line arguments
        int argc;
        char **argv;
        enum type {type_dbl, type_int, type_chr, type_vec};
        const char *get_typename(type type);
        struct data
        {
            data *next;
            const char *name;
            Header::type type;
            union {
                double dbl_val;
                int int_val;
                char *chr_val;
                double *vec_val;
            };
        } *data;
};

//extern class Header header;

#endif //_HEADER_H_
