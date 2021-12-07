# Makefile for RN
# Authors: Roland Haas and Peter Zimmerman
# Created: Sept 19, 2010

UNAME := $(shell uname)

TRAJECTORY_SRCS := cmdline.cc header.cc helper.cc hermite.cc particle.cc \
	timestamp.cc mathhelper.cc r_from_rstar.cc  trajectory.cc intersect.cc \
	geometry.cc
TRAJECTORY_OBJS = $(patsubst %.cc,%.o, $(TRAJECTORY_SRCS))

FIELD_SRCS := cmdline.cc header.cc helper.cc hermite.cc  \
	timestamp.cc field.cc extractor.cc generator.cc \
	particle.cc intersect.cc geometry.cc initial_data.cc regularize.cc \
 	r_from_rstar.cc transition.cc
FIELD_OBJS = $(patsubst %.cc,%.o, $(FIELD_SRCS))


ALL_SRCS = $(TRAJECTORY_SRCS) $(FIELD_SRCS)
ALL_OBJS = $(patsubst %.cc,%.o,$(ALL_SRCS)) 

OPTIMIZE = -O3 
DEBUG    = -g -DDEBUG 
PROFILE  = -pg
INCFLAGS = -I./gsl/ -I. -I$(HOME)/local/incude


ifeq "$(UNAME)"  "Darwin"
# Mac 
OPENMP = 
AFLAGS = -mtune=core2
#LDFLAGS = -L. $(PROFILE) 
WERROR   =
HAVE_ISNAN  = 
else
# LINUX 
OPENMP = -fopenmp	
AFLAGS=-march=native -m64
#LDFLAGS = $(PROFILE) $(OPENMP) -L$(HOME)/local/lib -L./gsl
WERROR = -Werror
HAVE_ISNAN  = -g -DHAVE_IS_NAN 
endif


CXXFLAGS = -D_GNU_SROUCE \
	$(PROFILE) $(OPTIMIZE) $(ARCH) $(DEBUG) $(OPENMP) $(HAVE_ISNAN)\
	-Wno-missing-braces \
	-Wall -Wunused-variable \
	$(WERROR) \
	-Wno-deprecated \
	-Wno-format-extra-args \
       	-Wconversion -Wshadow \
       	-Wpointer-arith -Wcast-align \
       	-Wwrite-strings \
       	-fshort-enums -fno-common \
	-isystem./gsl $(INCFLAGS) \
	-std=gnu++98 -DM=1.0 

LDFLAGS = $(PROFILE) $(OPENMP) -L$(HOME)/local/lib -L./gsl 
LIBS = -lm -lgsl -lgslcblas -lz

CXX = g++ 
			
RM = rm -f

.PHONY: clean depend all

all: trajectory field

trajectory: $(TRAJECTORY_OBJS)
	$(CXX) -o trajectory $(LDFLAGS) $(TRAJECTORY_OBJS) $(LIBS)  

field: $(FIELD_OBJS)
	$(CXX) -o generator $(LDFLAGS) $(FIELD_OBJS) $(LIBS)

timestamp.o: $(TRAJECTORY_SRCS)

*.o: Makefile

tags: $(ALL_SRCS)
	ctags $(ALL_SRCS)

clean:
	$(RM) *.o *~ trajectory generator tags 

clean-data:
	$(RM) *.out *.plt *.csv *.dat gmon.dat *.bak

clean-plots:
	$(RM) *.png *.eps

.deps:
depend:
	$(CXX) $(CXXFLAGS) -M $(ALL_SRCS) >.deps

include .deps                                                              
