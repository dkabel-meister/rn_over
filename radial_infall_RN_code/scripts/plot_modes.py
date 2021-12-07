#!/usr/bin/env python

import pylab
from glob import glob
from math import fabs
from optparse import OptionParser

#parser = OptionParser()
#parser.add_option("--l-min",        type = "int",    dest = "l_min",    help = "minimum l value")

run     = "runs_staticQ0p9"
dirname = "./data/fields/from50M/"+run
lsplit = 5

#run = "runs_hubeney_chop128"
#dirname = "./data/"+run
#lsplit = 3

filenames = dirname + '/psi_[0-9]*.out' 
files = glob(filenames)

lvals     = []
fieldvals = {}
rvals     = []

rplus = 1.4359

files.sort()
for ff in files:
    # get the ell value from the filename
    lval = int( ff.split('/')[lsplit].split('_')[1].rsplit('.')[0])
    lvals.append(lval)

    fieldvals[lval] = []
    # get the data
    for line in open(ff, 'r').readlines():
        if line[0] == '#': # skip comments headers 
            continue

        if lval == 0:
            rvals.append( map(float, line.split())[0]/2 )

        fieldvals[lval].append( fabs(map(float, line.split())[3]) ) # get phi value at line


# plot the data 
for l in lvals:
    if l<11:
        pylab.loglog(rvals, fieldvals[l])
        pylab.xlabel('$\log_{10} (r/r_{+}) $', fontsize=16)
        pylab.ylabel('$\log_{10}{(F_{tr}^{\mathrm{R}})_\ell}$', fontsize=16)       
        pylab.title("Modes for $r_0=50$ $Q=0.9$ $v_0=0$ $q=1$ $m=1$", fontsize=15)
pylab.show()
    




