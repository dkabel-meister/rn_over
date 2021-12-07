#!/usr/bin/env python
#-------------------------------------------------------------
# scalar_convergence.py
#-------------------------------------------------------------
# Computes the convergence factor and re-scaled "regression" 
# curves of the phi_1(2) code and plots the result
#-------------------------------------------------------------
# File created on July 1st, 2010 
# Last Modified August 6, 2010
# Author: Peter Zimmerman
#--------------------------------------------------------------

import os
import sys
import math
import pylab 
from optparse import OptionParser

#params =  {'text.usetex': True }  # gives tex functionality for plotting
#pylab.rcParams.update(params)

# command line arguments 
parser = OptionParser()
parser.add_option("-e", "--error",  action = "store_true",  dest = "do_error", default=False, help = "computing the error")
parser.add_option("--dicings",      type = "string", dest = "chops",    help = "list of decimations (must be three even numbers)")
parser.add_option("--pow",          type = "int",    dest = "power",    help = "power of regression")  
parser.add_option("--file-stem",    type = "string", dest = "stem",     help = "filename basename")
parser.add_option("--field",    type = "string", dest = "field",    help = "specifies the field file header"+ "\n" "phi_field_fh for phi and field_fh for psi")
parser.add_option("--field-params", type = "string", dest = "sim_opts", help = "list of simulation options")  
#
# An example of a valid command line argument is
# python convergence.py --dicings 8,4,2 --pow 2 --file-stem psi_regression --field field_fh --field-params \
#  "--dt 0.125 --tmax 600 --el 2 --r0 100 --rdot0 0.0 --q 1 --observer 200"
#

# parse the arguments 
(opts, args) = parser.parse_args()

stem     = opts.stem      # filestem
sim_opts = opts.sim_opts  # simulation options
field_fh = opts.field     # field header 
power    = opts.power     # power of the regression
do_err   = opts.do_error  # boolean for error 

chops = [ int(c) for c in opts.chops.strip().split(',') ] # chops are our decimates

panic = lambda *x:sys.stdout.write(" ".join(map(str,x)))  # complaint function 

#
# Things to complain about
#

length_chops = chops.__len__()
# not 3 decimates 
if length_chops != 3:
    print "ERROR: %d is not a valid number of chops. Three chops are needed to compute the convergence factor", length_chops
    sys.exit(1)
# chops are odd
'''if filter(lambda c: c%2, chops) != []:
    panic('ERROR: --dicings must be even.', '\n')
    sys.exit(1)'''
# chops input in the wrong order
if chops[0] < chops[1] or chops[1] < chops[2]:
    panic('ERROR: --dicings must be arranged in order of decreasing even integer values.', '\n')
    sys.exit(1)

# storage for plotting
con_factor_data = [] 
d10_data        = [] 
d21_data        = [] 
err_data        = []

files = [] # stores the files
COL = 1    # data is in column 2 but zero offset makes it a 1

vals = {"f0":[], "f1":[], "f2":[]} # store the results
time = []

# helpers for string formating
bs = ' '
us = '_'

try:
    if do_err:
        outfile = open("psi_" + us.join(map(str,chops)) + "_err.dat", 'w') #output file
    else:
        outfile = open("psi_" + us.join(map(str,chops)) + "_convergence.dat", 'w') #output file
except IOError:
    print "The file " + outfile + " file does not exist, exiting"
    raise

# run method
run_cmd = os.system

#
# run the simulation code at different step decimations
#

for c in chops:
    cs = str(c)
    # string that will be executed in the shell  
    run_string = './generator --chop ' + cs + ' --' +  field_fh + bs + stem + us + cs + '.out' \
        + bs + sim_opts 
    # print out some information about the run
    panic(run_string, '\n', "Running at a step decimation of", c, '\n')
    # run the simulation at decimation chop
    run_cmd(run_string)
    # pack list of filenames
    files.append(stem + us + cs + '.out')

f_idx = 0 # file index 

for fname in files: 
    lines = open(fname,'r').readlines()
    for line in lines:
        if line[0] == '#': # skip comments headers 
            continue
        # put the data into the list 
        phi_val = map(float,line.split())[COL]
        rstar_val =  map(float,line.split())[0]
        if phi_val != 0:
            vals["f"+str(f_idx)].append(phi_val) 
            if f_idx == 2:
                time.append(rstar_val)
    f_idx += 1 

# rstar values [should be same for each file...if not, then you're screwed]
#time = [map(float,line.strip().split())[1] for line in open(files[0],'r').readlines() if not line[0] == '#']

# we don't want to overstep arrays
length = min(len(vals["f0"]), len(vals["f1"]), len(vals["f2"]))

#
# compute the convergence factor
# 

idx = 0
idx_max = length
# loop over all lines in the data column
while (idx < idx_max):
    try: 
        # 
        # con_frac = (f_finest - f_medium)/(f_medium - f_coarsest)
        #
        con_frac = (vals["f2"][idx] - vals["f1"][idx]) / (vals["f1"][idx]-vals["f0"][idx])  # convergence factor fraction
        d10      = (vals["f1"][idx] - vals["f0"][idx]) / (2**power - 1)                     # rescaled differences 
        d21      = (vals["f2"][idx] - vals["f1"][idx]) / (2**(2*power) - 2**power)          # rescaled differences

        if do_err:
            err  = (vals["f0"][idx] - vals["f1"][idx]) / (1.0 - 2.0**power)

    except ZeroDivisionError:
        print "\nDivision by zero in the convergence fraction occurred. Breaking at iteration number %d of %d.\n" %(idx, idx_max)
        break
    try: 
        con_factor = math.log(math.fabs(con_frac))/math.log(2)
    except OverflowError:
        print "\nOverflow in the convergence fraction occurred. Breaking at iteration number %d of %d.\n" %(idx, idx_max)
        break                                                                           

    # pack arrays for plotting
    con_factor_data.append(con_factor)
    d10_data.append(d10)
    d21_data.append(d21) 

    if do_err: 
        err_data.append(math.fabs(err))

    # print results to file
    if do_err:
        print >>outfile, " %22.17f\t %22.17f\t %22.17f\t %22.17f\t %22.17e " %(time[idx], con_factor, d10, d21, err)
    else:
        print >>outfile, "%22.17f\t %22.17f\t %22.17f\t %22.17f\t" %(time[idx], con_factor, d10, d21)

    # increment
    idx += 1

# close the file
outfile.close()  


#
# plot the data 
#

s = 'psi'

if do_err:
    if (len(time) == len(err_data)):
        pylab.semilogy(time, err_data)
        pylab.grid(True)
        pylab.xlabel('$t/M$', fontsize=20)
        pylab.ylabel('resolution error', fontsize=14)
        pylab.title(sim_opts, fontsize=12)
        pylab.show()
    else:
        panic('\n', 'ERROR: lengths of data arrays are not the same. Plotting failed', '\n', 'Lengths are: \t')
else: 
    if (len(time) == len(d21_data) and len(time) == len(d21_data) and len(time) == len(con_factor_data)):
        pylab.subplots_adjust(hspace=0.35)
        pylab.subplot(211)                  
        pylab.plot(time, d10_data)
        pylab.plot(time, d21_data)
        pylab.ylim(-0.04,0.04)
        pylab.grid(True)
        pylab.xlabel('$t/M$', fontsize=20)
        pylab.ylabel('re-scaled field values', fontsize=14)
        pylab.title(sim_opts, fontsize=12)
        pylab.legend( ['$8-4$', '$16-8$'] ) 
        pylab.subplot(212)        
        pylab.plot(time, con_factor_data)
        pylab.ylim(0,6)
        pylab.xlabel('$t/M$', fontsize=20)
        pylab.ylabel('convergence factor', fontsize=14)
        pylab.grid(True)
        #pylab.savefig("fig.png")
        pylab.show()
    else:
        panic('\n', 'ERROR: lengths of data arrays are not the same. Plotting failed', '\n', 'Lengths are: \t')
        print "len(time) = %d  len(d21_data) = %d len(d10_data) = %d len(con_factor_data = %d" %(len(time),len(d21_data), len(d10_data), len(con_factor_data))
