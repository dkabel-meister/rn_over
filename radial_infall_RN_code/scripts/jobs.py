#!usr/bin/env python
#---------------------------------------------------
# Runs simulation from l-min to l-max, saves
# a single line in each output file.
# Plots the results (ell vs. field)
#
# Example of command line arguments: 
# python jobs.py -r --l-min 0 --l-max 25 --file-stem psi --field multipole_fh \
#  --sim-metadata aug26a --line 405 --field-params "--dt 0.125 --tmax 257  \
# --r0 125 --rdot0 -0.01 --chop  2 --q 1"
#---------------------------------------------------
# Author: Peter 
# Date Created : Aug 23 2009
#---------------------------------------------------
import os 
import sys
import glob
import pylab 
from optparse import OptionParser
from math import fabs

panic = lambda *x : sys.stdout.write(" ".join(map(str,x)))  # complaint function 
                                                                          
parser = OptionParser()

#parser.add_option 
parser.add_option("-r", "--run",  action = "store_true",  dest = "run_and_plot", default=True, help = "boolean for run and plot")
parser.add_option("-p", "--plot", action = "store_false", dest = "run_and_plot", help = "boolean for plot only")

parser.add_option("--l-min",        type = "int",    dest = "l_min",    help = "minimum l value")
parser.add_option("--l-max",        type = "int",    dest = "l_max",    help = "maximum l value")
parser.add_option("--file-stem",    type = "string", dest = "stem",     help = "filename basename")
parser.add_option("--field",        type = "string", dest = "field",    help = "specifies the field file header"+ "\n")
parser.add_option("--field-params", type = "string", dest = "sim_opts", help = "list of simulation options")  
parser.add_option("--sim-metadata", type = "string", dest = "sim_id",   help = "date of simulation")  
parser.add_option("--line",         type = "int",    dest = "line",     help = "line number for the data")

# parse the arguments 
(opts, args) = parser.parse_args()

run_and_plot = opts.run_and_plot

stem     = opts.stem      # filestem
sim_opts = opts.sim_opts  # simulation options
sim_id   = opts.sim_id    # simulation id
field_fh = opts.field     # field header 
l_min    = opts.l_min     # min ell
l_max    = opts.l_max     # max ell
l_dump   = opts.line      # dump line
     
if (os.path.exists('data') == False):
   panic("A directory named data is required to exist", "\n")
   sys.exit(1)
if (field_fh != 'multipole_fh'):
    panic("ERROR: Invalid data stream to jobs.py. Regularization only valid on data at the particle", '\n')
    sys.exit(1)

#dirname = "data/runs_" + sim_id
dirname = "data/" + sim_id

if (os.path.exists(dirname) == True):
    pass
else:
    os.mkdir(dirname)

l_vals = [l_min + l for l in range(0,l_max-l_min+1)]

# run method
run_cmd = os.system

# helpers for string formating
bs = ' '
us = '_'   

if (run_and_plot == True):
    for l in l_vals:
        ls = str(l)
        # add a zero before to make filenames sortable
        if l < 10:
            ls = '0'+ls
        name_tmp = dirname + '/' + stem + us + ls 
        
        run_string = './generator --el ' + ls + ' --' +  field_fh  \
                +  bs + name_tmp + '.out' +  bs + sim_opts 

        # output useful information
        panic("Running for l=", ls, '\n', run_string, '\n',)
        # run the simulation 
        run_cmd(run_string)     

filenames = dirname + '/psi_[0-9]*.out' 
files = glob.glob(filenames)

lvals       = []
phivals     = []
phivals_A   = []
phivals_AB  = []
phivals_ABD = []

for ff in files:
    # get the ell value from the filename
    line_idx = 0
    lval = int( ff.split('/')[2].split('_')[1].rsplit('.')[0])
    lvals.append(lval)
    # get the data
    for line in open(ff, 'r').readlines():
        if line[0] == '#': # skip comments headers 
            line_idx += 1 
            continue
        # append only if line number is the dump line number
        phi_val     = map(float, line.split())[1] # get phi value at line
        phi_val_A   = map(float, line.split())[2] # get phi^R (A term) value at line
        phi_val_AB  = map(float, line.split())[3] # get phi^R (A and B terms) value at line
        phi_val_ABD = map(float, line.split())[4] # get phi^R (A, B and D terms) value at line
        if (line_idx == l_dump):
            phivals.append(fabs(phi_val))
            phivals_A.append(fabs(phi_val_A))
            phivals_AB.append(fabs(phi_val_AB))
            phivals_ABD.append(fabs(phi_val_ABD))
        line_idx += 1 

# plot the data 
if (len(lvals) == len(phivals)):
    pylab.semilogy(lvals, phivals, 'bo', ms=8)
    pylab.semilogy(lvals, phivals_A, 'r^', ms=8)
    pylab.semilogy(lvals, phivals_AB, 'gs', ms=8)
    pylab.semilogy(lvals, phivals_ABD, 'm+', ms=8)
    pylab.xlabel('$\ell$', fontsize=16)
    pylab.ylabel('Field', fontsize=16)
    pylab.legend(['$r^2 (F_{tr})_\ell$','$r^2 (F_{tr})_\ell - A$', '$r^2(F_{tr})_\ell - A - B$','$r^2(F_{tr})_\ell - A - B -D$'], 'lower left')
    '''pylab.legend(['$(F_{tr})_\ell$','$(F_{tr})_\ell-A(\ell+1/2)$', '$(F_{tr})_\ell-A(\ell+1/2)-B$', 
    r'$(F_{tr})_\ell-A(\ell+1/2)-B-D/((\ell-1/2)(\ell+3/2))$'], 'lower left')'''
    pylab.title(sim_opts, fontsize=10)
    pylab.grid(True)
    pylab.show()
else:
   panic('\n', 'ERROR: lengths of data arrays are not the same. Plotting failed', '\n', 'Lengths are: \t')  
   print "len(lvals) = %d  len(phivals) = %d" %(len(lvals), len(phivals))
