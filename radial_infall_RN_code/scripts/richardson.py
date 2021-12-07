#!usr/bin/env python
#-----------------------------------------------------
# Runs simulation from l-min to l-max and computes 
# the Richardson Extrapolation of the retarded field
# 
#
# Example of command line arguments: 
# python richardson.py --l-min 0 --l-max 25 --file-stem psi --field multipole_fh \
#  --sim-metadata aug26a --field-params "--dt 0.125 --tmax 2046  \
# --r0 125 --rdot0 -0.01  --q 1 --Q 0.5"
#---------------------------------------------------
# Author: Peter 
# Date Created : Mon Jan 10, 2011
#---------------------------------------------------
import os 
import sys
import glob
from optparse import OptionParser
import commands
from math import fabs
import pylab

parser = OptionParser()

parser.add_option 
parser.add_option("-r", "--run",  action = "store_true",  dest = "run_and_plot", default=True, help = "boolean for run and plot")
parser.add_option("-p", "--plot", action = "store_false", dest = "run_and_plot", help = "boolean for plot only")


parser.add_option("--l-min",        type = "int",    dest = "l_min",    help = "minimum l value")
parser.add_option("--l-max",        type = "int",    dest = "l_max",    help = "maximum l value")
parser.add_option("--file-stem",    type = "string", dest = "stem",     help = "filename basename")
parser.add_option("--field",        type = "string", dest = "field",    help = "specifies the field file header"+ "\n" "phi_field_fh for phi")
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
     
def richardson(ell, power, files, chops):

    COL = 1 # column for the F_tr

    get_col = lambda x,y : map(float,x.split())[y]

    dirnameTemp = "data/runs_" + sim_id + "_richardson_Ftr"

    if (os.path.exists(dirnameTemp) == True):
        pass
    else:
        os.mkdir(dirnameTemp)   
    vals = {"f0":[], "f1":[]} # store the results

    regs = {"A": {"a0":[],"a1":[]}, "B":{"b0":[],"b1":[]}, "D":{"d0":[],"d1":[]}} # store the results
    time = []    
    el_files = []

    for fname in files:
        el_file = fname.split('/')[2].split('_')[1]
        if el_file == str(ell):
            el_files.append(fname)
    f_idx=0


    for ff in el_files:
        lines = open(ff,'r').readlines()
        for line in lines:
            if line[0] == '#': # skip comments headers 
                continue
            # put the data into the list 
            t_val =  get_col(line,0)
            vals["f"+str(f_idx)].append(get_col(line,1))
            regs["A"]["a"+str(f_idx)].append(get_col(line,2))
            regs["B"]["b"+str(f_idx)].append(get_col(line,3))
            regs["D"]["d"+str(f_idx)].append(get_col(line,4))

            if f_idx == 0:
                time.append(t_val)
        f_idx +=1      

    # helpers for string formating
    bs = ' '
    us = '_'

    try:
        outfile = open(dirnameTemp+"/Ftr_" + str(ell) + ".dat", 'w') # output file
    except IOError:
        print "The file " + outfile + " file does not exist, exiting"
        raise      

    # we don't want to overstep arrays
    if len(vals["f0"])!=len(vals["f1"]):
        panic("File lengths differ for fields", "\n")
        sys.exit(1)
    
    length = min(len(vals["f0"]), len(vals["f1"]))

    #
    # compute the extrapolation
    # 

    idx = 0
    idx_max = length                           
    extrapolatedField = []

    # loop over all lines in the data column
    while (idx < idx_max):
        try: 
            res  = (chops[0]**power*vals["f0"][idx] - chops[1]**power*vals["f1"][idx]) / (chops[0]**power-chops[1]**power)

            resA = (chops[0]**power*regs["A"]["a0"][idx] - chops[1]**power*regs["A"]["a1"][idx]) / (chops[0]**power-chops[1]**power)
            resB = (chops[0]**power*regs["B"]["b0"][idx] - chops[1]**power*regs["B"]["b1"][idx]) / (chops[0]**power-chops[1]**power)
            resD = (chops[0]**power*regs["D"]["d0"][idx] - chops[1]**power*regs["D"]["d1"][idx]) / (chops[0]**power-chops[1]**power)
        except ZeroDivisionError:
            print "\nDivision by zero occurred. Breaking at iteration number %d of %d.\n" %(idx, idx_max)
            break
        except OverflowError:
            print "\nOverflow in occurred. Breaking at iteration number %d of %d.\n" %(idx, idx_max)
            break

        print >>outfile, "%.16e %.16e %.16e %.16e %.16e" %(time[idx], res, resA, resB, resD)
        idx+=1

    outfile.close()


panic = lambda *x : sys.stdout.write(" ".join(map(str,x)))  # complaint function 


if (os.path.exists('data') == False):
   panic("A directory named data is required to exist", "\n")
   sys.exit(1)
if (field_fh != 'multipole_fh'):
    panic("ERROR: Invalid data stream to jobs.py. Regularization only valid on data at the particle", '\n')
    sys.exit(1)

dirname = "data/runs_" + sim_id 

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

chops = [2,4]
if (run_and_plot == True):
    for l in l_vals:
        ls = str(l) 

        for c in chops:
            cs = str(c)
            name_tmp = dirname + '/' + stem + us + ls + us +cs 

            run_string = './generator --chop ' + cs + ' --el ' + ls + ' --' +  field_fh  \
                +  bs + name_tmp + '.out' +  bs + sim_opts 
            # output useful information
            panic("Running for l=", ls, '\n', run_string, '\n',)
            # run the simulation 
            run_cmd(run_string)     

filenames = dirname + '/'+ stem + us + '*[0-9]' + us + '[0-9].out' 
power=2

files = glob.glob(filenames)
files.sort()

lvals   = []
phivals = []
phivals_A   = []
phivals_AB  = []
phivals_ABD = []      


if (run_and_plot == True):
    for ff in files:
        # get the ell value from the filename
        lval = int(ff.split('/')[2].split('_')[1])
        # call the extrolation routine
        richardson(lval, power, files, chops)

os.chdir("data/runs_" + sim_id + "_richardson_Ftr")
#richard_loc = "data/runs_" + sim_id + "_richardson_Ftr/" + "*.dat"
richard_files = glob.glob('*.dat')
richard_files.sort() 

for ff in richard_files:
    # get the ell value from the filename
    line_idx = 0      
    lval = int( ff.split('_')[1].split('.')[0] )
    lvals.append(lval)
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

# 
# plot the data 
#
if (len(lvals) == len(phivals)):
    pylab.semilogy(lvals, phivals, 'bo', ms=8)
    pylab.semilogy(lvals, phivals_A, 'r^', ms=8)
    pylab.semilogy(lvals, phivals_AB, 'gs', ms=8)
    pylab.semilogy(lvals, phivals_ABD, 'm+', ms=8)
    pylab.xlabel('$\ell$', fontsize=16)
    pylab.ylabel('Field', fontsize=16)
    pylab.legend(['$(F_{tr})_\ell$','$(F_{tr})_\ell - A$', '$(F_{tr})_\ell - A - B$', r'$(F_{tr})_\ell - A - B - D$'], 'upper right')
    '''pylab.legend(['$(F_{tr})_\ell$','$(F_{tr})_\ell-A(\ell+1/2)$', '$(F_{tr})_\ell-A(\ell+1/2)-B$', 
    r'$(F_{tr})_\ell-A(\ell+1/2)-B-D/((\ell-1/2)(\ell+3/2))$'], 'lower left')'''
    pylab.title(sim_opts, fontsize=10)
    pylab.grid(True)
    pylab.show()
else:
   panic('\n', 'ERROR: lengths of data arrays are not the same. Plotting failed', '\n', 'Lengths are: \t')  
   print "len(lvals) = %d  len(phivals) = %d" %(len(lvals), len(phivals)),  
