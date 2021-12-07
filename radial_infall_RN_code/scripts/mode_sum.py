#!usr/bin/env python
#---------------------------------------------------
# Computes the regular field (residual)
# using the mode sum subraction technique.
# NOTE: all files MUST be the same length and 
# have the same r values for this to work!
#---------------------------------------------------
# Author: Peter and Roland Haas (gawk)
# Date Created : Sept 14, 2010
#---------------------------------------------------  
import sys
import pylab 
import operator
import commands
from optparse import OptionParser  

panic = lambda *x : sys.stdout.write(" ".join(map(str,x)))  # complaint function 

parser = OptionParser()
parser.add_option("--sim-loc", type = "string", dest = "sim_loc", help = "location of simulation data.")  
parser.add_option("--traj-loc", type = "string", dest = "traj_loc", help = "location of trajectory data.")  
parser.add_option("--tail-loc", type = "string", dest = "tail_loc", help = "location of tail data.")  
parser.add_option("--use-tail", action = "store_false", dest = "tail", default=False, help = "boolean for remainder")  


(opts, args) = parser.parse_args()  # parse the arguments

# data locations
# where to find data (an example is ./data/runs_oct6r0/) #
sim_loc  = "./"+opts.sim_loc              

tail = opts.tail
q = float( commands.getoutput('grep -o \'\--q[[:space:]][0-9]\{1\}\.[[:digit:]]\{1\}\' %s | gawk \'{print $2}\'' %(opts.traj_loc)) )
m = float( commands.getoutput('grep -o \'\--m[[:space:]][0-9]\{1\}\.[[:digit:]]\{1\}\' %s | gawk \'{print $2}\'' %(opts.traj_loc)) )

# open file
traj_loc = open("./"+opts.traj_loc, 'r')

# create a list of a list of all the lines in all the files
# the first index of the list correspond to the file
# the second index of the list corresponds to the lines of that file

file_list =  commands.getoutput('find  %s/psi_[0-9]*.out '%(sim_loc)).split()
file_list.sort()

val_list = [[g.strip() for g in open(f,"r").readlines() if g[0] != '#']\
                       for f in file_list] 

if tail:
    tail_loc = open("./"+opts.tail_loc, 'r') 
    tail_list = [x.strip() for x in tail_loc.readlines() if x[0] != '#']

traj_list = [x.strip() for x in traj_loc.readlines() if x[0] != '#']

# Hack to get the energy
#E = float(commands.getoutput('grep \'E=\' %s/psi_10.out | gawk \'BEGIN{FS="[ =]"} /^# val:/{print $4}\'' %(sim_loc)))

nfiles = len(val_list)    # number of files
nlines = len(val_list[0]) # number of lines in the files

# column for the field
field_col = 4  # \psi - (el+1/2)*A -  B - D
tail_col  = 5
r_col     = 0   # column for the radial coordinates 

field_vals = {} # storage for field values

r_vals     = [] 
F_tr_vals  = []
force_vals = []     # f^r the radial component of the force

f0_dtdtau_vals = [float(x.split()[5]) for x in traj_list]

if tail:
    remainder_vals = [float(x.split()[1]) for x in tail_list]

for j in range(nlines):                                                          
    # append to the time list a value of time using the first file only
    # THIS KLUDGE  ONLY WORKS IF ALL FILES HAVE THE SAME times and same spacings
    r_vals.append(float(val_list[0][j].split()[r_col])) 

    # key dictionay by the line and make the val associated with the each 
    # key a list of the fields for that line in each file
    #
    field_vals[j] = [val_list[i][j].split()[field_col] for i in range(nfiles)] 

outfile = open("force_" + sim_loc.split('/')[3] + ".dat", 'w') #output file

for k in field_vals.keys():
    # create a temporary list of the floating point representation of the field values 
    tmpL = map(float, field_vals[k]) 
    # sum the elements in the list for each key (line of the files)
    F_tr_vals.append(reduce(operator.add,tmpL))
    force_vals.append(reduce(operator.add,tmpL))



nr = len(r_vals) # number of r values to plot
nf = len(force_vals) # number of field vales to plot

if tail:
    F_tr_plus_tail = []
    for idx in range(len(nf)):              
        F_tr_plus_tail[idx] = F_tr_vals[idx] + remainder_vals[idx]

for i in range(nf):
    force_vals[i] = -q*f0_dtdtau_vals[i]*force_vals[i]
    if tail:
        force_plus_tail[i] = -q*f0_dtdtau_vals[i]*F_tr_plus_tail[i] 
# 
# write out and plot the data
#
if (nr==nf): 
    for i in range(nf):
        print >>outfile, " %22.17e\t %22.17e\t %22.17e " %(r_vals[i], F_tr_vals[i], force_vals[i])

    pylab.subplots_adjust(hspace=0.35)
    if (tail):
        pylab.subplot(211)                  
        pylab.plot(r_vals, F_tr_plus_tail) #, linestyle='None', marker='.',markersize=3)
        pylab.xlabel('$r/M$', fontsize=18)
        pylab.ylabel('$F_{tr}^\mathsf{R}$', fontsize=18)
        pylab.grid(True)
        
        pylab.subplot(212)                  
        pylab.plot(r_vals, force_plus_tail)
        pylab.xlabel('$r/M$', fontsize=18)
        pylab.ylabel('$f^{\ \ r}$', fontsize=18)
        pylab.grid(True)
        pylab.show()
    else:
        pylab.subplot(211)                  
        pylab.plot(r_vals, F_tr_vals) #, linestyle='None', marker='.',markersize=3)
        pylab.xlabel('$r/M$', fontsize=18)
        pylab.ylabel('$r^2 F_{tr}^\mathsf{R}$', fontsize=18)
        pylab.grid(True)
        
        pylab.subplot(212)                  
        pylab.plot(r_vals, force_vals)
        pylab.xlabel('$r/M$', fontsize=18)
        pylab.ylabel('$r^2 f^{\ \ r}$', fontsize=18)
        pylab.grid(True)
        pylab.show()
else:
    panic('\n', 'PANIC: data arrays have unequal lengths. Plotting failed.')
    panic('\n', 'Lengths are: \t %s %s' %(nr,nf))
