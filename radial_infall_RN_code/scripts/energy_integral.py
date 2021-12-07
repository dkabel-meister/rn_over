#!/usr/bin/env python 
from math import sqrt, fabs
import sys

# input file
infile = sys.argv[1]       
# line number where the data begins to look clean
good_data_start_idx = int(sys.argv[2])

# trajectory parameters
r0 = 3500.0
E  = 1.700005e-01

# spacetime parameters 
M = 1

# hubeny parameters
eps = 0.1
a   = 2.0
b   = 1.5
c   = 1.0

# derived physical parameters
Q = M*(1.0 - 2.0*eps*eps)
q = a*eps
m = c*eps


# read data from file
lines = [x.strip() for x in open(infile, 'r').readlines()]
rvals = [float(x.split()[0]) for x in lines]        # r
r2Ftrvals = [float(x.split()[1]) for x in lines]    # r^2 F_{tr}

Ftrvals = [] # unscaled regular field values

# the unscaling (multiplication by 1/r^2) 
for x in zip(rvals,r2Ftrvals):
    Ftrvals.append(x[1]/x[0]**2)

# create lists for results
trapVals = [x*0 for x in range(0,len(rvals))]  # integration values
R2       = [x*0 for x in range(0,len(rvals))]  # R(t)^2 (un-modified traj)
R2g      = [x*0 for x in range(0,len(rvals))]  # generalized R(t)^2 (modified traj)


outfile = open('corrTraj_' + str(eps) +'_'+ str(a) +'_'+ str(b) +'_'+ str(c) + '_r0'+str(r0)+'.dat', 'w')
#outfile = open('corrTraj_' + infile + '.dat', 'w')

# integration of F_tr along the trajectory
for i in range(1,len(rvals)):
    if i<good_data_start_idx:
        pass
    else:
        trapVals[i] = trapVals[i-1] + 1/2.*(rvals[i-1]-rvals[i])*(Ftrvals[i-1]+Ftrvals[i]) # trapazoidal rule
        general_E = E+q*trapVals[i]
        R2[i] = 1/m**2 * ((E - q*Q/rvals[i])**2 - m**2*(1-2*M/rvals[i]+Q**2/rvals[i]**2))
        R2g[i] = 1/m**2 * ((general_E - q*Q/rvals[i])**2 - m**2*(1-2*M/rvals[i]+Q**2/rvals[i]**2))
        print >>outfile, "%22.17f\t %22.17f\t %22.17f" %(rvals[i], R2[i], R2g[i])

# plot the results
import pylab

pylab.plot(rvals, R2)
pylab.plot(rvals, R2g, ls='--')
pylab.xlabel('$r/M$', fontsize=18)
pylab.ylabel('$\dot{r}^2$', fontsize=18)
title='Trajectory for Q='+str(Q)
pylab.title(title)
pylab.legend(('No SF','SF')) 

pylab.show()
