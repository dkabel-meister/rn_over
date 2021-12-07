#!/usr/bin/env python
import sys
from math import sqrt, fabs

eps = 0.1
M   = 1.0
Q   = M*(1.0 - 2.0*eps**2)

a = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])

q  = a*eps
E0 = a*eps-2*b*eps**2
m  = c*eps

def f(r):
    return 1.0 - 2*M/r + Q*Q/(r*r)    


def r2_F_tr_local(r):
    return -4*Q*pow(q,4)/(3.*pow(m,3)*pow(r,1)) * sqrt( (E0-q*Q/r)**2 - m**2*f(r) )  
    
r_plus = M + sqrt(M**2-Q**2)

chop = 256.0
      
r_max = 100.0        
i_max = int(100*chop)+1

outfile = open('Ftr_local' + str(eps) +'_'+ str(a) +'_'+ str(b) +'_'+ str(c) + 'r0'+str(r_max)+'.dat', 'w')

for i in range(1,i_max):
    r = r_max-i/chop
    if r > r_plus:
        print >> outfile, "%22.17f %22.17f" %(r, r2_F_tr_local(r))

    
