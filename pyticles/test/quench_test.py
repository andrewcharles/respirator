""" Set up a system with some inhomegeneities and look at 
    the forces.
"""

# Create a particle system, nlist and add a force.
import matplotlib
import math
import numpy as np
from pylab import *
import sys
from spkernel import *
import splib
import fkernel
import eos
import feos
import sphforce3d
from time import time
import particles
from mpl_toolkits.mplot3d import axes3d

import spam_complete_force
import forces

import pview
import profile
import neighbour_list
from properties import spam_properties
ion()

# Global variables
XMAX = 10 
YMAX = 10 
ZMAX = 10 
VMAX = 0.0
dt = 0.000001
SPACING = 1.2
SIDE = (5,5,5)
NP = SIDE[0]*SIDE[1]*SIDE[2]
TEMPERATURE = 1.0
RINIT = 'grid'
HLONG = 5.0
HSHORT = 2.5
A = 7.45e+04 #2.0
B = 5.84e-01 #0.5
KB = 3.29e+04 #1.0
MASS = 0.3 #1.386e-01
K = 9.2e06
R = 0
S = 10

cnt = 0
fps = 0
tstart = time()

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
    ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT
    ,thermostat_temp=TEMPERATURE,mass=MASS)
nl = neighbour_list.VerletList(p,cutoff=10.0)
p.nlists.append(nl)
p.nl_default = nl

nl.build()
nl.separations()

p.forces.append(spam_complete_force.SpamComplete(p,nl,rcoef=R,sigma=S,
    adash=A,bdash=B,kbdash=KB,thermalk=K,cgrad=0.0))

p.update(dt)

print 'rhomin',p.rho[0:p.n].min()
print 'rhomean',p.rho[0:p.n].mean()
print 'rhomax',p.rho[0:p.n].max()

ax = axes3d.Axes3D(figure())
#ax = gca()
x = p.r
a, ascl = p.vdot, 1.0
v, vscl = p.v, 100.0

for i in range(p.n):
    # plot position
    ax.plot(x[:,0],x[:,1],x[:,2],'o' )
    # plot velocity
    #ax.plot( (x[i,0],x[i,0] + vscl * v[i,0]),
    #         (x[i,1],x[i,1] + vscl * v[i,1]),
    #         (x[i,2],x[i,2] + vscl * v[i,2]),'b-')

    # plot acceleration
    ax.plot( (x[i,0],x[i,0] + ascl * a[i,0]),
             (x[i,1],x[i,1] + ascl * a[i,1]),
             (x[i,2],x[i,2] + ascl * a[i,2]),'r-')

ax.set_xlim3d(0,10)
ax.set_ylim3d(0,10)
ax.set_zlim3d(0,10)




