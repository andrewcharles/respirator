#! /usr/local/bin/python

""" 
    A very small system used for functional testing.

    A 10x10x10 scaling length box
    Scaling length 2.81nm
    Scaling time 1ns
    Scaling mass 7520 au
    a = 7.45e+04
    b = 5.84e-01
    kb = 3.29e+04
    Number of particles = 6x6x6
    Target density = M / V
    V = 1000


    Particle mass 0.1386
    Temperature 0.95

    ## Copyright Andrew Charles, RMIT 2010                ###
    ## All rights reserved.                               ###

    todo: split into spawn system and run system.

"""

import sys
from time import time
import particles
import spam_complete_force
import neighbour_list
import f_properties as properties
from spam_nc import create_sph_ncfile, write_step, write_box
import numpy as np
import spkernel, spdensity
import forces
import c_forces
from box import PeriodicBox

max_steps = 5
XMAX = 10
YMAX = 10
ZMAX = 10
NDIM = 3
NP = 125
SIDE = (5,5,5)
VMAX = 0.0
dt = 0.0001
SPACING = 1.0
TEMPERATURE = 0.95
HLONG = 4.0
HSHORT = 2.0
thermalk = 1.0
RINIT = 'load'
sfname = '../data/runs/mini_equil_t9.nc'
ascl = 7.45e+04
bscl = 5.84e-01
kbscl = 3.29e+04
pmass = 1.386e-01
cgrad = 1000
thermostat = False
set_temperature = False
sigma = 0.0
rcoef = 0.0
eta = 0.0
zeta = 0.0
ofname = 'data/toybox.nc'
write_frequency = 1
config = None
gravity = 0
gravk = 1.0

#NP = SIDE[0]*SIDE[1]*SIDE[2]
cnt = 0
fps = 0
print max_steps
print "Initialising"
box = PeriodicBox(xmax=XMAX,ymax=YMAX,zmax=ZMAX)
p = particles.SmoothParticleSystem(
        NP,maxn=NP,
        d=3,
        rinit=RINIT,
        vmax=VMAX,
        side=SIDE,
        spacing=SPACING,
        xmax=XMAX,
        ymax=YMAX,
        zmax=ZMAX,
        source=sfname,
        temperature=TEMPERATURE,
        hlong=HLONG,
        hshort=HSHORT,
        thermostat_temp=TEMPERATURE,
        set_temperature=set_temperature,
        thermostat=False,
        mass=pmass,
        simbox=box
    )
#nl = neighbour_list.VerletList(p,cutoff=HLONG)
nl = neighbour_list.FastVerletList(p,cutoff=HLONG)
nl.max_neighbours = 20
p.nlists.append(nl)
p.nl_default = nl
p.forces.append(spam_complete_force.SpamComplete(
    p,nl,adash=ascl,bdash=bscl,kbdash=kbscl,cgrad=cgrad,eta=eta,
    zeta=zeta,sigma=sigma,rcoef=rcoef))
p.forces.append(forces.FortranCollisionForce(p,nl,cutoff=0.65))
if gravity == 1:
    print 'GRAVITY'
    p.forces.append(forces.BodyForce(p,k=gravk,direction=[0,-1.0,0]))
tstart = time() 

nl.build()
nl.separations()
properties.set_vdw_props(ascl,bscl,kbscl)
properties.spam_properties(p,nl,p.h,p.hlr)
# This is a crap way to set these parameters
p.set_vdw_properties((ascl,bscl,kbscl))


print 'Built list and calc properties',time()-tstart
cnt = 0
attribs = {'creator':'Andrew', 'log':'functional test' , 'vdwa':ascl, 'vdwb':bscl}
create_sph_ncfile(ofname,attribs,NP,NDIM,config=config)
write_box(ofname,p)
print "--------------------------------------------"
print "STEP   INT  DERIV =  PAIR + SPAM +  FORCE   "
tstartrun = time()
print p.r
write_step(ofname,p)

for i in range(max_steps):
    tstart = time()
    p.update(dt)
    if np.isnan(p.r).any():
        print 'stopping due to nan'
        break

print 'Completed',i,'steps, in',time()-tstartrun

