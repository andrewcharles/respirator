#! /usr/local/bin/python

""" 
    Looking for the source of the instability in the fortran wrap.

    Worked through, eliminating each step in the calculation.

    Using cython forces, we don't get an instability that as a
    preferred direction, the way we do with the fortran code.

    There is either:
        1. a bug in sphlib
        2. a bug in the python interface
        3. a bug in the code that uses the python interface.

    Definitely having high viscosity introduces anisotropy.
    Core sigma doesn't seem to make much difference.
    The fsph run is subject to much more clumping.
    Even with the core, gradient and viscosity coeffs all set to zero

    Solved - the heat flux was the wrong sign. (Feb 2010)
    Seems it needs dwdx_iij after all.

"""

CFORCES = False
#CFORCES = True

import sys
from time import time
import particles
import spam_complete_force
import forces
if CFORCES:
    import c_forces as forces
import pyglet
from pyglet.window import mouse
import pview
import profile
import neighbour_list
from pyglet.gl import *
from properties import spam_properties
import numpy as np

# Global variables
MAX_STEPS = 100000
XMAX = 20 
YMAX = 20 
ZMAX = 20 
VMAX = 0.0
MASS = 1.0
dt = 0.05
SPACING = 2.0
SIDE = (3,3,3) 
NP = (SIDE[0]*SIDE[1]*SIDE[2])
TEMPERATURE = 0.2
HLONG = 4.0
HSHORT = 2.0
RINIT = 'grid'

p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX,
    side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX,
    integrator='rk4')
s = pview.ZPRView(p)
nl = neighbour_list.VerletList(p,cutoff=5.0)

cnt = 0
fps = 0
tstart = time()
rebuild_nl = 1

def update(t):
    global cnt,fps,rebuild_nl ,dt
    t = time()
    cnt += 1
    if cnt >= MAX_STEPS:
        pyglet.app.exit()
    else:
        p.update(dt)
        # Trying to find the cause of instability
        if (p.vdot > 1000000).any():
            print 'isbig'
            pyglet.clock.unschedule(update)
        if np.isnan(p.vdot).any():
            print 'isnan'
            pyglet.clock.unschedule(update)
        if np.isnan(p.v).any():
            print 'isnan'
            pyglet.clock.unschedule(update)
        if (nl.rij[0:nl.nip] < 0.1).any():
            print 'isclose'
            #pyglet.clock.unschedule(update)
            print nl.rij[0:nl.nip].min()
    print 'update',time() - t

def redraw(t):
    s.redraw(p)

@s.win.event
def on_draw():
    s.clear()
    s.redraw(p)

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

def initialise():
    global p,nl,cnt,buttons
    print "Restarting"
    p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
        ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
        ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT,mass=MASS
        ,thermostat_temp=TEMPERATURE,thermostat=False)

    print np.mean(p.t)
    nl = neighbour_list.NeighbourList(p)#,cutoff=HLONG)
    
    p.nlists.append(nl)
    p.nl_default = nl


    if CFORCES:
        p.forces.append(forces.SpamForce(p,nl))
        p.forces.append(forces.CohesiveSpamForce(p,nl))
        p.forces.append(forces.SpamConduction(p,nl))
        particles.SPROPS = True
    else:
        p.forces.append(spam_complete_force.SpamComplete(p,nl,
            eta=0.0,zeta=0.0,cgrad=0.0,
            sigma=0.0,rcoef=0.0))

    #p.forces.append(forces.FortranCollisionForce(p,nl,cutoff=0.5))

    nl.build()
    nl.separations()

    # Use the python spam props to initialise
    spam_properties(p,nl)
    print 'initial mean temperature',np.mean(p.t)
    print 'initial mean density',np.mean(p.rho)

    cnt = 0

def main():
    initialise()
    pyglet.clock.schedule_interval(update,0.05)
    pyglet.clock.schedule_interval(redraw,0.2)
    pyglet.app.run()

if __name__ == "__main__":
    main()




