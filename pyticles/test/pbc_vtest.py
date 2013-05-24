#! /usr/local/bin/python

""" 
    Visually test periodic boundary conditions using the
    optimised smooth particle implementation written in 
    Fortran.
    Copyright Andrew Charles 2009
    All rights reserved.
"""

SPAMCOMPLETE = True

import sys
import box
from time import time
import particles
import spam_complete_force
import forces
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
XMAX = 10 
YMAX = 10 
ZMAX = 10 
VMAX = 0.0
dt = 0.01
SPACING = 5.0
SIDE = (2,2,2)
NP = SIDE[0]*SIDE[1]*SIDE[2]
TEMPERATURE = 0.3
HLONG = 4.0
HSHORT = 2.0
RINIT = 'grid'

# The particle, neighbour list and box are global in scope
box = box.PeriodicBox(xmax=XMAX,ymax=YMAX,zmax=ZMAX)
p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
    ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX,simbox=box)
box.p = p
s = pview.ZPRView(p)
nl = neighbour_list.SortedVerletList(p,cutoff=4.0)

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
        # Any of these conditions indicates numerical instability
        if (p.vdot > 1000000).any():
            print 'isbig'
            pyglet.clock.unschedule(update)
        if np.isnan(p.vdot).any():
            print 'isnan'
            pyglet.clock.unschedule(update)
        if np.isnan(p.v).any():
            print 'isnan'
            pyglet.clock.unschedule(update)
        if (nl.rij[0:nl.nip] < 0.5).any():
            print 'isclose'
            pyglet.clock.unschedule(update)
            print nl.rij[0:nl.nip].min()
    print 'update',time() - t

def redraw(t):
    s.redraw(p)

@s.win.event
def on_draw():
    s.clear()
    s.redraw(p)
    #print 'draw',time() - t

@s.win.event
def on_key_press(symbol,modifiers):
    if symbol == pyglet.window.key.R:
        initialise()

def initialise():
    global p,nl,cnt,buttons
    print "Restarting"
    p = particles.SmoothParticleSystem(NP,maxn=NP,d=3,rinit=RINIT,vmax=VMAX
        ,side=SIDE,spacing=SPACING,xmax=XMAX,ymax=YMAX,zmax=ZMAX
        ,temperature=TEMPERATURE,hlong=HLONG,hshort=HSHORT
        ,thermostat_temp=TEMPERATURE,simbox=box)

    p.v[:] = 0.0

    print np.mean(p.t)
    nl = neighbour_list.SortedVerletList(p,cutoff=4.0)
    
    p.nlists.append(nl)
    p.nl_default = nl

    p.forces.append(spam_complete_force.SpamComplete(p,nl))
    p.forces.append(forces.CollisionForce3d(p,nl,cutoff=0.6))

    nl.build()
    nl.separations()

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




