import sys
import particles
import spam_complete_force
import forces
import pyglet
from pyglet.window import mouse
import pview
import neighbour_list
from properties import spam_properties
import numpy as np

# Global variables
XMAX = 10 
YMAX = 10 
ZMAX = 10 
dt = 0.01
NP = 50


# First test, for consistency
# Build the list, then compress it (which should limit the number of 
# neighbours.
p = particles.ParticleSystem(NP,maxn=NP,rinit='random'
    ,xmax=XMAX,ymax=YMAX,zmax=ZMAX)

limlist = neighbour_list.LimitedVerletList(p,cutoff=4.0)
limlist.max_neighbours=8

limlist.build()
limlist.compress()
print 'list compressed',max(limlist.pn)
print limlist.nip

limlist.build_pnlist(max_neighbours=10)
print 'list built',max(limlist.pn)
print limlist.nip
limlist.remove_pnlist_duplicates()
print 'dupes gone',max(limlist.pn)
limlist.map_pnlist_to_iap()
print limlist.nip

print limlist.nip
print(max(limlist.pn))

