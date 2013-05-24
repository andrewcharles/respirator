""" Compare different methods for building the neighbour list.

    1. Create a random configuration of particles
    2. Run each nlist method in turn
    3. Compare the total number of pairs, and histogram of number of pairs
        per particle.

"""

import sys
from time import time
import particles
import spam_complete_force
import forces
#import pyglet
#from pyglet.window import mouse
#import pview
import neighbour_list
from properties import spam_properties
import numpy as np

# Global variables
XMAX = 10 
YMAX = 10 
ZMAX = 10 
dt = 0.01
NP = 250

p = particles.ParticleSystem(NP,maxn=NP,rinit='random'
    ,xmax=XMAX,ymax=YMAX,zmax=ZMAX)

brute_list = neighbour_list.NeighbourList(p)
brute_list.build()
brute_list.build_pnlist()

verlet_list = neighbour_list.VerletList(p,cutoff=4.0)
verlet_list.build()
verlet_list.compress()
verlet_list.build_pnlist()

verlet_list = neighbour_list.FastVerletList(p,cutoff=4.0)
verlet_list.build()
verlet_list.compress()
verlet_list.build_pnlist()

fortran_list = neighbour_list.FortranVerletList(p,cutoff=4.0)
fortran_list.build()
fortran_list.compress()
fortran_list.build_pnlist()

limited_list = neighbour_list.LimitedVerletList(p,cutoff=4.0)
limited_list.max_neighbours=10
limited_list.build()
limited_list.compress()
limited_list.build_pnlist()
limited_list.remove_pnlist_duplicates()
print 'dupes gone'
print max(limited_list.pn)

k_list = neighbour_list.KDList(p,cutoff=4.0)
k_list.build()
k_list.compress()
k_list.build_pnlist()

p_list = neighbour_list.BruteScipy(p)
p_list.build()
p_list.compress()
p_list.build_pnlist()

print 'Brute,  Verlet,  Fortran,  Limited  KDTree Plist'
print brute_list.nip,verlet_list.nip,fortran_list.nip,limited_list.nip,k_list.nip,p_list.nip
print(max(brute_list.pn),
      max(verlet_list.pn),
      max(fortran_list.pn),
      max(limited_list.pn),
      max(k_list.pn),
      max(p_list.pn))

