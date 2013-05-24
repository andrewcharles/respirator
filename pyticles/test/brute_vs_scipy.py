""" Compare the scipy.spatial.distance.pdist with
    regular brute to make sure they are equivalent.
    Create random particle positions for NP particles
"""
import particles
import neighbour_list
from time import time
import numpy as np

NP = 100
n = NP

p = particles.ParticleSystem(n,d=3,maxn=n)
bl = neighbour_list.NeighbourList(p)
pl = neighbour_list.BruteScipy(p)
for i in range(n):
        p.r[i,:] = np.random.random(3) * 10

dt_tot = {}
def time_list(nl):
    dt_tot[nl] = 0
    nl.build()
    nl.compress()
    t = time()
    nl.separations()
    # execution time is in seconds
    dt = (time() - t)
    print nl.nip,'%7.6f' %(dt)
    dt_tot[nl] += dt

time_list(bl)
time_list(pl)

for key in dt_tot.keys():
    '%7.6f' %(dt_tot[key])



