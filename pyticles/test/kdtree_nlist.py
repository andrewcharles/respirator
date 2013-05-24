""" Does the scipy kd neighbour list do a better job? """
""" http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdista
"""

import particles
import numpy as np
from scipy.spatial import KDTree, cKDTree
from scipy.spatial.distance import pdist, squareform
from pairsep import fast_pairsepr
from time import time

# Global variables
XMAX = 10 
YMAX = 10 
ZMAX = 10 
dt = 0.01
NP = 1000

p = particles.ParticleSystem(NP,maxn=NP,rinit='random'
    ,xmax=XMAX,ymax=YMAX,zmax=ZMAX)

tsep = time()
kt = KDTree(p.r)
pairs = kt.query_pairs(4.0,eps=0.1)
print time() - tsep

# cKDTree does not support the pairs function
# The implementation below does not achieve what the
# above code does.
#kt = cKDTree(p.r)
#cpairs = kt.query(p.r,distance_upper_bound=4.0)
pa = np.array(list(pairs))
dr,ds = fast_pairsepr(pa[:,:],p.r[:,:])

# This line does a brute force using scipy.spatial.distance
dist = pdist(p.r[:,:])
rij = np.zeros(dist.shape)
rijs = squareform(dist)
k = 0
for i in range(NP):
    for j in range(i+1,NP):
        rij[k] = rijs[i,j]
        k += 1

#ridx = pa[:,0]
#rjdx = pa[:,1]
#ridxdx0 = ridx == ridx[0]
#ri0 = 

