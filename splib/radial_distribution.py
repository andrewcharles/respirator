""" Computes the radial distribution function of a system
    of particles.
    This is just a normalised histogram of particle
    separations.

    http://www.physics.emory.edu/~weeks/idl/gofr2.html
    John Crocker, Eric Weeks

    Andrew Charles 2013
    All rights reserved.
    Research use with permission is permitted.

"""

from configuration import grid3d, random3d, fcc3d
from neighbours import build_list, separations, distances, minimum_image3d
import matplotlib.pyplot as plt
import numpy as np
import math
plt.ion()

def rdfplot(r,(x,y,z),rmax=10.):
    """ Plot the radial distribution of a 
        set of positions.
        (x,y,z) box dimensions
    """
    binw = rmax/100.
    bins = np.arange(0,rmax,binw)
    iap = build_list(r)
    rij = separations(r,iap)
    for i in range(rij.shape[0]):
        minimum_image3d(rij[i],x,y,z)
    rmax = np.sqrt(x*x+y*y+z*z)/2
    n = r.shape[0]
    r = distances(rij)
    V = x*y*z
    h,bins = np.histogram(r,bins=bins,normed=False)
    h = h*2 #required because pairs are counted only once
    # h now gives the number of particles in a shell of width binw
    # at a distance bins[i] to bins[i+1]
    rdf = np.zeros(h.shape)
    norm = np.zeros(h.shape)
    for i in range(h.size):
        # Normalise by the volume of each shell (4 pi rsquared dr)
        norm[i] = 4.0 * math.pi * binw * (0.5*(bins[i] + bins[i+1]))**2
        
        norm[i] = norm[i] * (n/V)
        if norm[i] != 0.0:
            rdf[i] =  h[i] / norm[i]
        else:
            rdf[i] = 0.0

        #norm = V / (2.0 * math.pi * (n**2) * binw**3  * ((i-.5)**2)) 

    plt.bar(bins[0:bins.size-1],rdf,width=binw)
    plt.xlim([0.0,rmax])
    return r

if __name__ == '__main__':
    """ If run as a script, run the code I used to get
        the first version working.
    """
    n = 1000
    x,y,z = 10.,10.,10.

    plt.figure()
    r = grid3d(n,(10,10,10),(0,0,0))
    rij = rdfplot(r,(x,y,z),rmax=10)
    
    plt.figure()
    r = fcc3d(n,(10,10,10),(0,0,0))
    rij = rdfplot(r,(x,y,z),rmax=10)
    
    plt.figure()
    r = random3d(n,0,x,0,y,0,z)
    rij = rdfplot(r,(x,y,z),rmax=10)



