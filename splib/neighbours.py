""" Light weight neighbour list
    Only depends on numpy.

    Not optimised, but very convenient.

    Computes lists of neighbouring pairs of points for arbitrary 
    applications.

"""

import numpy as np
import math
DIM = 2
MAXNEB = 20

def minimum_image(dr,xmax,ymax):
    """ Applies the minimum image convention to the distance
        between two particles. 
    """
    # don't use the general shape pbs
    # just the 2d rectangle

    pxmax = xmax/2.0
    pymax = ymax/2.0

    if (dr[0] > pxmax ):
        dr[0] = dr[0] - xmax
    if (dr[0] < -pxmax ):
        dr[0] = xmax + dr[0]
    if (dr[1] > pymax ):
        dr[1] =  dr[1] - ymax 
    if (dr[1] < -pymax ):
        dr[1] = ymax  + dr[1]

def min_image(dr,xmax,ymax):
    """ Applies the minimum image convention to the
        pair displacement matrix.
    """
    pxmax = xmax/2.0
    pymax = ymax/2.0

    gxdx = (dr[:,0] > pxmax )
    dr[gxdx,0] = dr[gxdx,0] - xmax
    lxdx = (dr[:,0] < -pxmax )
    dr[lxdx,0] = xmax + dr[lxdx,0]
    gydx = (dr[:,1] > pymax )
    dr[gydx,1] =  dr[gydx,1] - ymax 
    lydx = (dr[:,1] < -pymax )
    dr[lydx,1] = ymax  + dr[lydx,1]

def minimum_image3d(dr,xmax,ymax,zmax):
    """ Applies the minimum image convention to the distance
        between two particles. 
    """
    # don't use the general shape pbs
    # just the 2d rectangle

    pxmax = xmax/2.0
    pymax = ymax/2.0
    pzmax = zmax/2.0

    if (dr[0] > pxmax ):
        dr[0] = dr[0] - xmax
    if (dr[0] < -pxmax ):
        dr[0] = xmax + dr[0]
    if (dr[1] > pymax ):
        dr[1] =  dr[1] - ymax 
    if (dr[1] < -pymax ):
        dr[1] = ymax  + dr[1]
    if (dr[2] > pzmax ):
        dr[2] =  dr[2] - zmax 
    if (dr[2] < -pzmax ):
        dr[2] = zmax  + dr[2]


class Nlist:
    """ This is a neighbour list for an individual point.
        
        Given a point, find the neighbouring points. What neighbouring means is
        defined by the is_neighbour function.

    """

    def __init__(self,nmax,cutoff):
        """ position of the point
            nrefs: indices of neighbours
            r: position of its neighbours
            dr: displacement to each neighbour
            ds: distance to each neighbour
            h: smoothing length of each point
        """
        self.nrefs = []
        self.r = np.zeros((nmax,DIM))
        self.dr = np.zeros((nmax,DIM))
        self.ds = np.zeros(nmax)
        self.a = np.zeros(nmax)
        self.rho = np.zeros(nmax)
        self.m = np.zeros(nmax)
        self.h = np.zeros(nmax)
        self.cutoff = cutoff

    def is_neighbour(self,ds,cutoff):
        """ Checks if points seperated by dr neighbours
        """
        if ds < cutoff: 
            return True
        else:
            return False

    def assign_properties(self,mpts,rhopts,apts,hpts):
        """ Assigns properties to the neighbours of a point.
            Designed to work with build such that you build,
            and then run assign to scatter the densities and
            other properties to the appropriate members.

            nrefs: the indices of the neighbours
            mpts,rhopts,apts are arrays holding the values
        """

        k = 0
        for i in self.nrefs:
            self.h[k] = hpts[int(i)]
            self.m[k] = mpts[int(i)]
            self.rho[k] = rhopts[int(i)]
            self.a[k] = apts[int(i)]
            k = k + 1

        self.h = self.h[0:k]
        self.a = self.a[0:k]
        self.rho = self.rho[0:k]
        self.m = self.m[0:k]

    def build(self,x,rpts,hpts,periodic=False,xmax=0.0,ymax=0.0):
        """ x - the point for which we are finding neighbours
            rpts - np array of all points to consider
            h - smoothing length of each point
                
            uses the smoothing length as the interaction range

            periodic - if this is true, use xmax and ymax to apply a
            rectangular minimum image

        """

        # for each point in rtps
        # check the neighbourness criteria
        k = 0
        for i in range(len(rpts)):
            self.dr[k,0] = rpts[i,0] - x[0]
            self.dr[k,1] = rpts[i,1] - x[1]
            if periodic:
                minimum_image(self.dr[k,:],xmax,ymax)
            self.ds[k] = math.sqrt(self.dr[k,0]**2 + self.dr[k,1]**2)
            
            if self.is_neighbour(self.ds[k],hpts[i]):
                self.nrefs.append(i)
                self.r[k] = rpts[i,:]
                k+=1
            #if k >= MAXNEB:
            #    break
        #truncate
        self.r = self.r[0:k,:]
        self.dr = self.dr[0:k,:]
        self.ds = self.ds[0:k]

def build_list(r):
    """ Build a brute force neighbour list given an array of positions. """
    n = r.shape[0]
    npair = (n * (n-1)) / 2
    iap = np.zeros([npair,2])
    i,j,k = 0,0,0
    nip = 0
    for i in range(n):
        for j in range(i+1,n):
            iap[k,0] = i
            iap[k,1] = j
            k += 1
    return iap

def separations(r,nl):
    """ Computes the distance between pairs in the list and stores
        the result in the array rij, indexed by the same k that
        indexes the interacting pairs array nl.
    """
    npair = nl.shape[0]
    ndim = r.shape[1]
    rij = np.zeros([npair,ndim])
    if ndim == 2:
        for k in range(npair):
            i = nl[k,0]
            j = nl[k,1]
            rij[k,0] = r[j,0] - r[i,0]
            rij[k,1] = r[j,1] - r[i,1]
    elif ndim == 3:
        for k in range(npair):
            i = nl[k,0]
            j = nl[k,1]
            rij[k,0] = r[j,0] - r[i,0]
            rij[k,1] = r[j,1] - r[i,1]
            rij[k,2] = r[j,2] - r[i,2]
    return rij

def compress(rij,nl,cutoffsq = 16.0):
    npair = nl.shape[0]
    ndim = rij.shape[1]
    q = 0
    if ndim == 3:
        for k in range(npair):
            i = nl[k,0]
            j = nl[k,1]
            drx = rij[k,0]
            dry = rij[k,1]
            drz = rij[k,2]
            rsquared = drx**2 + dry**2 + drz**2
            if (rsquared < cutoffsq):
                rij[q,0] = drx
                rij[q,1] = dry
                rij[q,2] = drz
                nl[q,0] = i
                nl[q,1] = j
                q += 1
        nl = nl[0:q,:]
        rij = rij[0:q,:]
    elif ndim == 2:
        for k in range(npair):
            i = nl[k,0]
            j = nl[k,1]
            drx = rij[k,0]
            dry = rij[k,1]
            rsquared = drx**2 + dry**2
            if (rsquared < cutoffsq):
                rij[q,0] = drx
                rij[q,1] = dry
                nl[q,0] = i
                nl[q,1] = j
                q += 1
        nl = nl[0:q,:]
        rij = rij[0:q,:]
    else:
        print 'I crashed'
    return rij,nl

def distances(rij):
    """ np is the number of pairs
    """
    npair = rij.shape[0]
    ndim = rij.shape[1]
    dist = np.zeros([npair])
    if ndim == 3:
        for k in range(npair):
            drx = rij[k,0]
            dry = rij[k,1]
            drz = rij[k,2]
            rsquared = drx**2 + dry**2 + drz**2
            dist[k] = np.sqrt(rsquared)
    elif ndim == 2:
        for k in range(npair):
            drx = rij[k,0]
            dry = rij[k,1]
            rsquared = drx**2 + dry**2
            dist[k] = np.sqrt(rsquared)
    return dist

class nlist_light:
    """ A substitute object for working with Force. """
    def __init__(self,iap,drij,rij):
        self.iap = iap
        self.drij = drij
        self.rij = rij
        self.nip = iap.shape[0]

if __name__ == '__main__':
    n = 8
    side = (2,2,2)
    spacing = 1
    import configuration
    r = configuration.grid3d(n,side,(0,0,0)
        ,spacing=spacing)
    nl = build_list(r)
    rij = separations(r,nl)

