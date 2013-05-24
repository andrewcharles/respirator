""" Fast three dimensional lightweight neighbour list. """
""" Only depends on numpy.

    Computes lists of neighbouring pairs of points for arbitrary 
    applications.

    Andrew Charles 2013

"""

import numpy as np
cimport numpy as cnp
from math import *
import math
cimport cython

ctypedef cnp.float_t DTYPE_t
ctypedef cnp.int_t DTYPE_int

cdef extern from "math.h":
    float sqrtf(float a)
    float powf(float a, float b)

DTYPE = np.float
DTYPE_i = np.int

NDIM = 3
MAXNEB = 20

def minimum_image(
    dr,
    xmax,
    ymax,zmax):
    """ Applies the minimum image convention to the distance
        between two particles. 
    """

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


def build_list(
    cnp.ndarray[DTYPE_t,ndim=2] r
    ):
    """ Build a brute force neighbour list given an array of positions. """
    cdef unsigned int i,j,k,npair,n
    cdef cnp.ndarray[cnp.int_t,ndim=2] iap
    n = r.shape[0]
    npair = (n * (n-1)) / 2
    iap = np.zeros([npair,2],dtype=DTYPE_i)
    i,j,k = 0,0,0
    nip = 0
    for i in range(n):
        for j in range(i+1,n):
            iap[k,0] = i
            iap[k,1] = j
            k += 1
    return iap

def separations(
    cnp.ndarray[DTYPE_t,ndim=2] r,
    cnp.ndarray[DTYPE_int,ndim=2] nl
    ):
    """ Computes the distance between pairs in the list and stores
        the result in the array rij, indexed by the same k that
        indexes the interacting pairs array nl.
    """
    cdef cnp.ndarray[DTYPE_t,ndim=2] rij
    cdef unsigned int i,j,k,npair
    npair = nl.shape[0]
    rij = np.zeros([npair,3],dtype=DTYPE)
    for k in xrange(npair):
        i = nl[k,0]
        j = nl[k,1]
        rij[k,0] = r[j,0] - r[i,0]
        rij[k,1] = r[j,1] - r[i,1]
        rij[k,2] = r[j,2] - r[i,2]
    return rij

@cython.boundscheck(False)
def compress(
    cnp.ndarray[DTYPE_t,ndim=2] rij,
    cnp.ndarray[DTYPE_int,ndim=2] nl,
    float cutoffsq
    ):
    cdef float rsquared, drx, dry
    cdef unsigned int i,j,k,npair,ex,q
    npair = nl.shape[0]
    ex = 2
    q = 0
    for k in xrange(npair):
        i = nl[k,0]
        j = nl[k,1]
        drx = rij[k,0]
        dry = rij[k,1]
        drz = rij[k,2]
        rsquared = powf(drx,ex) + powf(dry,ex) + powf(drz,ex)
        if (rsquared < cutoffsq):
            rij[q,0] = drx
            rij[q,1] = dry
            rij[q,2] = drz
            nl[q,0] = i
            nl[q,1] = j
            q += 1
    nl = nl[0:q,:]
    rij = rij[0:q,:]
    return rij,nl

def distances(
    cnp.ndarray[DTYPE_t,ndim=2] rij
    ):
    """ np is the number of pairs
    """
    cdef cnp.ndarray[DTYPE_t,ndim=1] dist
    cdef float drx, dry
    cdef unsigned int ex, npair
    
    ex = 2
    npair = rij.shape[0]
    dist = np.zeros([npair],dtype=DTYPE)
    for k in range(npair):
        drx = rij[k,0]
        dry = rij[k,1]
        drz = rij[k,2]
        rsquared = powf(drx,ex) + powf(dry,ex) + powf(drz,ex)
        dist[k] = sqrtf(rsquared)
    return dist



