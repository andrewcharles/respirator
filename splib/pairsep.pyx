""" 
    The most expensive part of a neighbour list is the computation of pair
    separations.

    Andrew Charles 2013

"""

import numpy as np
cimport numpy as cnp
import math
import scipy
import neighbour_list
import particles
cimport cython

# CONSTANTS

ctypedef cnp.float_t DTYPE_t
ctypedef cnp.int_t DTYPE_int
DTYPE = np.float

cdef extern from "math.h":
    float sinf(float theta)
    float sqrtf(float a)
    float powf(float a, float b)

@cython.boundscheck(False)
def pairsep(nl):
    """ Computes the seperations between a set of points.
        r[n,3] = array of positions
        iap[nk,2] = integer array of interacting pairs
        dr[nk,3] = returned array of pair separations
        ds[nk] = returned array of pair separation distances

    """
        
    cdef cnp.ndarray[cnp.int_t,ndim=2,mode='c'] _pairs 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _r
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _v
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dr 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dv
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _ds
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
    cdef float drx, dry, drz
    cdef unsigned int i,j,k

    n = nl.particle.n
    nk = nl.nip
    dim = nl.particle.dim

    _pairs = nl.iap.astype(np.int)
    _r = nl.particle.r.astype(DTYPE)
    _v = nl.particle.v.astype(DTYPE)
    
    _dr = np.zeros((nk,dim))
    _ds = np.zeros(nk)
    _dv =  np.zeros((nk,dim))
    _rsq =  np.zeros(nk)

    # Which is faster?
    #_dr = nl.rij.astype(np.float)
    #_ds = nl.rij.astype(np.float)
    #_dv = nl.dv.astype(np.float)
    #_rsq = nl.rsq.astype(np.float)

    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        drz = _r[j,2] - _r[i,2]
        _rsq[k] = powf(drx,2) + powf(dry,2) + powf(drz,2)
        _ds[k] = sqrtf(_rsq[k])
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dr[k,2] = drz
        _dv[k,0] = _v[j,0] - _v[i,0]
        _dv[k,1] = _v[j,1] - _v[i,1]
        _dv[k,2] = _v[j,2] - _v[i,2]

    nl.drij[0:nk,:] = _dr[0:nk,:]
    nl.rij[0:nk] = _ds[0:nk]
    nl.rsq[0:nk] = _rsq[0:nk]
    nl.dv[0:nk,:] = _dv[0:nk,:]

@cython.boundscheck(False)
def pairsep2d(nl):
    """ Computes the seperations between a set of points.
        r[n,2] = array of positions
        iap[nk,2] = integer array of interacting pairs
        dr[nk,2] = returned array of pair separations
        ds[nk] = returned array of pair separation distances

    """
        
    cdef cnp.ndarray[cnp.int_t,ndim=2,mode='c'] _pairs 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _r
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _v
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dr 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dv
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _ds
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
    cdef float drx, dry, drz
    cdef unsigned int i,j,k

    n = nl.particle.n
    nk = nl.nip
    dim = nl.particle.dim

    _pairs = nl.iap.astype(np.int)
    _r = nl.particle.r.astype(DTYPE)
    _v = nl.particle.v.astype(DTYPE)
    
    _dr = np.zeros((nk,dim))
    _ds = np.zeros(nk)
    _dv =  np.zeros((nk,dim))
    _rsq =  np.zeros(nk)

    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        _rsq[k] = powf(drx,2) + powf(dry,2)
        _ds[k] = sqrtf(_rsq[k])
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dv[k,0] = _v[j,0] - _v[i,0]
        _dv[k,1] = _v[j,1] - _v[i,1]

    nl.drij[0:nk,:] = _dr[0:nk,:]
    nl.rij[0:nk] = _ds[0:nk]
    nl.rsq[0:nk] = _rsq[0:nk]
    nl.dv[0:nk,:] = _dv[0:nk,:]

def pairsep_minimage(nl):
    """ Computes the seperations between a set of points.
        r[n,3] = array of positions
        iap[nk,2] = integer array of interacting pairs
        dr[nk,3] = returned array of pair separations
        ds[nk] = returned array of pair separation distances

    """
        
    cdef cnp.ndarray[cnp.int_t,ndim=2,mode='c'] _pairs 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _r
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _v
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dr 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dv
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _ds
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
    cdef float drx, dry, drz, xmax, ymax, zmax
    cdef unsigned int i,j,k,n,nk,dim

    n = nl.particle.n
    nk = nl.nip
    dim = nl.particle.dim

    #_pairs = nl.iap.astype(np.int)
    #_r = nl.particle.r.astype(DTYPE)
    #_v = nl.particle.v.astype(DTYPE)

    _pairs = nl.iap[0:nk,:]
    _r = nl.particle.r[0:n,:]
    _v = nl.particle.v[0:n,:]

    xmax = nl.particle.box.xmax
    ymax = nl.particle.box.ymax
    zmax = nl.particle.box.zmax
    
    _dr = np.zeros((nk,dim))
    _ds = np.zeros(nk)
    _dv =  np.zeros((nk,dim))
    _rsq =  np.zeros(nk)

    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        drz = _r[j,2] - _r[i,2]
        if (drx > xmax/2.):
            drx = drx - xmax
        if (dry > ymax/2.):
            dry = dry - ymax 	
        if (drz > zmax/2.):
            drz = drz - zmax
        if (drx < -xmax/2.):
            drx = drx + xmax
        if (dry < -ymax/2.):
            dry = dry + ymax 	
        if (drz < -zmax/2.):
            drz = drz + zmax
        _rsq[k] = powf(drx,2) + powf(dry,2) + powf(drz,2)
        _ds[k] = sqrtf(_rsq[k])
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dr[k,2] = drz
        _dv[k,0] = _v[j,0] - _v[i,0]
        _dv[k,1] = _v[j,1] - _v[i,1]
        _dv[k,2] = _v[j,2] - _v[i,2]

    nl.drij[0:nk,:] = _dr[0:nk,:]
    nl.rij[0:nk] = _ds[0:nk]
    nl.rsq[0:nk] = _rsq[0:nk]
    nl.dv[0:nk,:] = _dv[0:nk,:]


@cython.boundscheck(False)
def fast_pairsepr(cnp.ndarray[DTYPE_int,ndim=2] iap, 
    cnp.ndarray[DTYPE_t,ndim=2] r
    ):
    """ Computes the seperations between a set of points.
        Just distance, not anything else
        r[n,3] = array of positions
        iap[nk,2] = integer array of interacting pairs
        dr[nk,3] = returned array of pair separations
        ds[nk] = returned array of pair separation distances
    """
        
    cdef cnp.ndarray[DTYPE_int,ndim=2,mode='c'] _pairs 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _r
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dr 
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _ds
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
    cdef float drx, dry, drz
    cdef unsigned int i,j,k,n,nk,dim

    n = r.shape[0]
    nk = iap.shape[0]
    dim = r.shape[1]

    _pairs = iap #iap.astype(np.int)
    _r = r #nl.particle.r.astype(DTYPE)
    
    _dr = np.zeros((nk,dim),dtype=DTYPE)
    _ds = np.zeros((nk),dtype=DTYPE)
    _rsq =  np.zeros((nk),dtype=DTYPE)

    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        drz = _r[j,2] - _r[i,2]
        _rsq[k] = powf(drx,2) + powf(dry,2) + powf(drz,2)
        _ds[k] = sqrtf(_rsq[k])
        _ds[k] = (_rsq[k])
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dr[k,2] = drz
    
    return _dr,_ds

@cython.boundscheck(False)
def fast_pairsep(cnp.ndarray[DTYPE_int,ndim=2] iap, 
    cnp.ndarray[DTYPE_t,ndim=2] r,
    cnp.ndarray[DTYPE_t,ndim=2] v
    ):
    """ Computes the seperations between a set of points.
        r[n,3] = array of positions
        iap[nk,2] = integer array of interacting pairs
        dr[nk,3] = returned array of pair separations
        ds[nk] = returned array of pair separation distances
    """
        
    cdef cnp.ndarray[DTYPE_int,ndim=2,mode='c'] _pairs 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _r
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _v
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dr 
    cdef cnp.ndarray[DTYPE_t,ndim=2,mode='c'] _dv
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _ds
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
    cdef float drx, dry, drz
    cdef unsigned int i,j,k,n,nk,dim

    n = r.shape[0]
    nk = iap.shape[0]
    dim = r.shape[1]

    _pairs = iap #iap.astype(np.int)
    _r = r #nl.particle.r.astype(DTYPE)
    _v = v #nl.particle.v.astype(DTYPE)
    
    _dr = np.zeros((nk,dim),dtype=DTYPE)
    _ds = np.zeros((nk),dtype=DTYPE)
    _dv =  np.zeros((nk,dim),dtype=DTYPE)
    _rsq =  np.zeros((nk),dtype=DTYPE)

    # Which is faster?
    #_dr = nl.rij.astype(np.float)
    #_ds = nl.rij.astype(np.float)
    #_dv = nl.dv.astype(np.float)
    #_rsq = nl.rsq.astype(np.float)

    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        drz = _r[j,2] - _r[i,2]
        _rsq[k] = powf(drx,2) + powf(dry,2) + powf(drz,2)
        _ds[k] = sqrtf(_rsq[k])
        _ds[k] = (_rsq[k])
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dr[k,2] = drz
        _dv[k,0] = _v[j,0] - _v[i,0]
        _dv[k,1] = _v[j,1] - _v[i,1]
        _dv[k,2] = _v[j,2] - _v[i,2]

    return _dr,_ds,_dv

@cython.boundscheck(False)
def fast_pairsep_minimage(
    cnp.ndarray[DTYPE_int,ndim=2] _pairs, 
    cnp.ndarray[DTYPE_t,ndim=2] _r,
    cnp.ndarray[DTYPE_t,ndim=2] _v,
    cnp.ndarray[DTYPE_t,ndim=1] _box,
    cnp.ndarray[DTYPE_t,ndim=2] _dr,
    cnp.ndarray[DTYPE_t,ndim=1] _ds,
    cnp.ndarray[DTYPE_t,ndim=2] _dv

    ):
    """ Computes the seperations between a set of points.
        r[n,3] = array of positions
        iap[nk,2] = integer array of interacting pairs
        v = array of velocities
        box = max of box dims
        
        dr[nk,3] = returned array of pair separations
        ds[nk] = returned array of pair separation distances

        This is about 20% faster than the one above
        NB drsq not returned from this one

    """
        
    cdef float drx, dry, drz, rsq
    cdef float xmax, ymax, zmax
    cdef unsigned int i,j,k,n,nk,dim,ex
    
    ex = 2

    n = _r.shape[0]
    nk = _pairs.shape[0]
    dim = _r.shape[1]

    xmax = _box[0]
    ymax = _box[1]
    zmax = _box[2]
    
    for k in xrange(nk):
        i = _pairs[k,0]
        j = _pairs[k,1]
        drx = _r[j,0] - _r[i,0]
        dry = _r[j,1] - _r[i,1]
        drz = _r[j,2] - _r[i,2]
        if (drx > xmax/2.):
            drx = drx - xmax
        if (dry > ymax/2.):
            dry = dry - ymax 	
        if (drz > zmax/2.):
            drz = drz - zmax
        if (drx < -xmax/2.):
            drx = drx + xmax
        if (dry < -ymax/2.):
            dry = dry + ymax 	
        if (drz < -zmax/2.):
            drz = drz + zmax
        _rsq = powf(drx,ex) + powf(dry,ex) + powf(drz,ex)
        _ds[k] = sqrtf(_rsq)
        _dr[k,0] = drx
        _dr[k,1] = dry
        _dr[k,2] = drz
        _dv[k,0] = _v[j,0] - _v[i,0]
        _dv[k,1] = _v[j,1] - _v[i,1]
        _dv[k,2] = _v[j,2] - _v[i,2]

    return

