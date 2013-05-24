""" 
    A generic force class and several force objects that operate on particle
    systems.

    Andrew Charles 2013

"""

import numpy as np
cimport numpy as np
import math
import scipy
import neighbour_list
import particles
import f_properties as properties
#cimport _particles

# CONSTANTS
CUTOFF_RADIUS = 10
DIM = 2
COLLISION_RADIUS_SQ = 1.0
VACUUM_VISCOSITY = 0.1
spring_k = 1.1
rest_distance =  0.2

# Uses this type for c_numpy arrays
ctypedef np.float_t DTYPE_t
DTYPE = np.float


cdef extern from "math.h":
    float sqrtf(float a)
    float powf(float a, float b)

class Force:
    """ A generic pairwise particle force, between two particle
        systems. The assumption is that the force is mutual. 
        
        p1: The first particle set
        p2: The second particle set
        nl: list of interacting pairs (p1,p2)

    """ 

    def __init__(self,particles1,particles2,nl):
        self.p1 = particles1
        self.p2 = particles2
        self.nl = nl
    
    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            self.apply_force(k)


class CollisionForce(Force):
    """ Not yet optimised to the point that it's faster than the
        fortran code inside a python for loop.
    """

    def __init__(self,particles,nl,cutoff=1.0):
        self.p = particles
        self.nl = nl
        self.cutoff = cutoff

    def apply(self):
        """Iterate over the neighbour list and apply the force to all
        particles.
        """

        p = self.p
        nl = self.nl

        cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap 
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _m
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _v
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dv
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _drij
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rsq
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rij
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _r

        # The underscore means this is to be propagated back to the
        # particle and nlist objects
        _iap = self.nl.iap.astype(np.int)
        _m = p.m.astype(np.float)
        _v = p.v.astype(np.float)
        _r = p.r.astype(np.float)
        _dv = nl.dv.astype(np.float)
        _drij = nl.drij.astype(np.float)
        _rsq = nl.rsq.astype(np.float)
        _rij = nl.rij.astype(np.float)

        cdef unsigned int nip = self.nl.nip
        cdef unsigned int i,j,k,n
        cdef float dvx, dvy, dvz, cut, ro, tmp, drx, dry, drz, drox, droy, droz
        cut = self.cutoff
        n = p.n

        for k in xrange(nip):
            i = _iap[k,0]
            j = _iap[k,1]
            
            #_drsq =  _drij[0],2) + powf(_drij[1],2) + powf(_drij[2],2)
            _drx = _drij[k,0] 
            _dry = _drij[k,1] 
            _drz = _drij[k,2] 
            # Separate
            # todo rollback entire system to the point of collision
            if _rsq[k] < (self.cutoff*self.cutoff):
                ro = (self.cutoff + 0.01 - _rij[k])/2.0
                drox = (_rij[k]/_drij[k,0]) * ro
                droy = (_rij[k]/_drij[k,1]) * ro
                droz = (_rij[k]/_drij[k,2]) * ro
                _r[i,0] -= drox
                _r[i,1] -= droy
                _r[i,2] -= droz
                _r[j,0] += drox
                _r[j,1] += droy
                _r[j,2] += droz
                _drx = _r[j,0] - _r[i,0]
                _dry = _r[j,1] - _r[i,1]
                _drz = _r[j,2] - _r[i,2]
                _rsq[k] = powf(_drx,2) + powf(_dry,2) + powf(_drz,2)
                _rij[k] = sqrtf(_rsq[k])
                p.r[i,0] = _r[i,0]
                p.r[i,1] = _r[i,1]
                p.r[i,2] = _r[i,2]
                p.r[j,0] = _r[j,0]
                p.r[j,1] = _r[j,1]
                p.r[j,2] = _r[j,2]
                nl.rsq[k] = _rsq[k]
                nl.rij[k] = _rij[k]
                nl.drij[k,0] = _drx
                nl.drij[k,1] = _dry
                nl.drij[k,2] = _drz

            # Divergence
            vidotr = _v[i,0]*drx + _v[i,1]*dry + _v[i,2]*drz
            vjdotr = _v[j,0]*drx + _v[j,1]*dry + _v[j,2]*drz

            # Do we not want to take the divergence of the two
            # velocity vectors?

            # If the particles are moving away from each other do nothing
            if (vidotr > 0) and (vjdotr <  0):
                # Calculate tangential and normal components
                # of velocity
                vtix = (vidotr/_rsq[k]) * drx
                vtiy = (vidotr/_rsq[k]) * dry
                vtiz = (vidotr/_rsq[k]) * drz
                vnix = _v[i,0] - vtix
                vniy = _v[i,1] - vtiy
                vniz = _v[i,2] - vtiz

                vtjx = (vjdotr/_rsq[k]) * drx
                vtjy = (vjdotr/_rsq[k]) * dry
                vtjz = (vjdotr/_rsq[k]) * drz
                vnjx = _v[j,0] - vtjx
                vnjy = _v[j,1] - vtjy
                vnjz = _v[j,2] - vtjz

                # Transfer tangential component of momentum
                tmp = vtix
                vtix = vtjx * (_m[j]/_m[i])
                vtjx = tmp * (_m[i]/_m[j])
                
                tmp = vtiy
                vtiy = vtjy * (_m[j]/_m[i])
                vtjy = tmp * (_m[i]/_m[j])

                tmp = vtiz
                vtiz = vtjz * (_m[j]/_m[i])
                vtjz = tmp * (_m[i]/_m[j])

                # Convert back to xy frame
                p.v[i,0] = vtix + vnix 
                p.v[i,1] = vtiy + vniy
                p.v[i,2] = vtiz + vniz
                p.v[j,0] = vtjx + vnjx
                p.v[j,1] = vtjy + vnjy
                p.v[j,2] = vtjz + vnjz

        #p.r[0:n,:] = _r[0:n,:]
        #p.v[0:n,:] = _v[:,:]


class SpamForce(Force):

    #cdef public _particles.SmoothParticleSystem p

    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl

    def apply(self):
        """Iterate over the neighbour list and apply the force to all
        particles.
        """
        properties.spam_properties(self.p,self.nl \
            ,self.p.h[0:self.p.n],self.p.hlr[0:self.p.n])

        p = self.p
        nl = self.nl

        cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap 
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dwdx
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _p
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rho
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _m
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _udot
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _vdot
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dv

        _iap = self.nl.iap.astype(np.int)
        _dwdx = nl.dwij.astype(np.float)
        _p = p.p.astype(np.float)
        _rho = p.rho.astype(np.float)
        _m = p.m.astype(np.float)
        _udot = p.udot.astype(np.float)
        _vdot = p.vdot.astype(np.float)
        _dv = nl.dv.astype(np.float)

        cdef unsigned int nip = self.nl.nip
        cdef unsigned int i,j,k,n
        cdef float dvx, dvy, dvz, du, ps, ax, ay, az

        n = p.n

        for k in xrange(nip):
            i = _iap[k,0]
            j = _iap[k,1]
            ps = (_p[i]/(_rho[i]**2) + _p[j]/(_rho[j]**2))
            ax = ps * _dwdx[k,0] 
            ay = ps * _dwdx[k,1]
            az = ps * _dwdx[k,2]
            du = 0.5 * (ax * _dv[k,0] + ay * _dv[k,1] + az * _dv[k,2])
            _vdot[i,0] += ax
            _vdot[i,1] += ay
            _vdot[i,2] += az
            _vdot[j,0] += -ax
            _vdot[j,1] += -ay
            _vdot[j,2] += -az
            _udot[j] +=  du * _m[i]
            _udot[i] +=  du * _m[j]

        p.udot[0:n] = _udot[:]
        p.rdot[0:n,:] = p.v[0:n,:]
        p.vdot[0:n,:] = _vdot[:,:]


    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        pass


class CohesiveSpamForce(Force):
    
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl

    def apply(self):
        """Iterate over the neighbour list and apply the force to all
        particles.
        """

        p = self.p
        nl = self.nl

        cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap 
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dwdx
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _p
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rho
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _m
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _udot
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _vdot
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dv

        _iap = self.nl.iap.astype(np.int)
        _dwdx = nl.dwij_lr.astype(np.float)
        _p = p.pco.astype(np.float)
        _rho = p.rho_lr.astype(np.float)
        _m = p.m.astype(np.float)
        _udot = p.udot.astype(np.float)
        _vdot = p.vdot.astype(np.float)
        _dv = nl.dv.astype(np.float)

        cdef unsigned int nip = self.nl.nip
        cdef unsigned int i,j,k,n
        cdef float dvx, dvy, dvz, du, ps, ax, ay, az

        n = p.n

        for k in xrange(nip):
            i = _iap[k,0]
            j = _iap[k,1]
            ps = (_p[i]/(_rho[i]**2) + _p[j]/(_rho[j]**2))
            ax = ps * _dwdx[k,0] 
            ay = ps * _dwdx[k,1]
            az = ps * _dwdx[k,2]
            du = 0.5 * (ax * _dv[k,0] + ay * _dv[k,1] + az * _dv[k,2])
            _vdot[i,0] += ax
            _vdot[i,1] += ay
            _vdot[i,2] += az
            _vdot[j,0] -= ax
            _vdot[j,1] -= ay
            _vdot[j,2] -= az
            _udot[j] +=  du * _m[i]
            _udot[i] +=  du * _m[j]

        p.udot[0:n] = _udot[:]
        p.rdot[0:n,:] = p.v[0:n,:]
        p.vdot[0:n,:] = _vdot[:,:]


class SpamConduction(Force):
    """ Heat conduction using the full heat flux vector
        which is probably not as efficient as JM's.
    """

    #cdef public _particles.SmoothParticleSystem p

    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl

    def apply(self):
        """Iterate over the neighbour list and apply the force to all
        particles.
        """

        p = self.p
        nl = self.nl

        cdef np.ndarray[np.int_t,ndim=2,mode='c'] _iap 
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _dwdx
        cdef np.ndarray[DTYPE_t,ndim=2,mode='c'] _q
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _rho
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _m
        cdef np.ndarray[DTYPE_t,ndim=1,mode='c'] _udot

        _iap = self.nl.iap.astype(np.int)
        _dwdx = nl.dwij.astype(np.float)
        _q = p.jq.astype(np.float)
        _rho = p.rho.astype(np.float)
        _m = p.m.astype(np.float)
        _udot = p.udot.astype(np.float)

        cdef unsigned int nip = self.nl.nip
        cdef unsigned int i,j,k,n
        cdef float dux, duy, duz, du, ps

        n = p.n

        for k in xrange(nip):
            i = _iap[k,0]
            j = _iap[k,1]
           
            dux = (_q[i,0]/_rho[i]**2 + _q[j,0]/_rho[j]**2) * _dwdx[k,0]
            duy = (_q[i,1]/_rho[i]**2 + _q[j,1]/_rho[j]**2) * _dwdx[k,1]
            duz = (_q[i,2]/_rho[i]**2 + _q[j,2]/_rho[j]**2) * _dwdx[k,2]

            _udot[i] -=  dux * _m[j]
            _udot[i] -=  duy * _m[j]
            _udot[i] -=  duz * _m[j]
            _udot[j] +=  dux * _m[i]
            _udot[j] +=  duy * _m[i]
            _udot[j] +=  duz * _m[i]

        p.udot[0:n] = _udot[:]
