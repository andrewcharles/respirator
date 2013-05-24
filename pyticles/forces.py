""" 
    A generic force class ,and several force objects that operate on particle
    PairwiseForce -- A generic pairwise force that can be instantiated but
            does nothing.

    HookesForce -- a pairwise force that is linear in particle distance

    Andrew Charles
"""

import numpy as np
import math
import scipy
import neighbour_list
import sys
sys.path.append('/Users/acharles/masters/active/fsph')
from collision import collision
import f_properties as properties
from controller import Controller

class Force(Controller):
    """ This kind of controller mutates particle systems by modifying
        their rates of change. Generic Forces are abstract, we don't
        ever expect to create them unless for testing.
    """

    def apply(self):
        print "Do nothing force"
        return

#class BodyForce(Force):
#    """ A constant one body mechanical force, with a direction and a
#        magnitude.
#    """
#
#    def __init__(self,p,dir=(0,1,0),mag=-1000):
#        """
#        p -- ParticleGroup
#        dir -- direction, either a 2 or 3 tuple depending on the particle
#                system's dimensionality
#
#        """
#        Controller.__init__(self)
#        self.dir = dir
#        self.mag = mag
#
#    def bind_particles(self,p):
#        self.groups.append(p)
##        if p.dim == 2:
#            self.dir = self.dir[0:1]
#
#    def apply(self):
#        """
#        """
#        for p in self.groups:
#            for i in range(p.n):
#                p.vdot[i,:] += (self.mag * numpy.array(self.dir))
#

class BodyForce2d:
    """ A constant one body mechanical force, with a direction and a
        magnitude.
    """
    """ Not a pairwise force. """
    def __init__(self,particles,direction=[0.0,1.0],k=1.0):
        self.p = particles
        self.k = k
        self.direction = direction

    def apply_force(self,i):
        # This is a bit inefficient
        p = self.p
        k = self.k
        p.vdot[i,0] += k * self.direction[0]
        p.vdot[i,1] += k * self.direction[1]

    def apply(self):
        """ Apply the force to all particles in the nlist """
        self.p.vdot[:,0] += self.k * self.direction[0]
        self.p.vdot[:,1] += self.k * self.direction[1]
        #for i in range(self.p.n):
        #    self.apply_force(i)


class BodyForce3d:
    """ Not a pairwise force. """
    def __init__(self,particles,direction=[0.0,1.0,0.0],k=1.0):
        self.p = particles
        self.k = k
        self.direction = direction

    def apply_force(self,i):
        p = self.p
        k = self.k
        p.vdot[i,0] += k * self.direction[0]
        p.vdot[i,1] += k * self.direction[1]
        p.vdot[i,2] += k * self.direction[2]

    def apply(self):
        """ Apply the force to all particles in the nlist """
        for i in range(self.p.n):
            self.apply_force(i)


class PlaneCollisionForce3d(Force):
   
    def __init__(self,particles,plane=(1.0,0.0,0.0),cutoff=5.0):
        print 'Initialising plane collision force'

    def apply_force(self,i):
        """ A hard collision is just an instantanous force.
            rollback and roll forward to get rid of sticky balls.
            Also implemented is a 'pushback' that ensures particles
            are no closer than the cutoff (collision) range.
            This has not been well tested, and there are other
            ways to implement this rollback that are more
            dynamically consistent.
            In fact when used on top of an sph force the pushback seems to 
            increase the instability.

            Start with a plane orthoganal to two axes

            An arbitrary plane will be a little tricky

        """
        debug = False
        p = self.p
        
class PairwiseForce(Force):
    """ A generic pairwise particle force
        systems. The assumption is that the force is mutual. 
        p1: The first particle set
        nl: list of interacting pairs (p1,p2)

        This will become a subclass of Controller

    """ 

    def __init__(self,particles,nl,cutoff=100.0):
        self.p = particles
        self.nl = nl
        self.cutoff = cutoff
        self.cutoffsq = cutoff*cutoff

    def apply_force(self,k):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.nl.p
        print 'Interacting particles i and j.'

    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)

    def apply_sorted(self):
        """ Apply the force to all particles in the sorted nlist """
        """ Todo: check that the list is of type sorted!
        """
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)
            else:
                return


class HookesForce2d(PairwiseForce):
    """
         Magnitude of force is k(|r-ro|)
         Assumes drij has been calculated and stored in nlist
         Takes a particle p and an nlist pair reference k
    """

    def __init__(self,particles,neighbour_list,cutoff=50.0,k=1.0,r0=5.0
        ,damping=0.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.r0 = r0
        self.k = k
        self.damping = damping
        self.cutoffsq = cutoff*cutoff

    def apply_force(self,k):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 
        fmag = (1.0 - self.damping) * abs (self.k * ( rdist - self.r0 ) )
            
        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy


class LennardJonesForce2d(PairwiseForce):
    """
        Magnitude of force is

        -24c [ 2 (sigma^12/r^13) - (sigma^6/r^7) ]
        
        rc -- force magnitude
        sig -- core size
        r -- pair distance

        Assumes drij has been calculated and stored in nlist
        Takes a particle p and an nlist
    """

    def __init__(self,particles,neighbour_list,cutoff=5.0,eps=1.0,sig=1.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.eps = eps
        self.sig = sig

    def apply_force(self,k):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        sig = self.sig
        eps = self.eps
        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsq = rdist**2 
        fmag = -24 * eps * ( 2 * (sig**12 / rdist**13)
                             - (sig**6  / rdist**7) )  
        #resolve into components
        dvx = fmag *  drx / rdist
        dvy = fmag *  dry / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy



class CoreForce2d(PairwiseForce):
    """
        Magnitude of force is

        2 * rc * ( abs(r) / (sig**2) ) * (-(r**2) / sig**2 + 1)**3
    forcemag = rcoef*((2*dx)/sigma_sq) * ( ((-drsq/sigma_sq) + 1)**3 )
        
        rc -- force magnitude
        sig -- core size
        r -- pair distance

        Assumes drij has been calculated and stored in nlist
        Takes a particle p and an nlist
    """

    def __init__(self,particles,neighbour_list,cutoff=5.0,rc=1.0,sig=5.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.rc = rc
        self.sig = sig

    def apply_force(self,k):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        sig = self.sig
        rc = self.rc
        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsq = rdist**2 
        fmag = 2 * rc * ( rdist / (sig**2) )  \
            * ((-rsq / (sig**2) + 1)**3)
        #resolve into components
        dvx = fmag *  drx / rdist
        dvy = fmag *  dry / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

class CoreForce3d(PairwiseForce):
    """
        Magnitude of force is

        2 * rc * ( abs(r) / (sig**2) ) * (-(r**2) / sig**2 + 1)**3
        
        rc -- force magnitude
        sig -- core size
        r -- pair distance

        Assumes drij has been calculated and stored in nlist
        Takes a particle p and an nlist
    """

    def __init__(self,particles,neighbour_list,cutoff=5.0,rc=1.0,sig=5.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.rc = rc
        self.sig = sig

    def apply_force(self,k):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        sig = self.sig
        rc = self.rc
        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        drz =  self.nl.drij[k,2]
        rdist = self.nl.rij[k]
        rsq = rdist**2 
        fmag = 2 * rc * ( rdist / (sig**2) )  \
            * ((-rsq / (sig**2) + 1)**3)
        #resolve into components
        dvx = fmag *  drx / rdist
        dvy = fmag *  dry / rdist
        dvz = fmag *  drz / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[i,2] += dvz
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy
        p.vdot[j,2] += -dvz


class CollisionForce2d(PairwiseForce):
   
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ A hard collision is just an instantanous force.

        """
        debug = False
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        
        dr = self.nl.drij[k]
        drsq = dr[0]**2 + dr[1]**2
      
        # Divergence
        vidotr = p.v[i,0]*dr[0] + p.v[i,1]*dr[1]
        vjdotr = (p.v[j,0]*dr[0] + p.v[j,1]*dr[1])
            
        # Do we not want to take the divergence of the two
        # velocity vectors?

        # If the particles are moving away from each other do nothing
        if (vidotr < 0) and (vjdotr >  0):
            return

        # Calculate tangential and normal components
        # of velocity
        vtix = (vidotr/drsq) * dr[0] 
        vtiy = (vidotr/drsq) * dr[1]
        vnix = p.v[i,0] - vtix
        vniy = p.v[i,1] - vtiy

        vtjx = (vjdotr/drsq) * dr[0] 
        vtjy = (vjdotr/drsq) * dr[1]
        vnjx = p.v[j,0] - vtjx
        vnjy = p.v[j,1] - vtjy

        if debug:
            vmagsqi = p.v[i,0]**2 + p.v[i,1]**2
            vmagsqj = p.v[j,0]**2 + p.v[j,1]**2
            vmagsqrot = (vtix+vnix)**2 + (vtiy+vniy)**2
            mom_i = math.sqrt(vmagsqi)
            mom_j = math.sqrt(vmagsqj)
            print "Before"
            print mom_i,mom_j,mom_i+mom_j

        # Transfer tangential component of momentum
        tmp = vtix
        vtix = vtjx * (p.m[j]/p.m[i])
        vtjx = tmp * (p.m[i]/p.m[j])
        
        tmp = vtiy
        vtiy = vtjy * (p.m[j]/p.m[i])
        vtjy = tmp * (p.m[i]/p.m[j])

        # Convert back to xy frame
        p.v[i,0] = vtix + vnix 
        p.v[i,1] = vtiy + vniy
        p.v[j,0] = vtjx + vnjx
        p.v[j,1] = vtjy + vnjy

        if debug:
            vmagsqi = p.v[i,0]**2 + p.v[i,1]**2
            vmagsqj = p.v[j,0]**2 + p.v[j,1]**2
            mom_i = math.sqrt(vmagsqi)
            mom_j = math.sqrt(vmagsqj)
            print "After"
            print mom_i,mom_j,mom_i+mom_j
            exit()


class CollisionForce3d(PairwiseForce):
   
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ A hard collision is just an instantanous force.
            rollback and roll forward to get rid of sticky balls.
            Also implemented is a 'pushback' that ensures particles
            are no closer than the cutoff (collision) range.
            This has not been well tested, and there are other
            ways to implement this rollback that are more
            dynamically consistent.
            In fact when used on top of an sph force the pushback seems to 
            increase the instability.
        """
        debug = False
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        
        dr = self.nl.drij[k]
        drsq = dr[0]**2 + dr[1]**2 + dr[2]**2
      
        # Divergence
        vidotr = p.v[i,0]*dr[0] + p.v[i,1]*dr[1] + p.v[i,2]*dr[2]
        vjdotr = p.v[j,0]*dr[0] + p.v[j,1]*dr[1] + p.v[j,2]*dr[2]

        # If the particles are moving away from each other do nothing
        if (vidotr < 0) and (vjdotr >  0):
            return

        # Calculate tangential and normal components
        # of velocity
        vtix = (vidotr/drsq) * dr[0] 
        vtiy = (vidotr/drsq) * dr[1]
        vtiz = (vidotr/drsq) * dr[2]
        vnix = p.v[i,0] - vtix
        vniy = p.v[i,1] - vtiy
        vniz = p.v[i,2] - vtiz

        vtjx = (vjdotr/drsq) * dr[0] 
        vtjy = (vjdotr/drsq) * dr[1]
        vtjz = (vjdotr/drsq) * dr[2]
        vnjx = p.v[j,0] - vtjx
        vnjy = p.v[j,1] - vtjy
        vnjz = p.v[j,2] - vtjz

        # Transfer tangential component of momentum
        tmp = vtix
        vtix = vtjx * (p.m[j]/p.m[i])
        vtjx = tmp * (p.m[i]/p.m[j])
        
        tmp = vtiy
        vtiy = vtjy * (p.m[j]/p.m[i])
        vtjy = tmp * (p.m[i]/p.m[j])

        tmp = vtiz
        vtiz = vtjz * (p.m[j]/p.m[i])
        vtjz = tmp * (p.m[i]/p.m[j])

        # Convert back to xy frame
        p.v[i,0] = vtix + vnix 
        p.v[i,1] = vtiy + vniy
        p.v[i,2] = vtiz + vniz
        p.v[j,0] = vtjx + vnjx
        p.v[j,1] = vtjy + vnjy
        p.v[j,2] = vtjz + vnjz

        if drsq < (self.cutoff*self.cutoff):
            print 'separating',self.nl.rij[k]
            ro = (self.cutoff + 0.01 - self.nl.rij[k])/2.0
            dro = (dr/self.nl.rij[k]) * ro
            p.r[i,0] -= dro[0]
            p.r[i,1] -= dro[1]
            p.r[i,2] -= dro[2]
            p.r[j,0] += dro[0]
            p.r[j,1] += dro[1]
            p.r[j,2] += dro[2]
            self.nl.drij[k] = p.r[j,:] - p.r[i,:]
            self.nl.rij[k] = np.linalg.norm(self.nl.drij[k,:])
            print 'separated',self.nl.rij[k],ro,np.linalg.norm(dro)
            if self.nl.rij[k] < self.cutoff:
                print 'wtf'
                print ro

def SpamForceBasic2d(r,p,nl):
    """ Purely functional definition of the force. """
    """ r -- positions
        p -- pressures
        rho -- density
        nl -- neighbour list
    """
    pass

class SpamForce2d(PairwiseForce):
    
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
            The pressure too...
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.p[i]
        prj = self.p.p[j]
        dwdx = self.nl.dwij[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        dvx = ps * dwdx[0] 
        dvy = ps * dwdx[1] 

        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        du = 0.5 * (dvx * self.nl.dv[k,0] + dvy * self.nl.dv[k,1]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]


class CohesiveSpamForce2d(PairwiseForce):
   
    def __init__(self,particles,neighbour_list,cutoff=10.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        """p = self.p
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        pri = p.pco[i]
        prj = p.pco[j]
        dwdx = self.nl.dwij[k,:]
        dvx =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[0]
        dvy =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[1]
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy"""

        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.p[i]
        prj = self.p.p[j]
        dwdx = self.nl.dwij[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        dvx = ps * dwdx[0] 
        dvy = ps * dwdx[1] 

        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        du = 0.5 * (dvx * self.nl.dv[k,0] + dvy * self.nl.dv[k,1]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]


class SpamForce(PairwiseForce):
    
    def __init__(self,particles,neighbour_list,cutoff=5.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply(self):
        """ Apply the force to all particles in the nlist """
        #properties.spam_properties(self.p,self.nl \
        #    ,self.p.h[0:self.p.n],self.p.hlr[0:self.p.n])
        #print self.p.rho
        #1/0
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """

        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p

        pri = self.p.p[i]
        prj = self.p.p[j]
        dwdx = self.nl.dwij[k,:]
        dv = self.nl.dv[k,:]
       
        ps = (pri/p.rho[i]**2 + prj/p.rho[j]**2)

        ax = ps * dwdx[0] 
        ay = ps * dwdx[1] 
        az = ps * dwdx[2]

        # a is not the acceleration. It is the term in the sph
        # momentum equation sans the mass

        p.vdot[i,0] += p.m[j] * ax
        p.vdot[i,1] += p.m[j] * ay
        p.vdot[i,2] += p.m[j] * az
        p.vdot[j,0] -= p.m[i] * ax
        p.vdot[j,1] -= p.m[i] * ay
        p.vdot[j,2] -= p.m[i] * az

        du = 0.5 * (ax * dv[0] + ay * dv[1] + az * dv[2]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]

class CohesiveSpamForce(PairwiseForce):
    """ Moderately ridiculous that this is a seperate function, there is
        only one line difference..."""
   
    def __init__(self,particles,neighbour_list,cutoff=10.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p = self.p
        pri = self.p.pco[i]
        prj = self.p.pco[j]
        dwdx = self.nl.dwij_lr[k,:]
        dv = self.nl.dv[k,:]
       
        ps = (pri/p.rho_lr[i]**2 + prj/p.rho_lr[j]**2)

        ax = ps * dwdx[0] 
        ay = ps * dwdx[1] 
        az = ps * dwdx[2]

        p.vdot[i,0] += p.m[j] * ax
        p.vdot[i,1] += p.m[j] * ay
        p.vdot[i,2] += p.m[j] * az
        p.vdot[j,0] -= p.m[i] * ax
        p.vdot[j,1] -= p.m[i] * ay
        p.vdot[j,2] -= p.m[i] * az

        du = 0.5 * (ax * dv[0] + ay * dv[1] + az * dv[2]) 
        p.udot[i] += du * p.m[j]
        p.udot[j] += du * p.m[i]

class Gravity2d(PairwiseForce):
    """Controller that attracts other objects with an inverse square force.

       Acceleration of affected particles is computed as 

                dv/dt = (g*M)/r^2

       and directed towards the centre of the attractive domain.

       To do:
        - can we speed up the sqrt and vector ops?
        - look ahead one frame for position?

    """
    def __init__(self,particles,neighbour_list,cutoff=100.0,g=1.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.g = g

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]

        p = self.p

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 

        fmag = self.g*p.m[i]*p.m[j]/rsquared

        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy
       
class Gravity3d(PairwiseForce):
    """Controller that attracts other objects with an inverse square force.

       Acceleration of affected particles is computed as 

                dv/dt = (g*M)/r^2

       and directed towards the centre of the attractive domain.

    """
    def __init__(self,particles,neighbour_list,cutoff=100.0,g=1.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.g = g

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]

        p = self.p

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        drz =  self.nl.drij[k,2]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 

        fmag = self.g*p.m[i]*p.m[j]/rsquared

        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        dvz = fmag * ( drz ) / rdist
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[i,2] += dvz
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy
        p.vdot[j,2] += -dvz


class FortranCollisionForce(PairwiseForce):
    """ A wrapper for the fortran collide3d subroutine.
    """

    def __init__(self,particles,neighbour_list,cutoff=1.0):
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)
        self.nl = neighbour_list

    # todo put this loop in c or fortran
    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            if self.nl.rij[k]**2 <= self.cutoffsq:
                self.apply_force(k)

    def apply_force(self,k):
        """ A hard collision is just an instantanous force.
        """
        nl = self.nl
        i = nl.iap[k,0]
        j = nl.iap[k,1]
        p = self.p
        collision.collide3d(p.v[i],p.v[j],p.m[i],p.m[j],nl.drij[k],nl.rsq[k])


