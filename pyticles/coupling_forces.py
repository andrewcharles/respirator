""" 
    A generic force class ,and several force objects that operate on particle
    systems.

    Force -- A generic pairwise force that can be instantiated but
            does nothing.


    Andrew Charles
"""

import numpy
import math
import scipy
import neighbour_list
from forces import Force

# CONSTANTS
CUTOFF_RADIUS = 10
DIM = 2
COLLISION_RADIUS_SQ = 1.0
VACUUM_VISCOSITY = 0.1
spring_k = 1.1
rest_distance =  0.2

class CouplingForce:
    """ A generic pairwise particle force, between two particle
        systems. The assumption is that the force is mutual. 
        
        p1: The first particle set
        p2: The second particle set
        nl: list of interacting pairs (p1,p2)

        The design is wrong here...

    """ 

    def __init__(self,particles1,particles2,nl):
        self.p1 = particles1
        self.p2 = particles2
        self.nl = nl
    
    def apply(self):
        """ Apply the force to all particles in the nlist """
        for k in range(self.nl.nip):
            self.apply_force(k)


class HookesForce(Force):
    def __init__(self,particles1,particles2,nl):
        self.p1 = particles1
        self.p2 = particles2
        self.nl = nl
        self.nl.cutoff_radius_sq = CUTOFF_RADIUS**2 

    def apply_force(self,k):
        """ Takes a particle p and an nlist pair reference k
        """
        # magnitude of force is k(|r-ro|)
        # should calculate this in nlist allocation
        # and use the stored value
        # actually, should calc it in a seperate method so
        # we can get new distances without getting a new
        # nlist completely
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p1 = self.p1
        p2 = self.p2

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 
        fmag = abs (spring_k * ( rdist - rest_distance ) )
            
        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p1.vdot[i,0] += dvx
        p1.vdot[i,1] += dvy
        p2.vdot[j,0] += -dvx
        p2.vdot[j,1] += -dvy


class CollisionForce(Force):
    
    def __init__(self,particles1,particles2,nl):
        self.p1 = particles1
        self.p2 = particles2
        self.nl = nl
        #self.nl.build_nl_brute_force()
        #self.nl.cutoff_radius_sq = COLLISION_RADIUS_SQ 

    def apply_force(self,k):
        """ A hard collision is just an instantanous force. More
            like an impulse maybe. Anyhow, it changes v directly.
            TODO: rollback and roll forward to get rid of sticky balls.
        """

        if self.nl.rij[k]**2 > self.nl.cutoff_radius_sq:
            #print "toofar"
            return
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        p1 = self.p1
        p2 = self.p2
        debug = False 
        # Calc v dot r, magrhosq, and others
        dr = self.nl.drij[k]
        
        drsq = dr[0]**2 + dr[1]**2
        vidotr = p1.v[i,0]*dr[0] + p1.v[i,1]*dr[1]
        vjdotr = (p2.v[j,0]*dr[0] + p2.v[j,1]*dr[1])

        # If the particles are moving away from each other do nothing
        if (vidotr < 0) and (vjdotr >  0):
            return

        # Calculate tangential and normal components
        # of velocity
        vtix = (vidotr/drsq) * dr[0] 
        vtiy = (vidotr/drsq) * dr[1]
        vnix = p1.v[i,0] - vtix
        vniy = p1.v[i,1] - vtiy

        vtjx = (vjdotr/drsq) * dr[0] 
        vtjy = (vjdotr/drsq) * dr[1]
        vnjx = p2.v[j,0] - vtjx
        vnjy = p2.v[j,1] - vtjy

        if debug:
            vmagsqi = p1.v[i,0]**2 + p2.v[i,1]**2
            vmagsqj = p1.v[j,0]**2 + p2.v[j,1]**2
            vmagsqrot = (vtix+vnix)**2 + (vtiy+vniy)**2
            mom_i = math.sqrt(vmagsqi)
            mom_j = math.sqrt(vmagsqj)
            print "Before"
            print mom_i,mom_j,mom_i+mom_j
            #print vmagsqreg,vmagsqrot

        # Transfer tangential component of momentum
        tmp = vtix
        vtix = vtjx * (p2.m[j]/p1.m[i])
        vtjx = tmp * (p1.m[i]/p2.m[j])
        
        tmp = vtiy
        vtiy = vtjy * (p2.m[j]/p1.m[i])
        vtjy = tmp * (p1.m[i]/p2.m[j])

        # Convert back to xy frame
        p1.v[i,0] = vtix + vnix 
        p1.v[i,1] = vtiy + vniy
        p2.v[j,0] = vtjx + vnjx
        p2.v[j,1] = vtjy + vnjy
        
        if debug:
            vmagsqi = p1.v[i,0]**2 + p1.v[i,1]**2
            vmagsqj = p2.v[j,0]**2 + p2.v[j,1]**2
            mom_i = math.sqrt(vmagsqi)
            mom_j = math.sqrt(vmagsqj)
            print "After"
            print mom_i,mom_j,mom_i+mom_j
            exit()


class SpamForce(Force):
    
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl

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
        dvx =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[0]
        dvy =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[1]
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        ### derived properties
        # force accumulator
        # internal energy
        # smoothing length
        # udot accumulator (sphforce(rho,q,p)
        # Q heat flux tensor
        # P pressure tensor ( eos(rho,t) )


class CohesiveSpamForce(Force):
    
    def __init__(self,particles,nl):
        self.p = particles
        self.nl = nl
        #self.nl.cutoff_radius_sq = CUTOFF_RADIUS**2

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        #def spam(p,i,j,dwdx):
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]
        pri = self.p.pco[i]
        prj = self.p.pco[j]
        p = self.p
        dwdx = self.nl.dwij[k,:]
        dvx =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[0]
        dvy =  (pri/p.rho[i]**2 + prj/p.rho[j]**2) * dwdx[1]
        p.vdot[i,0] += dvx
        p.vdot[i,1] += dvy
        p.vdot[j,0] += -dvx
        p.vdot[j,1] += -dvy

        ### derived properties
        # force accumulator
        # internal energy
        # smoothing length
        # udot accumulator (sphforce(rho,q,p)
        # Q heat flux tensor
        # P pressure tensor ( eos(rho,t) )


class Gravity(Force):
    """Controller that attracts other objects with an inverse square force.

       Acceleration of affected particles is computed as 

                dv/dt = (g*M)/r^2

       and directed towards the centre of the attractive domain.

       To do:
        - can we speed up the sqrt and vector ops?
        - look ahead one frame for position?

    """
    def __init__(self,particles1,particles2,nl):
        self.p1 = particles1
        self.p2 = particles2
        self.nl = nl
        self.g = 03.

    def apply_force(self,k):
        """ Calculates spam interaction between two particles.
            The spam density must have already been calculated.
        """
        i = self.nl.iap[k,0]
        j = self.nl.iap[k,1]

        p1 = self.p1
        p2 = self.p2

        drx =  self.nl.drij[k,0] 
        dry =  self.nl.drij[k,1]
        rdist = self.nl.rij[k]
        rsquared = rdist**2 

        fmag = self.g*p1.m[i]*p2.m[j]/rsquared

        #resolve into components
        dvx = fmag * ( drx ) / rdist
        dvy = fmag * ( dry ) / rdist
        p1.vdot[i,0] += dvx
        p1.vdot[i,1] += dvy
        p2.vdot[j,0] += -dvx
        p2.vdot[j,1] += -dvy
        

