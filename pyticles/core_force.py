""" A wrapper around the fortran core force routine.
    Computes some smooth particle properties first.
    np.set_printoptions(precision=5,suppress=True)

"""

import sys
import numpy as np
from time import time
from forces import PairwiseForce
import core_potential


class CoreForce(PairwiseForce):
    """ This is a simple pairwise core force, see the fortran
        code for the formula

        Parameters:
            sigma -- repulsive core size (fortran module variable)
            rcoef -- repulsive core strength (fortran module variable)
            cutoff

    """

    def __init__(self,particles,neighbour_list,
        sigma=1.0,
        rcoef=1.0,
        cutoff=5.0):
        self.sigma = sigma
        self.rcoef = rcoef
        # Call the generic Force init method
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

    def apply(self):
        """ Call a fortran loop over the nlist for all pairs.
        """
        p = self.p
        nl = self.nl

        n = p.n
        d = p.dim
        ni = nl.nip

        t = time()
        # Copy / reference particle properties
        x = np.asfortranarray(p.r[0:n,:])   
        v = np.asfortranarray(p.v[0:n,:])
        m = np.asfortranarray(p.m[0:n])
        
        # Copy / reference Neighbourly properties
        ilist = np.asfortranarray((nl.iap[0:ni,:]+1))
        rij = np.asfortranarray(nl.rij[0:ni])
        drij = np.asfortranarray(nl.drij[0:ni,:])
        
        # New allocations
        a = np.zeros([n,d],order='F')

        if ni == 0:
            return

        t = time()
        
        # This is where the majority of time is spent
        core_potential.core_potential.calc_coreforce3d( 
            ilist,x,a,m,rij,drij,self.sigma,self.rcoef)
        p.timing['SPAM time'] = time() - t
        # Resend data to python object
        p.vdot[0:n,:] = a[0:n,:]

