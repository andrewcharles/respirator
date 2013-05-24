""" A wrapper around the fortran sphforce3d routine.
    Computes some smooth particle properties first.
    np.set_printoptions(precision=5,suppress=True)

"""

import sys
import fkernel
import fsplib
import feos
import numpy as np
import sphforce3d
from time import time
from forces import PairwiseForce

class SpamComplete(PairwiseForce):
    """ Compute all the smooth particle properties and forces
        in a single function that uses the fortran routines.

        Parameters:
            adash -- van der Waals attraction (sets in feos module)
            bdash -- van der Waals repulsion (sets in feos module)
            kbdash -- van der Waals k (sets in feos module)
            sigma -- repulsive core size (fortran module variable)
            rcoef -- repulsive core strength (fortran module variable)
            cgrad -- density gradient coefficient
            eta -- shear viscosity
            zeta -- bulk viscosity
            thermalk -- thermal conductivity
            kernel_type -- 1: Gauss, 2: Lucy, 3: Debrun

    """

    def __init__(self,particles,neighbour_list,
        adash=2.0,
        bdash=0.5,
        kbdash=1.0,
        sigma=0.0,
        rcoef=0.0,
        cgrad=1.0,
        eta=0.0,
        zeta=0.0,
        thermalk=5.0,
        kernel_type=2,
        art_viscosity=1,
        inner_cutoff=0.2,
        cutoff=5.0):

        # Call the generic Force init method
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

        # Set the parameters in this module and imported modules
        feos.eos.adash = adash
        feos.eos.bdash = bdash
        feos.eos.kbdash = kbdash
        feos.eos.eos_set_vdw_params(adash,bdash,kbdash)
        self.adash = adash
        self.bdash = bdash
        self.kbdash = kbdash
        self.art_viscosity = art_viscosity
        sphforce3d.sphforce3d.set_vdw_params(adash,bdash,kbdash)
        sphforce3d.sphforce3d.core_sigma = sigma
        sphforce3d.sphforce3d.core_rcoef = rcoef
        sphforce3d.sphforce3d.inner_cutoff = rcoef
        sphforce3d.sphforce3d.cgrad = cgrad
        sphforce3d.sphforce3d.thermalk = thermalk
        sphforce3d.sphforce3d.print_vdw_params()
        self.eta = eta
        self.zeta = zeta
        self.kernel_type = 2 

    def apply(self):
        """ Calculates spam interaction between all particles.

            2010-07-31
            Performance may be impacted by array copies? Not likely. 
            Currently the only copy in the fortran calls is the ilist.  Having
            to copy the particle variables is a tiny fraction of the cost of
            this function call.
            For 8^3 particles we have (roughly)
            0.05 s for kernels
            0.01 s for density
            0.17 s for force

            print is 10 microseconds

            So the bulk of the expense is in the fortran 3d force subroutine.
        """
        p = self.p
        nl = self.nl

        #sphforce3d.sphforce3d.print_vdw_params()
        # Size parameters
        n = p.n
        d = p.dim
        ni = nl.nip

        t = time()
        # Copy / reference particle properties
        x = np.asfortranarray(p.r[0:n,:])   
        v = np.asfortranarray(p.v[0:n,:])
        T = np.asfortranarray(p.t[0:n])
        mass = np.asfortranarray(p.m[0:n])
        rho = np.zeros((n),order='F')
        rho_lr = np.zeros((n),order='F')
        u = np.asfortranarray(p.u[0:n])
        sml = np.asfortranarray(p.h[0:n])
        sml_lr = np.asfortranarray(p.hlr[0:n])
        
        # Copy / reference Neighbourly properties
        ilist = np.asfortranarray((nl.iap[0:ni,:]+1))
        dv =  np.asfortranarray(nl.dv[0:ni,:])
        rij = np.asfortranarray(nl.rij[0:ni])
        drij = np.asfortranarray(nl.drij[0:ni,:])
        
        # New allocations
        a = np.zeros([n,d],order='F')
        #x = np.zeros([n,d],order='F')
        jq = np.zeros([n,d],order='F')
        grad_rho = np.zeros((n,d),order='F')
        grad_v = np.zeros((n,d,d),order='F')
        grad_rho_lr = np.zeros((n,d),order='F')
        #print time()-t,'init'

        # This next section is an application of the Fortran library
        # we have created.
        # Velocity diff
        if ni == 0:
            return

        fsplib.splib.calc_dv(dv,ilist,v)
        t = time()
        
        # Kernels and kernel gradients
        w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml,
            self.kernel_type)
        w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij,drij,ilist,
            sml_lr,self.kernel_type)
        
        # TIMING
        #print time()-t,'kernels'

        t = time()
        # Density summation
        fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,
            self.kernel_type)
        fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass,w_lr,
            dwdx_lr,self.kernel_type)
        #print time()-t,'density'

        # Pressure tensor components
        p_rev = np.zeros((n,d,d),order='F')
        p_rev_lr = np.zeros((n,d,d),order='F')
        pi_irr = np.zeros((n,d,d),order='F')

        # Speed of sound
        c = np.ones([n],dtype=float,order='F')
        
        # Constants
        eta = np.zeros([n],order='F')
        zeta = np.zeros([n],order='F')
        eta[:] = self.eta 
        zeta[:] = self.zeta

        dedt = np.zeros([n],dtype=float,order='F')

        # Internal energy is the dynamical quantity that is
        # updated, so we set the temperature accordingly
        # There is an inverse method with the signature
        # feos.eos.calc_vdw_energy(u,T,rho)
        feos.eos.calc_vdw_temp(u,T,rho)
        T [T < 0.0] = 0.0
        avisc_on = self.art_viscosity
        t = time()
        # Call the force subroutine
        # This is where the majority of time is spent
        sphforce3d.sphforce3d.calc_sphforce3d( 
            ilist,x,v,a,  
            p_rev,p_rev_lr,pi_irr,          
            grad_rho,grad_v,               
            u,dedt,mass,rho,rho_lr,T,jq,     
            c,eta,zeta,               
#            self.adash,self.bdash,self.kbdash,
            dv,rij,drij,  
            sml,sml_lr,w,dwdx,dwdx_lr,avisc_on)
        #print time()-t,'force'
        p.timing['SPAM time'] = time() - t
        
        #feos.eos.calc_vdw_temp(u,T,rho)

        # Resend data to python object
        p.grad_rho[0:n,:] = grad_rho_lr[0:n,:]
        p.rho[0:n] = rho[0:n]
        p.rho_lr[0:n] = rho_lr[0:n]
        p.udot[0:n] = dedt[0:n]
        nl.wij[0:ni] = w[0:ni]
        nl.dwij[0:ni,:] = dwdx[0:ni,:]
        nl.wij_lr[0:ni] = w_lr[0:ni]
        nl.dwij_lr[0:ni,:] = dwdx_lr[0:ni,:]
        p.P[0:n] = p_rev + p_rev_lr + pi_irr
        p.p[0:n] = p_rev[:,0,0]
        # should this not be the normal of [0,0],[1,1],[2,2] ?
        p.pco[0:n] = p_rev_lr[:,0,0]
        p.vdot[0:n,:] = a[0:n,:]
        p.t[0:n] = T[0:n]
        p.jq[0:n,:] = jq[0:n,:]


class SpamHeat(PairwiseForce):
    """ Like spam_complete but only the heat flux is computed,
        there is no kinetic force.

        Parameters:
            adash -- van der Waals attraction (sets in feos module)
            bdash -- van der Waals repulsion (sets in feos module)
            kbdash -- van der Waals k (sets in feos module)
            sigma -- repulsive core size (fortran module variable)
            rcoef -- repulsive core strength (fortran module variable)
            cgrad -- density gradient coefficient
            eta -- shear viscosity
            zeta -- bulk viscosity
            thermalk -- thermal conductivity
            kernel_type -- 1: Gauss, 2: Lucy, 3: Debrun

    """

    def __init__(self,particles,neighbour_list,
        adash=2.0,
        bdash=0.5,
        kbdash=1.0,
        sigma=0.0,
        rcoef=0.0,
        cgrad=1.0,
        eta=0.0,
        zeta=0.0,
        thermalk=5.0,
        kernel_type=2,
        art_viscosity=True,
        cutoff=5.0):

        # Call the generic Force init method
        PairwiseForce.__init__(self,particles,neighbour_list,cutoff=cutoff)

        self.art_viscosity = art_viscosity
        # Set the parameters in this module and imported modules
        feos.eos.adash = adash
        feos.eos.bdash = bdash
        feos.eos.kbdash = kbdash
        self.adash = adash
        self.bdash = bdash
        self.kbdash = kbdash
        sphforce3d.sphforce3d.set_vdw_params(adash,bdash,kbdash)
        sphforce3d.sphforce3d.core_sigma = sigma
        sphforce3d.sphforce3d.core_rcoef = rcoef
        sphforce3d.sphforce3d.cgrad = cgrad
        sphforce3d.sphforce3d.thermalk = thermalk
        self.eta = eta
        self.zeta = zeta
        self.kernel_type = 2 

    def apply(self):
        """ Calculates spam interaction between all particles.
        """
        p = self.p
        nl = self.nl

        # Size parameters
        n = p.n
        d = p.dim
        ni = nl.nip

        t = time()
        # Copy / reference particle properties
        x = np.asfortranarray(p.r[0:n,:])   
        v = np.asfortranarray(p.v[0:n,:])
        T = np.asfortranarray(p.t[0:n])
        mass = np.asfortranarray(p.m[0:n])
        rho = np.zeros((n),order='F')
        rho_lr = np.zeros((n),order='F')
        u = np.asfortranarray(p.u[0:n])
        sml = np.asfortranarray(p.h[0:n])
        sml_lr = np.asfortranarray(p.hlr[0:n])
        
        # Copy / reference Neighbourly properties
        ilist = np.asfortranarray((nl.iap[0:ni,:]+1))
        dv =  np.asfortranarray(nl.dv[0:ni,:])
        rij = np.asfortranarray(nl.rij[0:ni])
        drij = np.asfortranarray(nl.drij[0:ni,:])
        
        # New allocations
        a = np.zeros([n,d],order='F')
        #x = np.zeros([n,d],order='F')
        jq = np.zeros([n,d],order='F')
        grad_rho = np.zeros((n,d),order='F')
        grad_v = np.zeros((n,d,d),order='F')
        grad_rho_lr = np.zeros((n,d),order='F')
        #print time()-t,'init'

        # This next section is an application of the Fortran library
        # we have created.
        # Velocity diff
        if ni == 0:
            return

        fsplib.splib.calc_dv(dv,ilist,v)
        t = time()
        
        # Kernels and kernel gradients
        w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml,
            self.kernel_type)
        w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij,drij,ilist,
            sml_lr,self.kernel_type)
        
        # TIMING
        #print time()-t,'kernels'

        t = time()
        # Density summation
        fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,
            self.kernel_type)
        fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass,w_lr,
            dwdx_lr,self.kernel_type)
        #print time()-t,'density'

        # Pressure tensor components
        p_rev = np.zeros((n,d,d),order='F')
        p_rev_lr = np.zeros((n,d,d),order='F')
        pi_irr = np.zeros((n,d,d),order='F')

        # Speed of sound
        c = np.ones([n],dtype=float,order='F')
        
        # Constants
        eta = np.zeros([n],order='F')
        zeta = np.zeros([n],order='F')
        eta[:] = self.eta 
        zeta[:] = self.zeta

        dedt = np.zeros([n],dtype=float,order='F')

        # Internal energy is the dynamical quantity that is
        # updated, so we set the temperature accordingly
        # There is an inverse method with the signature
        # feos.eos.calc_vdw_energy(u,T,rho)
        feos.eos.calc_vdw_temp(u,T,rho)
        T [T < 0.0] = 0.0

        t = time()
        # Call the force subroutine
        # This is where the majority of time is spent

        sphforce3d.sphforce3d.calc_sphforce3d( 
            ilist,x,v,a,  
            p_rev,p_rev_lr,pi_irr,          
            grad_rho,grad_v,               
            u,dedt,mass,rho,rho_lr,T,jq,     
            c,eta,zeta,               
            dv,rij,drij,  
            sml,sml_lr,w,dwdx,dwdx_lr,self.art_viscosity)
        
        p.timing['SPAM time'] = time() - t

        # Resend data to python object
        p.grad_rho[0:n,:] = grad_rho_lr[0:n,:]
        p.rho[0:n] = rho[0:n]
        p.rho_lr[0:n] = rho_lr[0:n]
        p.udot[0:n] = dedt[0:n]
        nl.wij[0:ni] = w[0:ni]
        nl.dwij[0:ni,:] = dwdx[0:ni,:]
        nl.wij_lr[0:ni] = w_lr[0:ni]
        nl.dwij_lr[0:ni,:] = dwdx_lr[0:ni,:]
        #p.P[0:n] = p_rev + p_rev_lr + pi_irr
        #p.p[0:n] = p_rev[:,0,0]
        #p.pco[0:n] = p_rev_lr[:,0,0]
        p.vdot[0:n,:] = 0.0
        #p.t[0:n] = T[0:n]
        p.jq[0:n,:] = jq[0:n,:]

