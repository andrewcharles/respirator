""" Wrapper around Fortran SPH code to compute SPH properties """
import sys
sys.path.append('/Users/acharles/masters/active/fsph')
import fkernel
import fsplib
import feos
import numpy as np
import properties

ADASH = 2.0
BDASH = 0.5
KBDASH = 1.0
thermalk = 1.0

def set_vdw_props(adash,bdash,kbdash):
    """ Sets the Fortran module attributes. """
    feos.eos.eos_set_vdw_params(adash,bdash,kbdash)

def print_vdw_props():
    print feos.eos.adash,feos.eos.bdash,feos.eos.kbdash

def energy_from_temp(p,nl,hs,hl):
    """ Set the internal energy from the temperature.
        This is just terrible. Don't do it like this. """
    n = p.n
    d = p.dim
    ni = nl.nip
    x = np.asfortranarray(p.r[0:n,:])   
    v = np.asfortranarray(p.v[0:n,:])
    T = np.asfortranarray(p.t[0:n])
    mass = np.asfortranarray(p.m[0:n])
    rho = np.zeros((n),order='F')
    rho_lr = np.zeros((n),order='F')
    U = np.asfortranarray(p.u[0:n])
    sml = np.asfortranarray(p.h[0:n])
    sml_lr = np.asfortranarray(p.hlr[0:n])
    
    # Copy / reference Neighbourly properties
    ilist = np.asfortranarray((nl.iap[0:ni,:]+1))
    dv =  np.asfortranarray(nl.dv[0:ni,:])
    rij = np.asfortranarray(nl.rij[0:ni])
    drij = np.asfortranarray(nl.drij[0:ni,:])
    
    grad_rho = np.zeros((n,d),order='F')
    grad_rho_lr = np.zeros((n,d),order='F')

    if ni == 0:
        # set values to isolated particle values
        # and issue a warning
        print 'Warning! No interactions!'
        return

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni,:] \
        ,ilist,sml,2)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni] \
        ,ilist,sml_lr,2)
   
    # Density summation
    fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,2)
    fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass  \
        ,w_lr,dwdx_lr,2)

    feos.eos.calc_vdw_energy(U,T,rho)
    p.rho[0:n] = rho[0:n]
    p.rho_lr[0:n] = rho_lr[0:n]
    p.u[0:n] = U[0:n]

def temp_from_energy(p,nl,hs,hl):
    """ Set the temperature from the internal energy """
    n = p.n
    d = p.dim
    ni = nl.nip
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
    
    grad_rho = np.zeros((n,d),order='F')
    grad_rho_lr = np.zeros((n,d),order='F')

    if ni == 0:
        # set values to isolated particle values
        # and issue a warning
        print 'Warning! No interactions!'
        return

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni,:] \
        ,ilist,sml,2)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni] \
        ,ilist,sml_lr,2)

    # Not that here dwdx is \Nabla_j W_ij
   
    # Density summation
    # The density summation is expecting \Nabla_j W_ij
    fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,2)
    fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass,w_lr,dwdx_lr,2)

    feos.eos.calc_vdw_temp(u,T,rho)

    p.rho[0:n] = rho[0:n]
    p.rho_lr[0:n] = rho_lr[0:n]
    p.t[0:n] = T[0:n]

def spam_properties(p,nl,hs,hl):
    """ Calculates and assigns:

        * kernel values
        * kernel gradient values
        * smoothed particle summation densities
        * velocity gradient

        Particle arrays will not be in fortran order.
        See if we can just do copies ($$$$$).
        Despite the copies it's still loads faster than straight python
        and seems to have an advantage over cython.

        nl.iap -- integer list of interactions
        nl.dv -- velocity difference

        p.rho -- particle density
        p.gradv -- particle velocity gradient
       
        I'm writing this to work with the long and short smoothing length

        Doing:
        * Heat flux

        To Do:
        * Long range density
        * Full short range reversible and irreversible tensors
        * Full long range reversible and irreversible tensors

    """

    n = p.n
    d = p.dim
    ni = nl.nip

    # If an array is input/output, it should have the array order set to fortran
    # v = np.reshape(p.v[0:n,:].copy(),(n,d))#,order='F')

    v = np.reshape(p.v[0:n,:].copy(),(n,d))
    dv =  np.reshape(nl.dv[0:ni,:].copy(),(ni,d))
    nlist = nl.iap[0:ni,:].copy()+1
    
    #nlist = nlist.transpose()
    #dv = dv.transpose()
    #v = v.transpose()
    #nlist = np.array(nlist[:,:].transpose(),order='F')

    T = p.t[0:n].reshape((n))#,order='F')
    sml = np.reshape(p.h[0:n].copy(),(n,1))#,order='F')
    sml_lr = np.reshape(p.hlr[0:n].copy(),(n,1))#,order='F')

    rij = np.reshape(nl.rij[0:ni].copy(),(ni,1))#,order='F')
    drij = np.reshape(nl.drij[0:ni,:].copy(),(ni,d))#,order='F')
   
    w = np.zeros((ni),order='F') 
    dwdx = np.zeros((ni,d),order='F') 
   
    w_lr = np.zeros((ni),order='F') 
    dwdx_lr = np.zeros((ni,d),order='F') 

    rho = np.zeros((n))#,order='F')
    rho_lr = np.zeros((n))#,order='F')
    #u = np.reshape(p.u[0:n].copy(),(n,1))#,order='F')
    u = p.u[0:n].reshape((n))#,order='F')
    gradv = np.reshape(p.gradv[0:n,:,:].copy(),(n,d,d))#,order='F')
    jq = np.reshape(p.jq[0:n,:].copy(),(n,d))#,order='F')
    # Listen up this is important.
    # order='F' means the array can be passed as in.out
    # if order is not F and you do an in,out pass
    # you get no output and no warnings?
    grad_rho = np.zeros((n,d),order='F')
    grad_rho_lr = np.zeros((n,d),order='F')
    mass = np.reshape(p.m[0:n].copy(),(n,1))

    if ni == 0:
        # set values to isolated particle values
        # and issue a warning
        print 'Warning! No interactions!'
        return

    #print nlist.flags

    # Velocity diff
    #dv = sphlib.sphlib.calc_dv(dv,nlist,v)
    #print dv.flags
    #print v.flags
    fsplib.splib.calc_dv(dv,nlist,v)

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni,:] \
        ,nlist,sml,2)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni] \
        ,nlist,sml_lr,2)

    #zerokern = spkernel.lucy_kernel(0.0,(0.0,0.0,0.0),h[0])[0]
    #p.rho[0:p.n] = zerokern

    # Density summation
    fkernel.kernel.density_sum(rho,grad_rho,nlist,sml,mass,w,dwdx,2)
    fkernel.kernel.density_sum(rho_lr,grad_rho_lr,nlist,sml_lr,mass,w_lr,dwdx_lr,2)

    # Grad v
    grad_v = np.zeros((n,d,d),order='F')
    fsplib.splib.calc_grad_v(grad_v,nlist,-dwdx,dv,mass,rho)

    phc = np.reshape(p.p[0:n].copy(),(n),order='F')
    pco = np.reshape(p.pco[0:n].copy(),(n),order='F')
 
    # we consider that the internal energy is primary
    feos.eos.calc_vdw_temp(u,T,rho)
    #feos.eos.eos_print_vdw_params()
    feos.eos.calc_vdw_energy(u,T,rho)
    T [T < 0.0] = 0.0
    feos.eos.calc_vdw_hc_pressure(phc,rho,u)
    feos.eos.calc_vdw_cohesive_pressure(pco,rho_lr,u)
    feos.eos.calc_vdw_temp(u,T,rho)

    # Capillary pressure
     
    # Heat Flux
    # No idea why this needs to be order F, and can't be the copy of the
    # particle's jq.
    jq = np.zeros((n,d),order='F')
    
    # subroutine calc_heat_flux_3d(q,ilist,rho,m,tmp,dwdx_jij,thermalk,n,ni)
    fsplib.splib.calc_heat_flux_3d(jq,nlist,rho,mass,T,-dwdx,thermalk)

    # Python implementation
    #phc = vdw_hc(rho,T)
    #pco = vdw_co(rho,T)
    # viscosity

    # Resend data to python object
    p.rho[0:n] = rho[0:n]
    p.rho_lr[0:n] = rho_lr[0:n]
    p.grad_rho[0:n,:] = grad_rho[0:n,:]
    p.grad_rho_lr[0:n,:] = grad_rho_lr[0:n,:]
    nl.wij[0:ni] = w[0:ni]
    nl.dwij[0:ni,:] = dwdx[0:ni,:]
    nl.wij_lr[0:ni] = w_lr[0:ni]
    nl.dwij_lr[0:ni,:] = dwdx_lr[0:ni,:]
    p.p[0:n] = phc[0:n]
    p.pco[0:n] = pco[0:n]
    p.t[0:n] = T[0:n]
    p.jq[0:n,:] = jq[0:n,:]

def density(p,nl,hs,hl):
    """ Set the particle density. """
    n = p.n
    d = p.dim
    ni = nl.nip
    x = np.asfortranarray(p.r[0:n,:])   
    mass = np.asfortranarray(p.m[0:n])
    rho = np.zeros((n),order='F')
    rho_lr = np.zeros((n),order='F')
    sml = np.asfortranarray(p.h[0:n])
    sml_lr = np.asfortranarray(p.hlr[0:n])
    
    # Copy / reference Neighbourly properties
    ilist = np.asfortranarray((nl.iap[0:ni,:]+1))
    rij = np.asfortranarray(nl.rij[0:ni])
    drij = np.asfortranarray(nl.drij[0:ni,:])
    
    grad_rho = np.zeros((n,d),order='F')
    grad_rho_lr = np.zeros((n,d),order='F')

    if ni == 0:
        # set values to isolated particle values
        # and issue a warning
        print 'Warning! No interactions!'
        return

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni,:] \
        ,ilist,sml,2)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij[0:ni],drij[0:ni] \
        ,ilist,sml_lr,2)
   
    # Density summation
    fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,2)
    fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml_lr,mass  \
        ,w_lr,dwdx_lr,2)

    p.rho[0:n] = rho[0:n]
    p.rho_lr[0:n] = rho_lr[0:n]
    






