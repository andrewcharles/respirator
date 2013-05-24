""" Calculation of sph densities, interparticle kernel functions.

    Relies on the SmoothParticleSystem object.

    See also f_properties.py which is an implementation of the
    same thing using the fsph Fortran modules.

    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.

"""

import sys
#sys.path.append("../vasp")
import spdensity
import spkernel
import math
import numpy as np

from eos import ideal_isothermal, art_water, vdw, vdw_energy, vdw_temp

#calc_pressure = art_water
calc_pressure = vdw

def hamiltonian(p):
    """ Compute the Hamiltonian of a smooth particle system. 
        Untested.
    """
    H = 0.0
    for i in range(p.n):
        H = H + 0.5 * p.m[i] * (p.v[i,0]**2 + p.v[i,1]**2 + p.v[i,2]**2)
        H = H + p.u[i]
    print H
    return H

def spam_properties2d(p,nl,adke=False):
    """ Calculates and assigns:

        1. kernel, kernel gradients
        2. smoothed particle summation densities
        * velocity gradient
        * pressure
        * internal energy

    """
    h = p.h
    hlong = p.hlr
    p.gradv[0:p.n] = 0.0
    n = p.n
    ni = nl.nip

    # calc the kernels, velocities and densities
    for k in range(nl.nip):
        i = nl.iap[k,0]
        j = nl.iap[k,1]
        nl.wij[k], nl.dwij[k] = spkernel.lucy_kernel(nl.rij[k],nl.drij[k,:],p.h[i])

    p.rho[:],p.h[:] = spdensity.sum_density(p.m[0:n],p.h[0:n],nl.iap[0:ni,:]
        ,nl.rij[0:ni],nl.drij[0:ni,:])
    
        #dv = p.v[j,:] - p.v[i,:]
        #for a in range(p.dim):
        #    for b in range(p.dim):
        #        p.gradv[i,a,b] += (p.m[j]/p.rho[i])*dv[a]*nl.dwij[k,b]
        #        p.gradv[j,a,b] += (p.m[i]/p.rho[j])*dv[a]*nl.dwij[k,b]

    for i in range(p.n):
        p.p[i],p.pco[i] = calc_pressure(p.rho[i],p.t[i])

    p.u = vdw_energy(p.rho,p.t)
    p.t = vdw_temp(p.rho,p.u)

def spam_properties(p,nl,adke=False):
    """ Calculates and assigns:

        * kernel values
        * kernel gradient values
        * smoothed particle summation densities
        * velocity gradient
        * pressure
        * internal energy

        adke -- adaptive density kernel estimation
    
    """
    # self contribution to density
    h = p.h
    hlong = p.hlr
    zerokern = spkernel.lucy_kernel(0.0,(0.0,0.0,0.0),h[0])[0]
    p.rho[0:p.n] = zerokern
    p.gradv[0:p.n] = 0.0
    n = p.n
    ni = nl.nip

    # calc the kernels, velocities and densities
    for k in range(nl.nip):
        i = nl.iap[k,0]
        j = nl.iap[k,1]

        nl.wij[k], nl.dwij[k] = spkernel.lucy_kernel(nl.rij[k],nl.drij[k,:],p.h[i])

        p.rho[i] += nl.wij[k] * p.m[j]
        p.rho[j] += nl.wij[k] * p.m[i]
    
        dv = p.v[j,:] - p.v[i,:]

        for a in range(p.dim):
            for b in range(p.dim):
                p.gradv[i,a,b] += (p.m[j]/p.rho[i])*dv[a]*nl.dwij[k,b]
                p.gradv[j,a,b] += (p.m[i]/p.rho[j])*dv[a]*nl.dwij[k,b]

    if adke:
        # We are using adaptive density kernel estimation
        # the density calculated above was just a pilot
        # the smoothing length above is the reference length
        KSC = 1.0
        SENS = 0.5
        rhoav = np.mean(p.rho)
        p.h = H * KSC * ((p.rho/rhoav)**SENS)
       
    for i in range(p.n):
        # todo add some logic to determine whether we have a one or two part
        # pressure
        p.p[i],p.pco[i] = calc_pressure(p.rho[i],p.t[i])

    rho2 = spdensity.sum_density(p.m[0:n],p.h[0:n],nl.iap[0:ni,:]
        ,nl.rij[0:ni],nl.drij[0:ni,:])

    p.u = vdw_energy(p.rho,p.t)
    p.t = vdw_temp(p.rho,p.u)


def spam_properties_ls(p,nl,adke=False):
    """ Calculates and assigns smooth properties for long and
        short smoothing lengths.

        * kernel values
        * kernel gradient values
        * smoothed particle summation densities
        * velocity gradient
        * pressure
        * internal energy
        
        adke -- adaptive density kernel estimation
    
    """
    # self contribution to density
    h = p.h
    hlr = p.hlr
    zerokern = spkernel.lucy_kernel(0.0,(0.0,0.0,0.0),h[0])[0]
    zerokern_lr = spkernel.lucy_kernel(0.0,(0.0,0.0,0.0),hlr[0])[0]
    p.rho[0:p.n] = zerokern
    p.rho_lr[0:p.n] = zerokern_lr
    p.gradv[0:p.n] = 0.0
    n = p.n
    ni = nl.nip

    # calc the kernels, velocities and densities
    for k in range(ni):
        i = nl.iap[k,0]
        j = nl.iap[k,1]

        nl.wij[k], nl.dwij[k] = spkernel.lucy_kernel(nl.rij[k],nl.drij[k,:],
            h[i])

        nl.wij_lr[k], nl.dwij_lr[k] = spkernel.lucy_kernel(nl.rij[k],
            nl.drij[k,:],hlr[i])

        p.rho[i] += nl.wij[k] * p.m[j]
        p.rho[j] += nl.wij[k] * p.m[i]
        
        p.rho_lr[i] += nl.wij_lr[k] * p.m[j]
        p.rho_lr[j] += nl.wij_lr[k] * p.m[i]
    
        dv = p.v[j,:] - p.v[i,:]

    for k in range(ni):
        for a in range(p.dim):
            for b in range(p.dim):
                p.gradv[i,a,b] += (p.m[j]/p.rho_lr[i])*dv[a]*nl.dwij_lr[k,b]
                p.gradv[j,a,b] += (p.m[i]/p.rho_lr[j])*dv[a]*nl.dwij_lr[k,b]

    if adke:
        # We are using adaptive density kernel estimation
        # the density calculated above was just a pilot
        # the smoothing length above is the reference length
        KSC = 1.0
        SENS = 0.5
        rhoav = np.mean(p.rho)
        p.h = H * KSC * ((p.rho/rhoav)**SENS)
       
    for i in range(p.n):
        p.p[i],p.pco[i] = calc_pressure(p.rho[i],p.t[i])

    rho2 = spdensity.sum_density(p.m[0:n],p.h[0:n],nl.iap[0:ni,:]
        ,nl.rij[0:ni],nl.drij[0:ni,:])

    p.u = vdw_energy(p.rho,p.t)
    p.t = vdw_temp(p.rho,p.u)


def test_vdw_energy_and_temp():
   rho = 1.0
   t = 5.0
   u = vdw_energy(rho,t)
   print u
   t2 = vdw_temp(rho,u)
   print t2

if __name__ == '__main__':
    test_vdw_energy_and_temp()



