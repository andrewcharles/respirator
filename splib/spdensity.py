""" A set of functions for computing smoothed particle summation
    densities.
"""
import sys
sys.path.append("../pyticles")
import spkernel
import numpy as np

def sum_pdensity(r,m):
    """ This computes summation density.
        pairs is an array giving the i,j
        indices of interacting pairs
        r is the particle positions 
        m is the particle masses
        We make a neighbour list and calculate the kernel
        in house.
    """
    #for i in range(n):
    print 'Not implemented'
    #    nlist = neighbours.Nlist(n,n)
    #    nlist.build(r[i],r,h,periodic=True,xmax=xmax,ymax=ymax)
    #    w = numpy.zeros(len(nlist.nrefs))
    #    for j in range(len(nlist.nrefs)): 
    #        w[j] = spkernel.kernel(nlist.ds[j],nlist.dr[j],h[j],'lucy')[0]
    #    rho[i] = spdensity.summation_density(w,m)
    #return rho

def summation_density(w,m):
    """ 
        w - a numpy array of kernel values
        m - a numpy array of masses
        function returns the value of the density
    """
    n = w.size
    rho = 0.0
    for i in range(n):
        rho += m[i] * w[i]
    return rho

def summation_density_gradient(dw,m):
    """ 
        w - a numpy array of kernel gradient values
        m - a numpy array of masses
        function returns the value of the density
    """
    n = dw.size
    drho = 0.0
    for i in range(n):
        drho -= m[i] * dw[i]
    return drho

def summation_density_laplacian(dw,m,drho,rho):
    """ 
        w - a numpy array of kernel gradient values
        m - a numpy array of masses
        function returns the value of the density
    """
    n = dw.size
    ddrho = 0.0
    for i in range(n):
        ddrho -= (m[i]/rho[i]) * dw[i] * drho[i]
    return ddrho

def sum_density(m,h,nl,rij,drij,adke=False):
    zdist = np.zeros(drij[0,:].shape)
    n = h.size
    ni = nl.shape[0]
    rho = np.zeros(n)
    
    for i in range(n):
        rho[i] = spkernel.lucy_kernel(0.0,zdist,h[i])[0] * m[i]

    # calc the kernels, velocities and densities
    for k in range(ni):
        i = nl[k,0]
        j = nl[k,1]
        wij = spkernel.lucy_kernel(rij[k],drij[k,:],h[i])[0]
        rho[i] += wij * m[j]
        rho[j] += wij * m[i]

    if adke:
        # We are using adaptive density kernel estimation
        # the density calculated above was just a pilot
        # the smoothing length above is the reference length
        h = adapt_smoothing_exp(rho,h)
        #h = adapt_smoothing_lap(rho,h)
        #h = adapt_smoothing(rho,h)

        for i in range(n):
            rho[i] = spkernel.lucy_kernel(0.0,zdist,h[i])[0] * m[i]
        
        for k in range(ni):
            i = nl[k,0]
            j = nl[k,1]
            wij = spkernel.lucy_kernel(rij[k],drij[k,:],h[i])[0]
            rho[i] += wij * m[j]
            rho[j] += wij * m[i]

    return rho,h

def adapt_smoothing(rho,h):
    H = h.mean()
    KSC = 1.0
    SENS = -0.5
    rhoav = np.mean(rho)
    return( H * KSC * ((rho/rhoav)**SENS) )

def adapt_smoothing_exp(rho,h):
    H = h.mean()
    KSC = 1.0
    SENS = -0.5
    n = rho.size
    rhoav = np.exp(np.mean(np.log(rho)))
    return( H * KSC * ((rho/rhoav)**SENS) )

def adapt_smoothing_grad(mgrho,h):
    H = h.mean()
    KSC = 1.0
    JSC = 4.5
    SENS = -0.5
    n = mgrho.size
    #laprhoav = np.exp(np.mean(np.log(laprho)))
    #laprhoav = np.mean(laprho)
    # deviation from mean curvature
    #return( H * KSC * (( (np.abs(laprho)+1.0)/(np.abs(laprhoav+1.0)))**SENS) )
    # Just based on curvature
    return( H * KSC * ( ( (JSC*mgrho+1.0) )**SENS) )

def sum_density_gradient(m,h,nl,rij,drij):
    n = h.size
    ni = nl.shape[0]
    ni,ndim = drij.shape
    drho = np.zeros([n,ndim])
    
    # calc the kernels and densities
    for k in range(ni):
        i = nl[k,0]
        j = nl[k,1]
        w,dwij = spkernel.lucy_kernel(rij[k],drij[k,:],h[i])[1]
        drho[i] -= dwij * m[j]
        drho[j] += dwij * m[i]

    return drho

def sum_density_laplacian(m,h,nl,rij,drij):
    """ Product rule symmetrised laplacian. """
    n = h.size
    ni = nl.shape[0]
    drho = sum_density_gradient(m,h,nl,rij,drij)
    gdotdrho = np.zeros(n)
    
    # calc the kernels and densities
    for k in range(ni):
        i = nl[k,0]
        j = nl[k,1]
        rhoij = (rho[i]+rho[j])/2.
        mij = ((m[i]+m[j])/2.)
        w,dwij = spkernel.lucy_kernel(rij[k],drij[k,:],h[i])[1]
        gdotdrho[i] -= (mij/rhoij) * (drho(i) - drho(j)) * dwij
        gdotdrho[j] += (mij/rhoij) * (drho(i) - drho(j)) * dwij

    return gdotdrho

