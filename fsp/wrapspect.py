""" Inspects the python wrapped Fortran.

    This is where we print and plot stuff.

    wraptest.py is where we codify what we learn here
    into a compliance test.

    du/dt  = 1/rho ( -div(J) - P.transpose dot dot grad(v)a

    dv/dt =  1/rho div(P)

    P = phc + pco + pi_os + pi_one

"""

import unittest
import sys
# These paths are now environment variables
#sys.path.append('\Users\acharles\masters\active\pyticles\trunk')
#sys.path.append('\Users\acharles\masters\active\vasp')
#sys.path.append('\home\ac\pyticles')
import numpy as np
import pylab as plt
plt.ion()
import fkernel
import fsplib as splib
from eos import vdw_co,vdw_hc
import feos
from tensor import symmetric_traceless
import spdensity
import spkernel
import renderspam
import interpolate #imports pyticles version of interpolate which has 2d and 3d
import sphforce3d
from numpy import sqrt

def set_force3d_vdw():
    """ In this example we try to set the value of an imported fortran 
        module's variable. While it might look like it works, in fact the
        variable is not accessible to change.

    """
    print sphforce3d.adash
    sphforce3d.adash = 20.0    
    print sphforce3d.adash

def set_force3d_vdw_params():
    """  We created a setter method to set the vdw parameters
    """
    sphforce3d.sphforce3d.print_vdw_params()
    sphforce3d.sphforce3d.set_vdw_params(20.0,0.4,100.0)
    sphforce3d.sphforce3d.print_vdw_params()

def force3d():
    """ Test sphforce3d. """
    """ 
        Four is the magic number for making sure the neighbour list calls
        don't have subtle errors.

    subroutine calc_sphforce3d(ilist,x,v,a,  &
      &   p_rev,p_rev_lr,pi_irr,          &
      &   grad_rho,grad_v,                &
      &   u,dedt,m,rho,rho_lr,temp,q      &
      &   c,eta,zeta,                     &
      &   dv,dx,
      &   sml,sml_long,dwdx,dwdx_long,n,ni)

        1. Compute velocity differences
            -- calc_dv(dv,nlist,v,n,ndim)

        2. Compute kernel weights and gradients

        1. Compute velocity gradient
            -- calc_grad_v(grad_v,ilist,dwdx,dv,m,rho,n,ndim,ni)
        
        2. Compute velocity divergence
            -- calc_div_v(grad_v,div_v,n,ndim)

        3. Compute the symmetric traceless velocity gradient
            -- calc_grad_v_os(grad_v_os,grad_v,div_v,n,ndim)

        4. Compute a symmetric traceless part of the pressure tensor
            -- calc_pi_os(pi_os,grad_v_os,eta,n,ndim)

        5. Compute the isotropic irreversible pressure
            -- calc_pi_one(pi_one,div_v,zeta,n)

        6. Now for the reversible pressure
            -- calc_vdw_hc_pressure(p,rho,u,n)
            -- calc_vdw_cohesive_pressure(pco,rho,u,n)

        7. And voila we have the full pressure tensor
        
    """

    # Data for a system of three particles
    v = np.array([     [ 1.0,  0.0, 0.0],
                       [ 0.0, 3.14, 0.0],
                       [20.0,  0.0, 0.0],
                       [ 0.0,  0.0, 7.0]
    ])


    x = np.array([     [ 1.0, 0.0, 0.0],
                       [ 0.0, 3.1, 0.0],
                       [20.0, 0.0, 0.0],
                       [ 0.0, 7.0, 0.0]
    ])


    n = 4
    ni = 6
    d = 3

    sml = np.array( [ [2.0],[2.0],[2.0],[2.0] ] )
    sml_lr = 2 * sml
    temp = np.ones([n])
    
    nlist = np.zeros((ni,2))
    rij = np.zeros(ni)
    drij = np.zeros([ni,d])
    dv = np.zeros([ni,d])
  
    k = 0
    for i in range(n):
        for j in range(i+1,n):
            print k
            drx = x[j,0] - x[i,0]
            dry = x[j,1] - x[i,1]
            drz = x[j,2] - x[i,2]
            rsquared = drx**2 + dry**2 + drz**2
            drij[k,0] = drx
            drij[k,1] = dry
            drij[k,2] = drz
            rij[k] = sqrt(rsquared)
            dv[k,:] = v[j,:] - v[i,:]
            nlist[k,0] = i
            nlist[k,1] = j
            k += 1
    #ni = k-1
    ilist = nlist + 1
    # Velocity diff
    print 'calc dv'
    splib.splib.calc_dv(dv,nlist+1,v)

    # Kernels and kernel gradients
    w = np.zeros((ni))
    dwdx = np.zeros((ni,3))
    w_lr = np.zeros((ni))
    dwdx_lr = np.zeros((ni,3))
    print 'calc kernels'
    w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml,1)
    w_lr,dwdx_lr = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml*2,1)

    # Density summation
    rho = np.zeros((n))
    rho_lr = np.zeros((n))
    grad_rho = np.zeros((n,d),order='F')
    grad_rho_lr = np.zeros((n,d),order='F')
    a = np.zeros((n,d),order='F')
    mass = np.ones((n))
    print 'calc density'
    fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,1)
    fkernel.kernel.density_sum(rho_lr,grad_rho_lr,ilist,sml*2,mass,w_lr,dwdx_lr,1)

    # Pressure tensor components
    p_rev = np.zeros((n,d,d),order='F')
    p_rev_lr = np.zeros((n,d,d),order='F')
    pi_irr = np.zeros((n,d,d),order='F')
    grad_v = np.zeros((n,d,d),order='F')
    u = np.zeros([n],dtype=float,order='F')
    q = np.zeros([n,d],dtype=float,order='F')
    c = np.ones([n],dtype=float,order='F')
    cgrad = np.ones([n],dtype=float,order='F')
    eta = np.zeros([n],order='F')
    zeta = np.ones([n],order='F')
    dedt = np.zeros([n],dtype=float,order='F')
    
    feos.eos.adash = 25.0
    feos.eos.bdash = 0.5
    feos.eos.kbdash = 1.0
    print 'calc energy'
    feos.eos.calc_vdw_energy(u,temp,rho)

    sphforce3d.sphforce3d.cgrad = 8.35e02
    print 'before'
    # Call the force subroutine
    sphforce3d.sphforce3d.calc_sphforce3d( 
        ilist,x,v,a,  
        p_rev,p_rev_lr,pi_irr,          
        grad_rho,grad_v,               
        u,dedt,mass,rho,rho_lr,temp,q,     
        c,eta,zeta,               
        dv,rij,drij,  
        sml,sml_lr,w,dwdx,dwdx_lr,True)#,n,ni)

    print 'middle'
    np.set_printoptions(precision=5,suppress=True)
    print p_rev
    print a
    print dedt
    print 'x',x
    print 'u',u
    print 'a',a

    x = x + 0.1 * v
    v = v + 0.1 * a
    u = u + 0.1 * dedt

    sphforce3d.sphforce3d.calc_sphforce3d( 
        ilist,x,v,a,  
        p_rev,p_rev_lr,pi_irr,          
        grad_rho,grad_v,               
        u,dedt,mass,rho,rho_lr,temp,q,     
        c,eta,zeta,               
        dv,rij,drij,  
        sml,sml_lr,w,dwdx,dwdx_lr,True)#,n,ni)
    
    print 'after'
    print 'x',x
    print 'u',u
    print 'a',a


def heat_flux():
    """ Let's see if the 3D heat flux routine works
    """
    n = 2
    ndim = 3
    ni = 1
    r = np.array([ [0.0,0.0,0.0], [1.0,0.0,0.0] ])
    t = np.array([ [1.0], [0.2] ])
    ilist = np.array([ [0,1] ])
    q = np.zeros((n,ndim),dtype=np.float,order='F')
    rho = np.zeros((n),dtype=float)
    grad_rho = np.zeros((n,ndim),dtype=float)
    mass = np.array([[1],[1]])
    sml = np.array([[2.0],[2.0]])
    drij = np.zeros([ni,ndim])
    rij = np.zeros(ni)

    for k in range(ni):
        i = ilist[k,0]
        j = ilist[k,1]
        drij[k,:] = r[j,:] - r[i,:]
        rij[k] = np.linalg.norm(drij[k,:])
    
    w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist+1,sml,2)
    fkernel.kernel.density_sum(rho,grad_rho,ilist+1,sml,mass,w,dwdx,2)
    splib.splib.calc_heat_flux_3d(q,ilist+1,rho,mass,t,-dwdx,1.0)
    print q


def density():
    """ Test the exposed density summation routine. Should give answers
        consistent with Python density summation script.
    """
    n = 2
    ni = 1 
    d = 3
    r = np.array( [ [0.0,0.0,0.0],
                    [1.0,0.0,0.0]   ]) #,
                    #[0.0,1.0,0.0],
                    #[0.0,0.0,1.0] ])

    nlist = np.zeros((ni,2))

    k = 0
    for i in range(n):
        for j in range(i+1,n):
            nlist[k,0] = i
            nlist[k,1] = j
            k += 1

    sml = np.zeros(n)
    rij = np.zeros(ni)
    sml[:] = 3.0
    w = np.zeros(ni)
    wp = np.zeros(ni)
    dwdx = np.zeros([ni,d])
    drij = np.zeros([ni,d])

    for k in range(ni):
        i = nlist[k,0]
        j = nlist[k,1]
        drij[k,:] = r[j,:] - r[i,:]
        rij[k] = np.linalg.norm(drij[k,:])
                        
    rho = np.zeros((n))
    grad_rho = np.zeros((n,d),order='F')

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,nlist+1,sml,2)

    for k in range(ni):
        i = nlist[k,0]
        j = nlist[k,1]
        wp[k] = spkernel.lucy_kernel(rij[k],drij[k,:],sml[i])[0]

    #wtf = spkernel.lucy_kernel(1.0,(1.0,1.0,0.0),2.0)[0]
    #bbq = fkernel_wrap(fkernel.kernel.lucy_kernel,1.0,h=2.0,ndim=3)
    #print wtf,bbq

    print 'fortran kernel',w
    print 'python kernel',wp

    mass = np.ones((n),dtype=float)
    fkernel.kernel.density_sum(rho,grad_rho,nlist+1,sml,mass,w,dwdx,2)
    rhopy = spdensity.sum_density(mass,sml,nlist,rij,drij)

    print 'fortran density',rho
    print 'python density',rhopy

    # What about the density gradient
    # Is it in the right direction!!???


    """
    xmin,xmax = -4,4
    ymin,ymax = -4,4
    res =0.2 
    cutoff = 2
    x,y = renderspam.get_grid_map(xmin,xmax,ymin,ymax,res)
    bounds = xmin,xmax,ymin,ymax
    Z = interpolate.splash_map_3d(x,y,0.0,r,mass,rho,rho,sml,bounds,cutoff) 
    Z = Z.transpose()
    plt.imshow(Z, origin='lower', extent=(xmin,xmax,ymin,ymax),
        interpolation='bilinear',vmin=0.0,vmax=2.0 )
    plt.plot(r[:,0],r[:,1],'ko')
    plt.axis([xmin,xmax,ymin,ymax])
    """

def density_gradient_test():
    # Initialise particles at an increasing x distance apart
    pass

    n = 125

    # Compute density and density gradient


    # Do a visual check


    # Then write an automated check




def tensors():
    """ Test subroutines in splib related to velocity gradient 
        and pressure tensor for a system of four particles.

        Four is the magic number for making sure the neighbour list calls
        don't have subtle errors.

        I've moved this code into a script inspect.py for a bit more of an
        indepth investigation aimed at ironing out all the
        inconsistencies between my math and the code.

        Starting data:
            -- a small array representing particle velocities
            -- a neighbour list integer array

        1. Compute velocity differences
            -- calc_dv(dv,nlist,v,n,ndim)

        2. Compute kernel weights and gradients

        1. Compute velocity gradient
            -- calc_grad_v(grad_v,ilist,dwdx,dv,m,rho,n,ndim,ni)
        
        2. Compute velocity divergence
            -- calc_div_v(grad_v,div_v,n,ndim)

        3. Compute the symmetric traceless velocity gradient
            -- calc_grad_v_os(grad_v_os,grad_v,div_v,n,ndim)

        4. Compute a symmetric traceless part of the pressure tensor
            -- calc_pi_os(pi_os,grad_v_os,eta,n,ndim)

        5. Compute the isotropic irreversible pressure
            -- calc_pi_one(pi_one,div_v,zeta,n)

        6. Now for the reversible pressure
            -- calc_vdw_hc_pressure(p,rho,u,n)
            -- calc_vdw_cohesive_pressure(pco,rho,u,n)

        7. And voila we have the full pressure tensor
        
    """

    # Data for a system of three particles
    v = np.array([     [1.0,0.0],
                       [0.0,3.14],
                       [20.0,0.0],
                       [0.0,7.0]
    ])

    n = 4
    ni = 5
    d = 2
    T = np.zeros([n])
    T[:] = 1.0
    nlist = np.array( [ [0,1], [0,2], [1,2], [0,3], [2,3] ],order='F' )
    sml = np.array( [ [2.0],[2.0],[2.0],[2.0] ] )
    w = np.array([[0.0],[0.0],[0.0],[0.0],[0.0]])
    dwdx = np.array([[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]])
    dv = np.zeros(dwdx.shape,order='F')
    rij = np.array([ [1.0], [1.0], [1.41], [3.0], [3.0] ])
    drij = np.array([ [1.0,0.0], [0.0,1.0], [1.0,1.0], [1.0,2.0], [2.0,1.0] ])
    rho = np.zeros((n))
    grad_rho = np.zeros((n,d))

    ilist = nlist + 1

    # Velocity diff
    splib.splib.calc_dv(dv,ilist,v)
    print 'dv'
    print dv

    for i in range(ni):
        dv[i,:] = v[nlist[i,0],:] - v[nlist[i,1],:]

    print dv

    # Kernels and kernel gradients
    w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,ilist,sml,1)
    print 'w'
    print w
    print 'dwdx'
    print dwdx

    # Density summation
    mass = np.ones((n))
    fkernel.kernel.density_sum(rho,grad_rho,ilist,sml,mass,w,dwdx,1)
    print 'rho'
    print rho
    print 'grad_rho'
    print grad_rho

    # Grad v
    grad_v = np.zeros((n,d,d),order='F')
    splib.splib.calc_grad_v(grad_v,ilist,dwdx,dv,mass,rho)
    print 'grad v'
    print grad_v

    # Div v
    div_v = np.zeros([n],order='F')
    splib.splib.calc_div_v(grad_v,div_v)
    print 'div v'
    print div_v

    #Symmetric traceless velocity gradient
    grad_v_os = np.zeros((n,d,d),order='F')
    splib.splib.calc_grad_v_os(grad_v_os,grad_v,div_v)

    print 'grad_v_os'
    print grad_v_os

    aos = symmetric_traceless(grad_v[0,:,:])
    print aos[:,:]

    #Compute a symmetric traceless part of the pressure tensor
    #    -- calc_pi_os(pi_os,grad_v_os,eta,n,ndim)
    pi_os = np.zeros((n,d,d),order='F')
    eta = np.zeros([n],order='F')
    eta[:] = 1
    splib.splib.calc_pi_os(pi_os,grad_v_os,eta)

    print 'pi_os'
    print pi_os

    #Compute the isotropic irreversible pressure
    #    -- calc_pi_one(pi_one,div_v,zeta,n)
    pi_one = np.zeros([n],order='F')
    zeta = np.zeros([n],order='F')
    zeta[:] = 1
    splib.splib.calc_pi_one(pi_one,div_v,zeta)

    print 'pi_one'
    print pi_one

    # Isotropic equilibrium (reversible) pressure
    # (python implementation)
    phc = np.zeros([n],dtype=float,order='F')
    pco = np.zeros([n],dtype=float,order='F')
    u = np.zeros([n],dtype=float,order='F')

    feos.eos.adash = 2.0
    feos.eos.bdash = 0.5
    feos.eos.kbdash = 1.0

    feos.eos.calc_vdw_energy(u,T,rho)
    feos.eos.calc_vdw_hc_pressure(phc,rho,u)
    feos.eos.calc_vdw_cohesive_pressure(pco,rho,u)

    #phc = vdw_hc(rho,T)
    #pco = vdw_co(rho,T)

    print 'phc'
    print phc

    print 'pco'
    print pco

    P = np.zeros([n,d,d],order='F')
    
    for i in range(n):
        for j in range(d):
            for k in range(d):
                P[i,j,k] = phc[i] + pco[i] +  pi_one[i] + pi_os[i,j,k]

    print 'pressure tensor'
    np.set_printoptions(precision=5,suppress=True)
    print P


def vel_grad():
    import splib
    gradv = np.array( [[1.0,2.0,3.0],
                        [4.0,5.0,6.0],
                        [7.0,8.0,9.0] ])
    gv = np.vstack((gradv,gradv))
    gv = np.vstack((gv,gradv))
    gv = np.vstack((gv,gradv))
    gv = gv.reshape((4,3,3))
    div = np.zeros([4])
    # Velocity divergence (really just vector trace...)
    splib.splib.calc_div_v(gv,div)
    return div


def fkernel_wrap(kfunc,r,h,ndim=1):
    """ Wraps fkernel subroutines such that they
        return a value like a function.
        
        The only argument is r, the distance.
        We are not concerned with the kernel derivative in this function.
        The returned value is w.

    """
    if ndim == 1:
        rr = np.atleast_1d(r)
        dr = np.atleast_1d(r)
        w = np.array([0.0])
        dwdx = np.array([0.0])
    if ndim==2:
        rr = np.atleast_1d(r)
        dr = np.array([r,0.0])
        w = np.array([0.0])
        dwdx = np.array([0.0,0.0])
    if ndim==3:
        rr = np.atleast_1d(r)
        dr = np.array([r,0.0,0.0])
        w = np.array([0.0])
        dwdx = np.array([0.0,0.0,0.0])
    kfunc(rr,dr,h,w,dwdx)
    return w

def fkernel_area(kfunc,range,h=2.0):
    """ Returns the normalised area under the kernel
        from zero to range[1].
    """
    from scipy import integrate
    kfwrap = fkernel_wrap(kfunc,r,h)
    area = integrate.quad(lambda r: kfwrap(r,h),0.0,h)
    pass

def kernel_areas():
    """ Compute normalisation and do some cross checking with the
        python implementation for the lucy kernel in one,two and three
        dimensions.
    """
    from scipy import integrate
    import numpy as np
    import spkernel
    h = 1.0

    # Lucy Kernel
    print 'Lucy Kernel'

    # One dimensional kernel areas
    kfc = fkernel.kernel.lucy_kernel
    area = integrate.quad(lambda r: fkernel_wrap(kfc,r,h=h),0.0,h)
    print 'Fortran 1D area:',2. * area[0]
    area2 = integrate.quad(lambda r: spkernel.lucy_kernel(r,r,h)[0],0.0,h)
    print 'Python 1D area:',2. * area2[0]

    # Two dimensional kernel areas
    area = integrate.quad(lambda r: r*fkernel_wrap(kfc,r,h=h,ndim=2),0.0,h)
    print 'Fortran 2D area:',2 * np.pi * (area[0]) 
    area2 = integrate.quad(lambda r: r*spkernel.lucy_kernel(r,[r,0.0],h)[0],0.0,h)
    print 'Python 2D area:', 2* np.pi * area2[0]

    # Three dimensional kernel volume
    area = integrate.quad(lambda r: r*r*fkernel_wrap(kfc,r,h=h,ndim=3),0.0,h)
    print 'Fortran 3D area:',4 * np.pi * (area[0]) 
    area2 = integrate.quad(lambda r: r*r*spkernel.lucy_kernel(r,[r,0.0,0.0],h)[0],0.0,h)
    print 'Python 3D area:', 4 * np.pi * area2[0]


    # Gaussian Kernel
    # Note that the Gaussian kernel needs to be greater than the smoothing
    # length

    print 'Gaussian Kernel'
    # One dimensional kernel areas
    kfc = fkernel.kernel.gauss_kernel
    area = integrate.quad(lambda r: fkernel_wrap(kfc,r,h=h),0.0,10.*h)
    print 'Fortran 1D area:',2. * area[0]
    area2 = integrate.quad(lambda r: spkernel.gauss_kernel(r,r,h)[0],0.0,5.*h)
    print 'Python 1D area:',2. * area2[0]

    # Two dimensional kernel areas
    area = integrate.quad(lambda r: r*fkernel_wrap(kfc,r,h=h,ndim=2),0.0,10.*h)
    print 'Fortran 2D area:',2 * np.pi * (area[0])
    area2 = integrate.quad(lambda r: r*spkernel.gauss_kernel(r,[r,0.0],h)[0]
        ,0.0,5.*h)
    print 'Python 2D area:',2 * np.pi * area2[0]

    # Three dimensional kernel areas
    area = integrate.quad(lambda r: r*r*fkernel_wrap(kfc,r,h=h,ndim=3),0.0,10.*h)
    print 'Fortran 3D area:',4 * np.pi * (area[0])
    area2 = integrate.quad(lambda r: r*r*spkernel.gauss_kernel(r,[r,0.0,0.0],h)[0]
        ,0.0,5.*h)
    print 'Python 3D area:',4 * np.pi * area2[0]


    print 'Debrun Kernel'
    # One dimensional kernel areas
    kfc = fkernel.kernel.debrun_kernel
    area = integrate.quad(lambda r: fkernel_wrap(kfc,r,h=h),0.0,10.*h)
    print 'Fortran 1D area:',2. * area[0]
    area2 = integrate.quad(lambda r: spkernel.debrun_kernel(r,r,h)[0],0.0,5.*h)
    print 'Python 1D area:',2. * area2[0]

    # Two dimensional kernel areas
    area = integrate.quad(lambda r: r*fkernel_wrap(kfc,r,h=h,ndim=2),0.0,10.*h)
    print 'Fortran 2D area:',2 * np.pi * (area[0])
    area2 = integrate.quad(lambda r: r*spkernel.debrun_kernel(r,[r,0.0],h)[0]
        ,0.0,5.*h)
    print 'Python 2D area:',2 * np.pi * area2[0]

    # Three dimensional kernel areas
    area = integrate.quad(lambda r: r*r*fkernel_wrap(kfc,r,h=h,ndim=3),0.0,10.*h)
    print 'Fortran 3D area:',4 * np.pi * (area[0])
    area2 = integrate.quad(lambda r: r*r*spkernel.debrun_kernel(r,[r,0.0,0.0],h)[0]
        ,0.0,5.*h)
    print 'Python 3D area:',4 * np.pi * area2[0]


def lucy_kernel_plot():
    """ Plot the lucy kernel, from several implementations. """
    from scipy import integrate
    import numpy as np
    import spkernel
    res = 100
    lims = 0,3
    h = 1.0

    w_lucy_fort = np.zeros(res,float)
    w_lucy_py = np.zeros(res,float)

    r = np.arange(lims[0],lims[1],float(lims[1]-lims[0])/res)
    for i in range(res):
        w_lucy_fort[i] = fkernel_wrap(fkernel.kernel.lucy_kernel,r[i],h=h)
        w_lucy_py[i] = spkernel.lucy_kernel(r[i],r[i],h)[0]

    plt.plot(r,w_lucy_fort,'r^-',label='Fortran Lucy')
    plt.plot(r,w_lucy_py,'bo-',label='Python Lucy')
    plt.legend()


def fkernel_plots():
    """ Plots of the three implemented kernels in two dimensions.
    """
    import fkernel
    res = 100
    lims = 0,3
    h = 1.0
    ndim = 2
    w_lucy = np.zeros(res,float)
    w_lucy_wrap = np.zeros(res,float)
    dwdx_lucy = np.zeros((res,ndim))
    w_gauss = np.zeros(res,float)
    dwdx_gauss = np.zeros((res,ndim))
    w_debrun = np.zeros(res,float)
    dwdx_debrun = np.zeros((res,ndim))

    r = np.arange(lims[0],lims[1],float(lims[1]-lims[0])/res)
    dx = np.zeros((r.size,ndim))
    dx[:,0] = r[:]

    # Compute the Lucy kernel
    for i in range(res):
        ww = np.zeros(1)
        fkernel.kernel.lucy_kernel(r[i],dx[i],h,ww,dwdx_lucy[i,:])
        w_lucy[i] = ww

    # Compute the Gaussian kernel
    for i in range(res):
        ww = np.zeros(1)
        fkernel.kernel.gauss_kernel(r[i],dx[i],h,ww,dwdx_gauss[i,:])
        w_gauss[i] = ww

    # Compute the Debrun kernel
    for i in range(res):
        ww = np.zeros(1)
        fkernel.kernel.debrun2d(r[i],dx[i],h,ww,dwdx_debrun[i,:])
        w_debrun[i] = ww
   
    plt.subplot(121)
    plt.plot(r,w_lucy,'r--',label='Lucy')
    plt.plot(r,w_gauss,'bs',label='Gaussian')
    plt.plot(r,w_debrun,'g^',label='Debrun')
    plt.legend()

    plt.subplot(122)
    plt.plot(r,dwdx_lucy[:,0],'r--',label='Lucy')
    plt.plot(r,dwdx_gauss[:,0],'bs',label='Gaussian')
    plt.plot(r,dwdx_debrun[:,0],'g^',label='Debrun')
    plt.legend()
    plt.show() 


def fkernel_norm():
    """ Compute area under the kernels to check normalisation.
    """
    from scipy import integrate
    import numpy as np

    r = np.arange(0,5,0.1)
    dr = np.array([ri])
    w = np.array(1.0)
    dwdx = np.array([1.0])

    # Need to write a function that returns the kernel value
    # it must be a wrapper of the fortran kernel function
    # lk = lambda r: fkernel.kernel.lucy_kernel(r,np.array([0.0]),2.0,w,dwdx)
    # integrate.quad(lk,0.0,3.0) 

    lk = lambda r,w,dwdx: kern(r,dr,2.0,w,dwdx)


def avisc_test():
    import art_viscosity
    avisc = 0.0
    dr = np.array([0.0,1.0])
    dv = np.array([0.0,-1.0])
    # These values are probably not consistent
    rho = 1.0
    sml = 1.0
    wij = 1.0
    rij = 1.0
    # Test 2 dimensional artificial viscosity
    print art_viscosity.art_viscosity.vnr_avisc(avisc,dr,dv,rho,sml,wij,rij)
    print art_viscosity.art_viscosity.mon_beta_avisc(avisc,dr,dv,rho,sml)
    # Test 3 dimensional artificial viscosity
    dr = np.array([0.0,1.0,0.0])
    dv = np.array([0.0,-1.0,0.0])
    print art_viscosity.art_viscosity.vnr_avisc(avisc,dr,dv,rho,sml,wij,rij)
    print art_viscosity.art_viscosity.mon_beta_avisc(avisc,dr,dv,rho,sml)

def collision_test():
    import collision

    # Test two dimensional collision
    va = np.array([0.0,1.0])
    vb = np.array([0.0,0.0])
    ma = 1.0
    mb = 1.0
    dx = np.array([0.0,1.0])
    rsq = 1.0
    print va,vb
    collision.collision.collide2d(va,vb,ma,mb,dx,rsq)
    print va,vb

    va = np.array([0.0,1.0,0.0])
    vb = np.array([0.0,0.0,0.0])
    ma = 1.0
    mb = 1.0
    dx = np.array([0.0,1.0,0.0])
    rsq = 1.0
    print va,vb
    from collision import collision
    collision.collide3d(va,vb,ma,mb,dx,rsq)
    print va,vb

def core_test():
    import core_potential
    #subroutine core_force(repulsion,dr,drsq,sigma,rcoef,ndim)
    # get the core force
    dr = np.array([0.0,1.0])
    rep = np.zeros([2])
    drsq = 1.0
    core_potential.core_potential.core_force(rep,dr,drsq,2.0,10.0)
    print rep 
    # ri = 0.0
    # rj = 1.0
    # dr = rj-ri = 1.0
    #a(i,1) = a(i,1) + repulsion(1)/m(i))
    #a(j,1) = a(j,1) - repulsion(1)/m(j))
    # so rcoef is negative for repulsion

if __name__ == '__main__':
    #tensors()
    # density()
    #fkernel_plots()
    #kernel_areas()
    #heat_flux()
    #core_test()
    force3d()
    #pairsep()
    #set_force3d_vdw()
    #set_force3d_vdw_params()






