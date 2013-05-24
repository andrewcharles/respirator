#!/usr/bin/python
""" 
    This module provides functions that calculate several versions of the sph
    smoothing kernel. 
    
    There should be two functions for each kernel. One just gets the kernel
    this is designed for doing interpolation. The other also gets the kernel
    gradient - this is designed for doing simulations. For simulation purposes we
    don't want to incur the function call overhead when doing this twice.
    
    For example: 
    
    gauss_kernel returns both the kernel and its derivative. 
    gauss_kernel_w returns only the kernel
    
"""


from math import *
import numpy as np
import fkernel 

ktable = {1:'gaussian',2:'lucy',3:'debrun'}

def kernel(r,dx,h,type):
    """ Calculates the kernel and kernel gradient for one pair of particles """
    """ 
    r -- distance
    x -- displacement vector
    h -- smoothing length
    type -- kernel type string: 'gaussian','lucy','debrun'
    """

    if type == 'gaussian':
        w,dwdx = gauss_kernel(r,dx,h)
    elif type == 'lucy':
        w,dwdx = lucy_kernel(r,dx,h)
    elif type == 'debrun':
        w,dwdx = debrun_kernel(r,dx,h)
    return w,dwdx

def gauss_kernel(r,dx,h):
    """ Ye old Gaussian SPH kernel. Returns w, dwdx """
    """ Works for as many dimensions as you care to inhabit
    """
    
    dx = np.atleast_1d(dx)
    dim = dx.size
    q = r/h
    factor = 1.e0 / ( (h**dim) * (np.pi**(dim/2.)) )
    w = factor * exp(-q*q)
    dwdx=[]
    for i in range(0,dim):
        dwdx[i:] = [w * ( (-2.e0) * (-1.e0) * dx[i]/(h*h))]
    return w,dwdx

def lucy_w(r,h):
    """ The Lucy Kernel. Returns W. Only 2d. """
    q = 5. / (pi * h**2)
    if r < 0:
        r = abs(r)
    if( r < h ):
        w = q * (1 + 3.*r/h)*((1.-r/h))**3
    else:
        w = 0
    return w

def lucy_w3d(r,h):
    """ The Lucy Kernel. Returns W. Only 3d
        This is a wrapper around my bonza fortran wrapped
        kernel. Why not start incorporating this stuff I
        says.

        I be worried that we allocate a new array here every time what a 
        flipping waste.

        r -- expect a just a distance
        kfc -- kernel function
        h -- smoothing length
    """
    kfc = fkernel.kernel.lucy_kernel
    w = np.array([0.0])
    dwdx = np.array([0.0,0.0,0.0])
    dr = np.array([r,0.0,0.0])
    kfc(r,dr,h,w,dwdx)
    return w


def lucy_kernel(r,dx,h):
    """ The Lucy Kernel. Returns W, dwdx"""
    """ Works for one, two and three dimensions
    """

    dx = np.atleast_1d(dx)
    dim = dx.size
    dwdx=[]

    if dim == 1:
        q = 5. / (4. * h)
    elif dim == 2:
        q = 5. / (pi * h**2)
    elif dim == 3:
        q = 105. / (pi * 16. * (h**3) )
    else:
        print 'ERROR'

    if r < 0:
        r = abs(r)
    if( r < h ):
        w = q * (1 + 3.*r/h) * ((1.-r/h))**3
        if(r == 0):
            for i in range(0,len(dx)):
                dwdx[i:] = [0]
        else:
            for i in range(0,dim):
                dwdx[i:] =  [ q * ( (-12./(h**4))*(r**3) +\
                             (24./(h**3))*(r**2) - (12.*r/(h**2)) )* dx[i]/r ]
    else:
        w = 0
        for i in range(0,len(dx)):
            dwdx[i:] = [0]
    return w,dwdx

def debrun_kernel(r,dx,h):
    """ The Debrun spiky kernel. Returns a tuple of w, dwdx """
    """ Works for two dimensions only
        Is not normalised correctly for one dimension?
    """

    #! gamma_g := 20/(Pi*h^6);
    #! w_spiky_2d := (r,h) -> piecewise((r<(h)), gamma_g * (h-r)^3, 0);
    dx = np.atleast_1d(dx)
    dim = dx.size
    dwdx=[] 
    if r < 0:
        r = abs(r)
    if( r < h ):
        if dim == 1:
            q = ( ( (-7.*(h**4)/4) + 2*(h**4) )**(-1) )/ 2  
            #(2.0/(3.0*(pi*(h**2)))) 
        elif dim == 2:
            q = (10.0/(pi*(h**5)))
        elif dim == 3:
            q = (15.0/(pi*(h**6)))
        w =  q * ((h - r)**3)
        if(r == 0):
            dwdx=0.0
        else:
            for i in range(0,dim):
                dwdx[i:] = [q * (-3.0*(h**2) + 6.0*h*r - 3.0*(r**2)) \
                              * (dx[i]/r) ]
    else:
        for i in range(0,dim):
            dwdx[i:] = [0]
        w = 0
    return w,dwdx

def test():
    """ Tests the module """
    from scipy import integrate

    print "Testing spkernel.py"

    print "Testing for one dimension"
    r,dx,h = 1.0,1.0,2.0
    type = "gaussian"
    
    # Gaussian kernel needs to be integrated further than smoothing length
    print "Gaussian"
    #print kernel(r,dx,h,"gaussian")
    #print gauss_kernel(r,dx,h)
    print 'area:',integrate.quad(lambda r: gauss_kernel(r,dx,h)[0],0.0,4.*h) 

    print "Debrun"
    #print kernel(r,dx,h,"debrun")
    #print debrun_kernel(r,dx,h)
    print 'area:',integrate.quad(lambda r: debrun_kernel(r,dx,h)[0],0.0,4.*h) 

    print "Lucy"
    #print kernel(r,dx,h,"lucy")
    #print lucy_kernel(r,dx,h)
    print 'area:',integrate.quad(lambda r: lucy_kernel(r,dx,h)[0],0.0,h) 


    # TWO DIMENSIONS
    print "Testing for two dimensions"
    r = 1.4
    dx = 1.0,1.0
    h = 3.0
    
    print "Gaussian"
    #print kernel(r,dx,h,"gaussian")
    #print gauss_kernel(r,dx,h)
    print 'area:',2*pi*integrate.quad(lambda r: r*gauss_kernel(r,dx,h)[0],0.0,10.*h)[0]

    print "Debrun"
    def dwrap(r,dr,h):
        rr = np.atleast_1d(r)
        dr = np.array([r,0.0])
        w = np.array([0.0])
        dwdx = np.array([0.0,0.0])
        fkernel.kernel.debrun_kernel(rr,dr,h,w,dwdx)
        return w
    print dwrap(r,dx,h)
    print debrun_kernel(r,dx,h)
    print 'area:',2*pi*integrate.quad(lambda r: r*debrun_kernel(r,dx,h)[0],0.0,10.*h)[0]
    print 'area:',2*pi*integrate.quad(lambda r: r*dwrap(r,dx,h),0.0,10.*h)[0]

    print "Lucy"
    #print kernel(r,dx,h,"lucy")
    #print lucy_kernel(r,dx,h)
    print 'area:',2*pi*integrate.quad(lambda r: r*lucy_kernel(r,dx,h)[0],0.0,h)[0]

    # Test for three dimensions
    # NB this isn't taking 3d integrals for the gaussian
    print "Testing for three dimensions"
    r = 1.73
    dx = 1.0,1.0,1.0
    h = 2.0
    
    print "Gaussian"
    #print kernel(r,dx,h,"gaussian")
    #print gauss_kernel(r,dx,h)
    print 'area:',4 * pi * \
        integrate.quad(lambda r: r*r*gauss_kernel(r,dx,h)[0],0.0,4.*h)[0] 

    print "Debrun"
    #print kernel(r,dx,h,"debrun")
    #print debrun_kernel(r,dx,h)
    print 'area:',4 * pi * \
        integrate.quad(lambda r: r*r*debrun_kernel(r,dx,h)[0],0.0,h)[0] 
    
    print "Lucy"
    #print kernel(r,dx,h,"lucy")
    #print lucy_kernel(r,dx,h)
    print 'area:',4 * pi * \
        integrate.quad(lambda r: r*r*lucy_kernel(r,dx,h)[0],0.0,2*h)[0]
        

if __name__ == "__main__":
    test()
