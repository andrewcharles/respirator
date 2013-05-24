""" Smooth particle rendering / field interpolation. 

    Andrew Charles 2013

"""

from math import *
import numpy
cimport numpy

import neighbours
import spkernel

import fkernel
cimport cython

ctypedef numpy.float_t DTYPE_t
ctypedef numpy.int_t DTYPE_int
DTYPE = numpy.float

cdef extern from "math.h":
    float sqrtf(float a)
    float powf(float a, float b)

def splash_map_2d(
	numpy.ndarray[DTYPE_t,ndim=1] x, 
	numpy.ndarray[DTYPE_t,ndim=1] y, 
	numpy.ndarray[DTYPE_t,ndim=2] r, 
	numpy.ndarray[DTYPE_t,ndim=1] m, 
	numpy.ndarray[DTYPE_t,ndim=1] rho, 
	numpy.ndarray[DTYPE_t,ndim=1] a, 
	numpy.ndarray[DTYPE_t,ndim=1] h,
	bounds,
	cutoff):
    """ Renders a 2D collection of particles to a 2D plane.

        x -- x positions of grid
        y -- y positions of grid
        r -- particle positions
        rho -- particle densities
        bounds -- box boundaries xmin,xmax,ymin,ymax
        cutoff -- smoothing cutoff

    """

    assert a.dtype == DTYPE

    cdef int nx,ny
    nx = x.size
    ny = y.size
    cdef numpy.ndarray[DTYPE_t,ndim=2,mode='c'] Z = numpy.zeros((nx,ny))
    cdef numpy.ndarray[DTYPE_t,ndim=1,mode='c'] dr, dwdx, w
    cdef float drx, dry, drz, xmax, ymax, xmin,ymin,rsq,res
    cdef int xpx,ypx
    # be careful trying to use unsigned ints for speed, as they will
    # not work if you make them negative!
    cdef int xpixmin,xpixmax,ypixmin,ypixmax,i

    if cutoff > nx: print "cutoff is too large this won't work"
    res = x[1]-x[0] #assume regular
    dr = numpy.zeros((2))
    dwdx = numpy.zeros((2))
    w = numpy.zeros((1))
    (xmin,xmax),(ymin,ymax) = bounds
    #loop over particle positions
    for i in xrange(m.size):
        # each particle contributes to pixels
        xpixmin = int(floor( (r[i,0] - cutoff - xmin)/res ) )
        ypixmin = int(floor( (r[i,1] - cutoff - ymin)/res ) )
        xpixmax = int(floor( (r[i,0] + cutoff - xmin)/res ) )
        ypixmax = int(floor( (r[i,1] + cutoff - ymin)/res ) )
        
        for ypx in xrange(ypixmin,ypixmax):
            for xpx in xrange(xpixmin,xpixmax):
               
                # where the pixel bounds cross the box bounds, need to 
                # wrap pixel coordinates
                if ypx >= ny: ypx -= ny
                if ypx < 0: ypx += ny
                if xpx >= nx: xpx -= nx
                if xpx < 0: xpx += nx
                 
                drx = r[i,0] - x[xpx]
                dry = r[i,1] - y[ypx]
                
				# need to apply a minimum image convention
                if (drx > xmax/2.):
                    drx = drx - xmax
                if (dry > ymax/2.):
                    dry = dry - ymax 	
                if (drx <= -xmax/2.):
                    drx = drx + xmax
                if (dry <= -ymax/2.):
                    dry = dry + ymax 	
                        
                rsq =  drx*drx + dry*dry
                dr[0] = drx
                dr[1] = dry
                if rsq < (cutoff*cutoff):
                    s = sqrtf(rsq)
                    fkernel.kernel.lucy_kernel(s,dr,h[i],w,dwdx,ndim=2)
                    # bug in cython means this doesn't work
                    #Z[xpx,ypx] +=  a[i]*(m[i]/rho[i])
                    Z[xpx,ypx] = Z[xpx,ypx] + a[i]*(m[i]/rho[i]) * w[0]
                                
    return Z


def splash_map_3d(
	numpy.ndarray[DTYPE_t,ndim=1] x, 
	numpy.ndarray[DTYPE_t,ndim=1] y, 
	numpy.float_t z, 
	numpy.ndarray[DTYPE_t,ndim=2] r, 
	numpy.ndarray[DTYPE_t,ndim=1] m, 
	numpy.ndarray[DTYPE_t,ndim=1] rho, 
	numpy.ndarray[DTYPE_t,ndim=1] a, 
	numpy.ndarray[DTYPE_t,ndim=1] h,
	bounds,
	cutoff): 
    """ I used Daniel Price's algorithm for this.
        res is the width of grid cells.
        This renders a 3D collection of particles to a 2D plane.

        x -- x positions of grid
        y -- y positions of grid
        z -- z coordinate of grid
        r -- particle positions
        rho -- particle densities
        bounds -- box boundaries xmin,xmax,ymin,ymax
        cutoff -- smoothing cutoff

    """

    assert a.dtype == DTYPE

    cdef int nx,ny
    nx = x.size
    ny = y.size
    cdef numpy.ndarray[DTYPE_t,ndim=2,mode='c'] Z = numpy.zeros((nx,ny))
    cdef numpy.ndarray[DTYPE_t,ndim=1,mode='c'] dr, dwdx, w
    cdef float drx, dry, drz, xmax, ymax, zmax, xmin,ymin,zmin,rsq,res
    cdef int xpx,ypx
    # be careful trying to use unsigned ints for speed, as they will
    # not work if you make them negative!
    cdef int xpixmin,xpixmax,ypixmin,ypixmax,i

    if cutoff > nx: print "cutoff is too large this won't work"
    res = x[1]-x[0] #assume regular
    dr = numpy.zeros((3))
    dwdx = numpy.zeros((3))
    w = numpy.zeros((1))
    (xmin,xmax),(ymin,ymax),(zmin,zmax) = bounds
    #loop over particle positions
    for i in range(m.size):
        # if the particle is too far in the z direction from the
        # plane then it is excluded
        #if (r[i,2] - z) > cutoff:
        #    continue
        #print i#,cutoff,xmin,res
        # each particle contributes to pixels
        xpixmin = int(floor( (r[i,0] - cutoff - xmin)/res ) )
        ypixmin = int(floor( (r[i,1] - cutoff - ymin)/res ) )
        xpixmax = int(floor( (r[i,0] + cutoff - xmin)/res ) )
        ypixmax = int(floor( (r[i,1] + cutoff - ymin)/res ) )
        for ypx in range(ypixmin,ypixmax):
            for xpx in range(xpixmin,xpixmax):
               
                # where the pixel bounds cross the box bounds, need to 
                # wrap pixel coordinates
                if ypx >= ny: ypx -= ny
                if ypx < 0: ypx += ny
                if xpx >= nx: xpx -= nx
                if xpx < 0: xpx += nx
               
                drx = r[i,0] - x[xpx]
                dry = r[i,1] - y[ypx]
                drz = r[i,2] - z
                
				# need to apply a minimum image convention
                if (drx > xmax/2.):
                    drx = drx - xmax
                if (dry > ymax/2.):
                    dry = dry - ymax 	
                if (drz > zmax/2.):
                    drz = drz - zmax
                if (drx < -xmax/2.):
                    drx = drx + xmax
                if (dry < -ymax/2.):
                    dry = dry + ymax 	
                if (drz < -zmax/2.):
                    drz = drz + zmax
                        
                rsq =  drx*drx + dry*dry + drz*drz
                dr[0] = drx
                dr[1] = dry
                dr[2] = drz
                if rsq < (cutoff*cutoff):
                    s = sqrtf(rsq)
                    fkernel.kernel.lucy_kernel(s,dr,h[i],w,dwdx,ndim=3)
                    # bug in cython means this doesn't work
                    #Z[xpx,ypx] +=  a[i]*(m[i]/rho[i])
                    Z[xpx,ypx] = Z[xpx,ypx] + a[i]*(m[i]/rho[i]) * w[0]
                                
    return Z


def splash_map_3dvol(
	numpy.ndarray[DTYPE_t,ndim=1] x, 
	numpy.ndarray[DTYPE_t,ndim=1] y, 
	numpy.ndarray[DTYPE_t,ndim=1] z, 
	numpy.ndarray[DTYPE_t,ndim=2] r, 
	numpy.ndarray[DTYPE_t,ndim=1] m, 
	numpy.ndarray[DTYPE_t,ndim=1] rho, 
	numpy.ndarray[DTYPE_t,ndim=1] a, 
	numpy.ndarray[DTYPE_t,ndim=1] h,
	bounds,
	cutoff): 
    """ I used Daniel Price's algorithm for this.
        res is the width of grid cells.
        This renders a 3D collection of particles to a 2D plane.

        x -- x positions of grid
        y -- y positions of grid
        z -- z coordinate of grid
        r -- particle positions
        rho -- particle densities
        bounds -- box boundaries xmin,xmax,ymin,ymax
        cutoff -- smoothing cutoff

    """

    assert a.dtype == DTYPE

    cdef int nx,ny,nz
    nx = x.size
    ny = y.size
    nz = z.size
    cdef numpy.ndarray[DTYPE_t,ndim=3,mode='c'] Z = numpy.zeros((nx,ny,nz))
    cdef numpy.ndarray[DTYPE_t,ndim=1,mode='c'] dr, dwdx, w
    cdef float drx, dry, drz, xmax, ymax, zmax, xmin,ymin,zmin,rsq,res
    cdef int xpx,ypx,zpx,i
    cdef int ypixmin,ypixmax

    if cutoff > nx: print "cutoff is too large this won't work"
    res = x[1]-x[0] #assume regular
    dr = numpy.zeros((3))
    dwdx = numpy.zeros((3))
    w = numpy.zeros((1))
    (xmin,xmax),(ymin,ymax),(zmin,zmax) = bounds
    #loop over particle positions
    for i in xrange(m.size):
        #print i,cutoff,xmin,res
        # each particle contributes to pixels
        xpixmin = int(floor( (r[i,0] - cutoff - xmin)/res ) )
        xpixmax = int(floor( (r[i,0] + cutoff - xmin)/res ) )
        
        ypixmin = int(floor( (r[i,1] - cutoff - ymin)/res ) )
        ypixmax = int(floor( (r[i,1] + cutoff - ymin)/res ) )
        
        zpixmin = int(floor( (r[i,2] - cutoff - zmin)/res ) )
        zpixmax = int(floor( (r[i,2] + cutoff - zmin)/res ) )
        for ypx in xrange(ypixmin,ypixmax):
            for xpx in xrange(xpixmin,xpixmax):
                for zpx in xrange(zpixmin,zpixmax):
                    
                    # where the pixel bounds cross the box bounds, need to 
                    # wrap pixel coordinates
                    if ypx >= ny: ypx -= ny
                    if ypx < 0: ypx += ny
                    if xpx >= nx: xpx -= nx
                    if xpx < 0: xpx += nx
                    if zpx >= nz: zpx -= nz
                    if zpx < 0: zpx += nz
               
                    drx = r[i,0] - x[xpx]
                    dry = r[i,1] - y[ypx]
                    drz = r[i,2] - z[zpx]
                    
                    # need to apply a minimum image convention
                    if (drx > xmax/2.):
                        drx = drx - xmax
                    if (dry > ymax/2.):
                        dry = dry - ymax 	
                    if (drz > zmax/2.):
                        drz = drz - zmax
                    if (drx < -xmax/2.):
                        drx = drx + xmax
                    if (dry < -ymax/2.):
                        dry = dry + ymax 	
                    if (drz < -zmax/2.):
                        drz = drz + zmax
                            
                    rsq =  drx*drx + dry*dry + drz*drz
                    dr[0] = drx
                    dr[1] = dry
                    dr[2] = drz
                    if rsq < (cutoff*cutoff):
                        s = sqrtf(rsq)
                        fkernel.kernel.lucy_kernel(s,dr,h[i],w,dwdx)
                        Z[xpx,ypx,zpx] = Z[xpx,ypx,zpx] + a[i]*(m[i]/rho[i]) * w[0]
                                
    return Z


