#! /usr/local/bin/python
""" Functions for rendering a spam field in two dimensions.

    Andrew Charles 2013.

"""

import math
import numpy
import spkernel
import splib.interpolate as interpolate
import neighbours
import pylab
import configuration

NPTS = 100
# resolution required
nx = 0
ny = 0
ngridpts = 0

def map_to_grid():
    """ Maps particle coordinates to grid coordinates """
    print "Not implemented"

def get_grid_map(xmin,xmax,ymin,ymax,res):
    """ Returns the coordinates of the centers of the
        grid points.
        res -- grid spacing (Smaller res = finer grid).
    """
    gridx = numpy.arange(0+res/2.0,xmax,res)
    gridy = numpy.arange(0+res/2.0,ymax,res)
    return gridx,gridy 

def get_grid_map3d((xmin,xmax),(ymin,ymax),(zmin,zmax),res):
    """ Returns a the coordinates of grid point centers for a 
        set of axes with uniform point spacing
        res -- grid spacing (Smaller res = finer grid).
    """
    gridx = numpy.arange(0+res/2.0,xmax,res)
    gridy = numpy.arange(0+res/2.0,ymax,res)
    gridz = numpy.arange(0+res/2.0,zmax,res)
    return gridx,gridy,gridz

def func2(x,y):
    return numpy.exp(x+y)

def renderspam2d(x,y,r,m,rho,a,h,xmin,xmax,ymin,ymax,cutoff):
    """ x: x coordinates of the grid
        y: y coordinates of the grid
        r: coordinates of points [x,y]
        rho: sph summation densities of the points
        m: mass of each point
        a: property we are rendering (of each point) 
        xmin,xmax,ymin,ymax
    """

    bounds = (xmin,xmax),(ymin,ymax)
   
    #  Z = interpolate.map_to_grid(x,y,r,m,rho,a,h,bounds,cutoff) 
    Z = interpolate.splash_map_2d(x,y,r,m,rho,a,h,bounds,cutoff) 
    Z = Z.transpose()
    return pylab.imshow( Z , origin='lower', extent=(xmin,xmax,ymin,ymax),interpolation='bilinear',vmin=0.0,vmax=2.0)# , aspect='auto' )


def main():
    """ Plots a regular 2d spam field
    """

    # initial positions
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 10
    n = 10
    cutoff = 3
    r = configuration.random(n,xmin,xmax,ymin,ymax)

    h = numpy.zeros((r.shape[0]))
    m = numpy.zeros((r.shape[0]))
    rho = numpy.zeros((r.shape[0]))
    h[:] = 2
    rho[:] = 1.0
    m[:] = 1.0
    #ry = numpy.array([1,1,1,2,2,2,3,3,3])
    gridx,gridy = get_grid_map(xmin,xmax,ymin,ymax,0.10)
    pylab.ion()
    pylab.axis([xmin,xmax,ymin,ymax])
    renderspam2d(gridx,gridy,r,m,rho,rho,h,xmin,xmax,ymin,ymax,cutoff)
    pylab.plot(r[:,0],r[:,1],'wo')
    pylab.axis([xmin,xmax,ymin,ymax])
    pylab.show()

if __name__ == "__main__":
    main()

