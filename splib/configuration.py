""" A library for generating configurations of particle positions.
"""

import numpy
from numpy import pi
from math import *

def randpt(xmin,xmax,ymin,ymax):
    """ returns a numpy array representing a random point """
    x = numpy.random.rand() * (xmax - xmin) + xmin
    y = numpy.random.rand() * (ymax - ymin) + ymin
    return numpy.array((x,y))

def grid(n,xside,yside,origin,spacing=1.0):
    """ Makes a regular spaced grid with the sides
        as specified.
        Origin is a tuple x,y coordinate.
    """
    n = int(n)
    x = numpy.zeros([n,2])
    j,k = origin
    for i in range(n):
        x[i,0]=j
        x[i,1]=k
        j=j+spacing
        if (i+1)%xside == 0:
            j=origin[0]
            k=k+spacing
    return x

def grid3d(n,side,origin,spacing=1.0):
    """ Makes a regular spaced grid with the sides
        as specified.
        side -- x,y,z tuple: number of points along each side
        origin -- a tuple x,y,z coordinate.
    """
    r = numpy.zeros([n,3])
    x,y,z = origin
    q = 0
    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = (x - (side[0]-1)*spacing/2.) + i * spacing
                r[q,1] = (y - (side[1]-1)*spacing/2.) + j * spacing
                r[q,2] = (z - (side[2]-1)*spacing/2.) + k * spacing
                q += 1
                if q >= n:
                    return r
    return r

def grid3d_uneven(n,side,origin,sx=1.0,sy=1.0,sz=1.0):
    """ Makes a regular spaced grid with the sides
        as specified.
        side -- x,y,z tuple: number of points along each side
        origin -- a tuple x,y,z coordinate.
    """
    r = numpy.zeros([n,3])
    x,y,z = origin
    q = 0
    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = (x - (side[0]-1)*spacing/2.) + i * sx 
                r[q,1] = (y - (side[1]-1)*spacing/2.) + j * sy
                r[q,2] = (z - (side[2]-1)*spacing/2.) + k * sz
                q += 1
                if q >= n:
                    return r
    return r


def fcc3d(n,side,origin,spacing=1.0):
    """ Create a 3 dimensional FCC grid. Assume periodicity,
        so that the interstices on one side are adjacent to
        the cell corners on the other side.

        s3 x 4 x (s1 x s2)

    """
    #if not n == side[2] * 4 * (side[0] + side[1]):
    #    print 'Error'
    #    print 'todo: throw exception'
    #    return None
    r = numpy.zeros([n,3])
    x,y,z = origin
    # Code for centering
    #x = x - (side[0]-1)*spacing/2.
    #y = y - (side[1]-1)*spacing/2.
    #z = z - (side[2]-1)*spacing/2.
    q = 0
    d = spacing
    for i in range(side[2]):
        # Make the first plane's corners
        for j in range(side[1]):
            for k in range(side[0]):
                r[q,0] = x + (k * d)
                r[q,1] = y + (j * d)
                r[q,2] = z + (i * d)
                q += 1
        # Make the first plane's interstices
        for j in range(side[1]):
            for k in range(side[0]):
                r[q,0] = x + ((k + 0.5) * d)
                r[q,1] = y + ((j + 0.5) * d)
                r[q,2] = z + (i * d)
                q += 1
        # Make a plane at 1/2 d k
        # Make the first plane's corners
        for j in range(side[1]):
            for k in range(side[0]):
                r[q,0] = x + ((k + 0.5) * d)
                r[q,1] = y + ((j) * d)
                r[q,2] = z + ((i + 0.5) * d)
                q += 1
        # Make the first plane's interstices
        for j in range(side[1]):
            for k in range(side[0]):
                r[q,0] = x + ((k) * d)
                r[q,1] = y + ((j + 0.5) * d)
                r[q,2] = z + ((i + 0.5) * d)
                q += 1
    return r

def twophase3d(n,nh,side,sideh,origin,spacingh=0.7,spacing=2.0):
    """ Distributes particles in a high density and low
        density rectangular phase.
        n -- number
        nh -- number of high density
        side -- total side lengths
        sideh -- side lengths of the high density phase

        e.g.

        box is 10x5x5
        First half is high density
        5x5x5 with spacing 0.25
        That's 5x4 cubed particles = 8000
        20x20x20 particles
        Second half is low density
        5x5x5 with spacing 1
        That's 5 cubed particles = 125

    """
    # First make the high density phase
    rhigh = grid3d(nh,sideh,origin,spacing=spacingh)
    # low phase is translated in the x direction
    origin2 = (origin[0] + sideh[0]*spacingh,origin[1],origin[2])
    rlow = grid3d(n-nh,side,origin2,spacing=spacing)
    r = numpy.concatenate((rhigh,rlow),axis=0)
    return r

def hotspotgrid3d(n,side,origin,spacing=1.0,temp=(0.2,1.9)):
    """ Makes a regular spaced grid with the sides
        as specified. The central third of points are hot.
        side -- x,y,z tuple: number of points along each side
        origin -- a tuple x,y,z coordinate.
    """
    r = numpy.zeros([n,3])
    t = numpy.zeros([n])
    x,y,z = origin
    q = 0
    # get the thirds
    xt = side[0]/3,2*side[0]/3
    yt = side[1]/3,2*side[1]/3
    zt = side[2]/3,2*side[2]/3
    if side[2] == 1: zt = 0,0
    print xt,yt,zt
    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = (x - (side[0]-1)*spacing/2.) + i * spacing
                r[q,1] = (y - (side[1]-1)*spacing/2.) + j * spacing
                r[q,2] = (z - (side[2]-1)*spacing/2.) + k * spacing
                if i>=xt[0] and i<=xt[1] and j>=yt[0] and j<=yt[1] \
                    and k>=zt[0] and k<=zt[1]:
                    t[q] = temp[1]
                else:
                    t[q] = temp[0]
                q += 1
                if q >= n:
                    return r,t
    return r,t
    
def random(n,xmin,xmax,ymin,ymax):
    """ 
    """
    r = numpy.zeros((n,2))
    width = xmax - xmin
    height = ymax - ymin
    for i in range(n):
       r[i,:] = numpy.random.random()*width + xmin,\
                numpy.random.random()*height + ymin 
    return r

def random3d(n,xmin,xmax,ymin,ymax,zmin,zmax):
    """ 
    """
    r = numpy.zeros((n,3))
    width = xmax - xmin
    height = ymax - ymin
    depth = zmax - zmin
    for i in range(n):
       r[i,:] = numpy.random.random()*width + xmin,\
                numpy.random.random()*height + ymin,\
                numpy.random.random()*depth + zmin
    return r

from integrator import rk4
import particles, neighbour_list, forces
import numpy as np
def disorder3d(n,xmin,xmax,ymin,ymax,zmin,zmax):
    """ Make a random configuration 
        with spatial dimension xmin,xmax,ymin,ymax,zmin,zmax
        Relax the configuration back a little by applying a
        repulsive force.
    """
    r = random3d(n,xmin,xmax,ymin,ymax,zmin,zmax)
    p = particles.ParticleSystem(n,3,maxn=n,integrator='rk4',vmax=0.0)
    p.r = r
    nl = neighbour_list.FastVerletList(p,cutoff=20.0)
    nl.build()
    nl.separations()
    p.nlists.append(nl)
    #p.forces.append(forces.CoreForce3d(p,nl,rc=1.0,sig=4.0))
    p.forces.append(forces.Gravity3d(p,nl,g=-1.0))
    for i in range(100):
        p.update(0.01)
    return p.r

def test_disorder2():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    #r = fcc3d(8*4,(2,2,2),(0.5,0.5,0.5),spacing=0.1)
    #r = twophase3d(20*4*4,14*4*4,(20,4,4),(14,4,4),
    #    (0.5,0.5,0.5),spacing=0.1)
    #r = disorder3d(10,0,10.,0,10.,0,10.)
    r = random3d(10,0,10.,0,10.,0,10.)
    p = particles.ParticleSystem(n,3,maxn=n,integrator='rk4',vmax=0.0)
    p.r = r
    nl = neighbour_list.FastVerletList(p,cutoff=20.0)
    nl.build()
    nl.separations()
    p.nlists.append(nl)
    #p.forces.append(forces.CoreForce3d(p,nl,rc=1.0,sig=4.0))
    p.forces.append(forces.Gravity3d(p,nl,g=-1.0))
    for i in range(100):
        p.update(0.01)

    ax = axes3d.Axes3D(figure())
    ax.plot3D(r[:,0],r[:,1],r[:,2],'o')
    plt.figure()
    from neighbours import build_list, separations, distances
    iap = build_list(r)
    rij = separations(r,iap)
    r = distances(rij)
    n,bins,patches = plt.hist(r,bins=25,normed=True)
    plt.xlim(0.0,2.0)
    print n[0],bins[0]
    ion()
    #rhigh = grid3d(nh,sideh,origin,spacing=spacingh)
    show()
    return r

def test_fcc():
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    #s3 x 4 x (s1 x s2)
    #r = fcc3d(N,(5,5,5),orign,spacing=0.7)
    #r = fcc3d(5*4*25,(5,5,5),(0.5,0.5,0.5),spacing=0.7)
    r = fcc3d(4*(9*9*9),(9,9,9),(0.5,0.5,0.5),spacing=1.4)
    ax = axes3d.Axes3D(figure())
    ax.plot3D(r[:,0],r[:,1],r[:,2],'o')
    plt.figure()
    from neighbours import build_list, separations, distances
    iap = build_list(r)
    rij = separations(r,iap)
    r = distances(rij)
    n,bins,patches = plt.hist(r,bins=25,normed=True)
    plt.xlim(0.0,2.5)
    print n[0],bins[0]
    ion()
    show()
    return r

def test_disorder():
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    #r = fcc3d(8*4,(2,2,2),(0.5,0.5,0.5),spacing=0.1)
    #r = twophase3d(20*4*4,14*4*4,(20,4,4),(14,4,4),
    #    (0.5,0.5,0.5),spacing=0.1)
    r = disorder3d(10,0,10.,0,10.,0,10.)
    ax = axes3d.Axes3D(figure())
    ax.plot3D(r[:,0],r[:,1],r[:,2],'o')
    plt.figure()
    from neighbours import build_list, separations, distances
    iap = build_list(r)
    rij = separations(r,iap)
    r = distances(rij)
    n,bins,patches = plt.hist(r,bins=25,normed=True)
    plt.xlim(0.0,2.0)
    print n[0],bins[0]
    ion()
    #rhigh = grid3d(nh,sideh,origin,spacing=spacingh)
    show()
    return r

def test_twophase():
    n=8125
    nh=8000
    sidel=(5,5,5)
    sideh=(20,20,20)
    origin=(0.,0.,0.)
    spacing=1.0
    spacingh=0.25
    # First make the high density phase
    rhigh = grid3d(nh,sideh,origin,spacing=spacingh)
    # low phase is translated in the x direction
    origin2 = (origin[0] + sideh[0]*spacingh,origin[1],origin[2])
    rlow = grid3d(n-nh,sidel,origin2,spacing=spacing)
    print rhigh.shape
    print rlow.shape
    r = numpy.concatenate((rhigh,rlow),axis=0)
    #return r
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    ax = axes3d.Axes3D(figure())
    ax.plot3D(r[:,0],r[:,1],r[:,2],'o')
    ax.plot3D(rhigh[:,0],rhigh[:,1],rhigh[:,2],'go')
    ax.plot3D(rlow[:,0],rlow[:,1],rlow[:,2],'ro')
    ion()
    show()

def cosline(n,xmin=0.0,xmax=1.0,scale=0.2):
    """ n -- number of particles

    """
    xp = np.linspace(0,2*pi,n)
    cxp = (1.0 - scale) + scale * np.cos(xp)
    xt = (xp - pi) * cxp + pi
    xt = xt * (xmax / (2*pi)) + xmin
    return xt

def cos3dx(side,scale=0.2,spacing=1.0,origin=(0,0,0)):
    """ Telescopes in the x dimension. """
    n = side[0]*side[1]*side[2]
    xr = side[0] * spacing
    r = numpy.zeros([n,3])
    x0,y,z = origin
    x = np.zeros([side[0]])
    dx = np.zeros([side[0]])
    tf = 0.5
    q = 0

    # Indices of each plane
    xcoords = cosline(side[0],xmin=origin[0]-(side[0]*spacing)/2.
                             ,xmax=origin[0]+(side[0]*spacing)/2.
                             ,scale=scale)

    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = xcoords[i]
                r[q,1] = (y - (side[1]-1)*spacing/2.) + j * spacing
                r[q,2] = (z - (side[2]-1)*spacing/2.) + k * spacing
                q += 1
                if q >= n:
                    return r
    return r

def cos3d(side,scale=0.2,spacing=1.0,origin=(0,0,0)):
    """ Telescopes in all three dimensions. """
    n = side[0]*side[1]*side[2]
    xr = side[0] * spacing
    r = numpy.zeros([n,3])
    x0,y,z = origin
    x = np.zeros([side[0]])
    dx = np.zeros([side[0]])
    tf = 0.5
    q = 0

    # Indices of each plane
    xcoords = cosline(side[0],xmin=origin[0]-(side[0]*spacing)/2.
                             ,xmax=origin[0]+(side[0]*spacing)/2.
                             ,scale=scale)

    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = xcoords[i]
                r[q,1] = xcoords[j]#(y - (side[1]-1)*spacing/2.) + j * spacing
                r[q,2] = xcoords[k]#(z - (side[2]-1)*spacing/2.) + k * spacing
                q += 1
                if q >= n:
                    return r
    return r

def sinline(n,k=1,xmin=0.0,xmax=1.0,scale=0.2):
    """ n -- number of particles
        k -- wavenumber

            1.  Create an equally spaced domain xp from 0 to k*Pi
            2.  Perturbation cxp is given by sin(xp), scaled and translated
                so that the perturbation wave sits just under y=1
            3.  Positions xt are given by xp from 0 to k*pi
            4.  Positions xt are scaled by k*pi

    """
    xp = np.linspace(0,k*pi,n)
    cxp = 1.0 + scale*np.sin(xp)
    xt = (xp) * cxp
    xt = xt * (xmax / (k*pi)) + xmin
    return xt

def sin2d(side,scale=0.2,spacing=1.0,k=1,origin=(0,0)):
    """ Telescopes in all three dimensions. """
    n = side[0]*side[1]*side[2]
    xr = side[0] * spacing
    r = numpy.zeros([n,3])
    x0,y,z = origin
    x = np.zeros([side[0]])
    dx = np.zeros([side[0]])
    tf = 0.5
    q = 0

    # Indices of each plane
    xcoords = sinline(side[0],xmin=origin[0]-(side[0]*spacing)/2.
                             ,xmax=origin[0]+(side[0]*spacing)/2.
                             ,scale=scale,k=k)

    for i in range(side[0]):
        for j in range(side[1]):
            r[q,0] = xcoords[i]
            r[q,1] = xcoords[j]
            q += 1
            if q >= n:
                return r
    return r

def sin3d(side,scale=0.2,spacing=1.0,k=(1,1,1),origin=(0,0,0)):
    """ Telescopes in all three dimensions. """
    n = side[0]*side[1]*side[2]
    xr = side[0] * spacing
    r = numpy.zeros([n,3])
    x0,y,z = origin
    x = np.zeros([side[0]])
    dx = np.zeros([side[0]])
    tf = 0.5
    q = 0

    # Indices of each plane
    xcoords = sinline(side[0],xmin=origin[0]
                             ,xmax=origin[0]+(side[0]*spacing)
                             ,scale=scale,k=k[0])

    ycoords = sinline(side[1],xmin=origin[1]
                             ,xmax=origin[1]+(side[1]*spacing)
                             ,scale=scale,k=k[1])

    zcoords = sinline(side[2],xmin=origin[2]
                             ,xmax=origin[2]+(side[2]*spacing)
                             ,scale=scale,k=k[2])

    for i in range(side[0]):
        for j in range(side[1]):
            for k in range(side[2]):
                r[q,0] = xcoords[i]
                r[q,1] = ycoords[j]
                r[q,2] = zcoords[k]
                q += 1
                if q >= n:
                    return r
    return r


def test_sin3d():
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    r = sin3d((20,20,20),k=(1,1,1),spacing=1.0,scale=0.1)
    ax = axes3d.Axes3D(figure())
    ax.plot3D(r[:,0],r[:,1],r[:,2],'o')
    plt.xlim(0.0,20.0)
    plt.ylim(0.0,20.0)
    #plt.zlim(0.0,10.0)
    ion()
    show()
    return r 


def test_cos3d():
    from mpl_toolkits.mplot3d import axes3d
    from pylab import figure, show, ion
    r = cos3d((10,10,10),k=(1,1,1),spacing=1.0,scale=0.5)
    ax = axes3d.Axes3D(figure())
    ax.plot3D(r[:,0],r[:,1],r[:,2],'o')
    from neighbours import build_list, separations, distances
    #iap = build_list(r)
    #rij = separations(r,iap)
    #rs = distances(rij)
    #n,bins,patches = plt.hist(rs,bins=25,normed=True)
    plt.xlim(0.0,2.5)
    #print n[0],bins[0]
    ion()
    show()
    return r 

from math import exp
def test_gauss1d():
    """ Distribute with gaussian adjusted spacing
        over a unit cell.
    """
    plt.clf()
    side = 10
    dx0 = 0.1
    tf = 0.9
    xr = side * dx0 
    dx = np.zeros(side)
    x = np.zeros(side)
    x[0] = 0.0
    dx[0] = dx0

    for i in range(1,side):
        q = exp( - (i*dx0 - (0.0))**2 / (xr/2.)**2)
        dx[i] = dx0 - tf * q 
        x[i] = x[i-1] + dx[i]

    print x
    print dx
    print xr
    x = x - (side*dx0)/2.
    zer = np.zeros(side)
#    plt.plot(x,zer,'ko')
    #plt.ylim(-0.1,0.1)
#    plt.xlim(-(side-1)*(dx0/2.),(side)*(dx0/2.),dx0)
    plt.plot(x,dx,'ko')


def smooth_particle_test():
    pass


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        if sys.argv[1] == 'twophase':
            r = test_twophase()
        if sys.argv[1] == 'sph':
            r = smooth_particle_test()
        if sys.argv[1] == 'fcc':
            r = test_fcc()
        if sys.argv[1] == 'gauss':
            r = test_gauss()
    else:
        #r = test_disorder2()
        #r = test_cos3d()
        r = test_sin3d()


