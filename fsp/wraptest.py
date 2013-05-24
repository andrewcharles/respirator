""" Tests the python wrappers around the fortran sph modules.

    The shell script wrapfort.py defines what python interfaces
    are generated.

    Some notes on usage:

    1. Derived types do not work.
    2. If a file has a module, that module is in the object path.

    ========================================================
    Python Library    |  Fortran Source     |   Test Present
    ========================================================
    fpairsep.so       |  pairsep.f90        |   Yes
    art_viscosity.so  |  art_viscosity.f90  |   Yes
    collision.so      |  collision.f90      |   Minimal
    fkernel.so        |  kernel.f90         |   Yes
    feos.so           |  eos.f90            |
    fsplib.so         |  splib.f90          |
    core_potential.so |  core_potential.f90 |   Minimal
    sphforce3d.so     |  sphforce3d.f90     |
    ========================================================

"""

import unittest
import sys
import numpy as np

class PairsepTest(unittest.TestCase):
    """ Test the neighbour list. """
    def test_one(self):
        import fpairsep
        from numpy import sqrt
        tolerance = 0.001
        x = np.array([ [ 1.0, 0.0, 0.0],
                       [ 0.0, 3.1, 0.0],
                       [20.0, 0.0, 0.0],
                       [ 0.0, 7.0, 0.0]])
        n = 4
        d = 3
        ni = (n*(n-1))/2
        nlist = np.zeros((ni,2))
        rij = np.zeros(ni)
        drij = np.zeros([ni,d])

        # Build a quick brute force neighbour list
        k = 0
        for i in range(n):
            for j in range(i+1,n):
                drx = x[j,0] - x[i,0]
                dry = x[j,1] - x[i,1]
                drz = x[j,2] - x[i,2]
                rsquared = drx**2 + dry**2 + drz**2
                drij[k,0] = drx
                drij[k,1] = dry
                drij[k,2] = drz
                rij[k] = sqrt(rsquared)
                nlist[k,0] = i
                nlist[k,1] = j
                k += 1
        ilist = nlist + 1
        frij = np.zeros(rij.shape)
        fpairsep.fpairsep.compute_pairsep(drij,frij,nlist+1,x)       
        self.failUnless(np.std(frij - rij) <= tolerance)


class AviscTest(unittest.TestCase):
    def test_one(self):
        """ Just tests that when we put numbers in we get numbers out."""
        import art_viscosity
        import numpy as np
        avisc = 0.0
        dr = np.array([0.0,1.0])
        dv = np.array([0.0,-1.0])
        # These values are probably not consistent
        rho = 1.0
        sml = 1.0
        wij = 1.0
        rij = 1.0
        # Test 2 dimensions
        art_viscosity.art_viscosity.vnr_avisc(avisc,dr,dv,rho,sml,wij,rij)
        art_viscosity.art_viscosity.mon_beta_avisc(avisc,dr,dv,rho,sml)
        # Test 3 dimensions
        dr = np.array([0.0,1.0,0.0])
        dv = np.array([0.0,-1.0,0.0])
        art_viscosity.art_viscosity.vnr_avisc(avisc,dr,dv,rho,sml,wij,rij)
        art_viscosity.art_viscosity.mon_beta_avisc(avisc,dr,dv,rho,sml)


class CollisionTest(unittest.TestCase):
    def test_collide2d(self):
        """ Simple collision test - two particles are one
            unit apart. The first has unit velocity.
            The expected result of the collision is for 
            all momentum to be transferred from the first
            to the second particle.
        """
        import collision
        import numpy as np
        tolerance = 0.001
        va = np.array([0.0,1.0])
        vb = np.array([0.0,0.0])
        ma = 1.0
        mb = 1.0
        dx = np.array([0.0,1.0])
        rsq = 1.0
        collision.collision.collide2d(va,vb,ma,mb,dx,rsq)
        self.failUnless(abs(np.linalg.norm(vb) - 1.0) <= tolerance )


class CoreTest(unittest.TestCase):
    def test_one(self):
        import core_potential


def fkernel_wrap(kfunc,r,h,ndim=1):
    """ Wraps fkernel subroutines such that they
        return a value like a function.
        
        The only argument is r, the distance.
        We are not concerned with the kernel derivative in this function.
        The returned value is w.

    """
    import numpy as np
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

class KernelTest(unittest.TestCase):
    def test_import(self):
        import fkernel
        from scipy.integrate import quad
        from numpy import pi

        h = 1.0
        tolerance = 0.001

        # Lucy Kernel area test
        kfc = fkernel.kernel.lucy_kernel

        # One dimension
        lucy_area_1d = 2 * quad(lambda r: fkernel_wrap(kfc,r,h=h),0.0,h)[0]
        self.failUnless((lucy_area_1d - 1.0) <= tolerance)

        # Two dimensions
        lucy_area_2d = 2 * pi \
            * quad(lambda r: r*fkernel_wrap(kfc,r,h=h,ndim=2),0.0,h)[0] 
        self.failUnless((lucy_area_2d - 1.0) <= tolerance)
        
        # Three dimensional kernel volume
        lucy_area_3d = 4 * pi * \
            quad(lambda r: r*r*fkernel_wrap(kfc,r,h=h,ndim=3),0.0,h)[0]
        self.failUnless((lucy_area_3d - 1.0) <= tolerance)

        # Gaussian Kernel area test
        # Note that the Gaussian kernel needs to be integrated further
        # out than the smoothing length
        kfc = fkernel.kernel.gauss_kernel
        # One dimension
        gauss_area_1d =  2 * \
            quad(lambda r: fkernel_wrap(kfc,r,h=h),0.0,10.*h)[0]
        self.failUnless((gauss_area_1d - 1.0) <= tolerance)
        # Two dimensions
        gauss_area_2d = 2 * pi \
            * quad(lambda r: r*fkernel_wrap(kfc,r,h=h,ndim=2),0.0,10.*h)[0] 
        self.failUnless((gauss_area_2d -1.0) <= tolerance)

        # Three dimensional kernel areas
        gauss_area_3d = 4 * pi * \
            quad(lambda r: r*r*fkernel_wrap(kfc,r,h=h,ndim=3),0.0,10.*h)[0]
        self.failUnless((gauss_area_3d -1.0) <= tolerance)

        # Debrun kernel area test
        # One dimensional kernel areas
        kfc = fkernel.kernel.debrun_kernel
        deb_area_1d = 2 * quad(lambda r: fkernel_wrap(kfc,r,h=h),0.0,10.*h)[0]
        self.failUnless((deb_area_1d -1.0) <= tolerance)

        # Two dimensional kernel areas
        deb_area_2d = 2 * pi * \
            quad(lambda r: r*fkernel_wrap(kfc,r,h=h,ndim=2),0.0,10.*h)[0]
        self.failUnless(abs(deb_area_2d -1.0) <= tolerance)
        print deb_area_2d

        # Three dimensional kernel areas
        deb_area_3d = 4 * pi * \
            quad(lambda r: r*r*fkernel_wrap(kfc,r,h=h,ndim=3),0.0,10.*h)[0]
        self.failUnless((deb_area_3d -1.0) <= tolerance)


class EosTest(unittest.TestCase):
    def test_import(self):
        import feos

class SplibTest(unittest.TestCase):
    def test_import(self):
        import fsplib

class Force3dTest(unittest.TestCase):
    def test_import(self):
        import sphforce3d

if __name__=='__main__':
	unittest.main()







