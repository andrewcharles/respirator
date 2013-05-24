"""
    Particle system model. Contains data structures and basic functions
    for a particle system and associated neighbour list.

    Designed so that the neighbour list and force
    modules are accessed through here.

    Classes
    -------
    ParticleSystem
    SmoothParticleSystem (inherits from ParticleSystem)

    Copyright Andrew Charles 2008
    All rights reserved.
    This module is new BSD licensed.

"""

import random
import numpy as np
import math
#import forces
import neighbour_list
#import properties
import eos
import netCDF4

# I don't want to have to select this here!
#import c_properties as cproperties
import f_properties as fproperties
import properties
#from properties import vdw_energy,spam_properties,hamiltonian
import scipy
import integrator
import configuration
import box
import sys
#from Numeric import *
#import pdb
from integrator import rk4, euler, imp_euler
from time import time
from spam_nc import read_step
import os
#dt = 0.05

# variables for the integrator - put these somewhere cleaver
verbose = False

XMAX = 500 # 64
YMAX = 500 # 48
ZMAX = 500
VMAX = 0.1
CUTOFF_RADIUS = 6 
VACUUM_VISCOSITY = 0.1
RAMAL = 0.5      # Amalgamation radius
VSPLIT = 10.0 
rhosplit = 0.0
ADKE = False     # Sigalotti style adaptive density kernel estimation
AMALGAMATE = False
SPLIT = False
ADVECTIVE = False
SPROPS = False

class ParticleInputReader:
    """ A module that reads data for the particle system. """
    def __init__(self):
        pass
    def read_netcdf(self,filename):
        """ Read position, velocity, internal energy, 
            mass from a netcdf file. """
        pass

class ParticleSystem:
    """ A group of similar particles with basic mechanical properties.
        This basic implementation is two dimensional.
    """
    def __init__(self,n,
            maxn=125,
            controllers=[],
            xmax=XMAX,
            ymax=YMAX,
            vmax=VMAX,
            simbox=None,
            mass=1.0,
            rinit=None,
            source=None,
            side=(5,5),
            integrator='rk4',
            twophase=None,
            spacing=0.1):
        """
        DIMENSIONS
        n -- initial number of particles.
        maxn -- maximum number of particles.
        steps -- number of steps taken
        nlists -- neighbour lists associated with this particle system.
        colour -- a 3 tuple giving the particles' RBG color.
        box -- the simulation box. Should replace this with constraint forces.
        nlists -- neighbour lists associated with this system.
        forces -- internal forces associated with this system.
        configuration -- initial positions.
        side -- side lengths (only for grid positions)
        rinit -- string code for initialisation strategy
               grid
               fcc
        """
        d=2
        self.n = n
        self.dim = d
        self.maxn = maxn
        self.dt = 0.0
        self.steps = 0

        # Basic mechanical properties
        if simbox is None:
            print 'No box, setting mirror bounds'
            self.box = box.MirrorBox2d(p=self,xmax=xmax,ymax=ymax)
        else:
            self.box = simbox
            self.box.p = self
       
        # Select integrator
        integrator_mapping = {'euler':euler,
                              'ieuler':imp_euler,
                              'rk4':rk4}

        self.step = integrator_mapping[integrator]

        # Start with a random configuration
        self.r = self.box.xmax * np.random.random([self.maxn,self.dim])
        self.m = np.zeros(self.maxn,dtype=float)
        vrand = vmax * np.random.random([self.maxn,self.dim])
        vnorm = np.zeros(vrand.shape)
        vnorm[:,0] = vrand[:,0] - vrand[:,0].mean()
        vnorm[:,1] = vrand[:,1] - vrand[:,1].mean()
        self.v = vnorm
        self.rdot = np.zeros(self.r.shape)
        self.vdot = np.zeros(self.v.shape)
        self.mdot = np.zeros(self.m.shape)
    
        if rinit == 'random':
            self.r[0:n,:] = configuration.random(n,0,xmax,0,ymax)
            self.m[:] = mass 

        if rinit == 'grid':
            self.r[0:n,:] = configuration.grid(
                n,side[0],side[1],(xmax/2.,ymax/2.),spacing=spacing)
            self.m[:] = mass 

        elif rinit == 'fcc':
            print 'UNSUPPORTED'

        elif rinit == 'twophase':
            print 'UNSUPPORTED'

        elif rinit == 'load':
            print 'UNSUPPORTED'

        else:
            self.r[0:n,:] = configuration.random(n,0,xmax,0,ymax)
            self.m[:] = mass

        # Initialise values
        self.colour = 1.0,0.0,0.0 

        # State vectors to pass to numerical integrators
        # 
        n_variables = 5
        self.x = np.zeros([n_variables,self.maxn])
        self.xdot = np.zeros([n_variables,self.maxn])

        self.nlists = []
        self.forces = []

        self.controllers = controllers
        for controller in self.controllers:
            controller.bind_particles(self)

        """ Variables for measuring performance. """
        self.timing = {}
        self.timing['force time'] = -1
        self.timing['pairsep time'] = -1
        self.timing['update time'] = -1
        self.timing['integrate time'] = -1

    def create_particle(self,r,v=(0.0,0.0,0.0)):
        """Adds a new particle to the system.
        """
        self.r[self.n] = r
        self.m[self.n] = self.m[self.n-1] 
        self.v[self.n] = v
        self.n = self.n+1
        self.rebuild_lists()

    def rebuild_lists(self):
        """ rebuilds all nlists """
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build()

    def update(self,dt):
        """ Update the particle system, using the
            neighbour list supplied.
        """
        self.rebuild_lists()
        t = time()
        self.step(self.gather_state,self.derivatives, \
                       self.gather_derivatives,self.scatter_state,dt)
        self.timing['integrate time'] = time() - t
        self.box.apply(self)
        self.steps += 1

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.v[0:self.n,0]
        self.x[4,0:self.n] = self.v[0:self.n,1]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.v[0:self.n:,0] = x[3,0:self.n]
        self.v[0:self.n:,1] = x[4,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = self.mdot[0:self.n] 
        self.xdot[1,0:self.n] = self.rdot[0:self.n,0]
        self.xdot[2,0:self.n] = self.rdot[0:self.n,1]
        self.xdot[3,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,1]
        return self.xdot

    def derivatives(self):
        """ Compute the rate of change of each variable 
            for every particle. The force.apply() call
            accumulates the forces.
            This is the expensive part of the program.
            
        """
        self.rdot[:,:] = self.v[:,:]
        self.vdot[:,:] = 0.0
    
        t = time()
        for nl in self.nlists: 
            nl.separations()
        self.timing['nlist sep time'] = time() - t
        
        t = time()
        for force in self.forces:
            force.apply()
        self.timing['force2 time'] = time() - t

        t = time()
        for controller in self.controllers:
            controller.apply()
        self.timing['controller time'] = time() - t


class ParticleSystem3D(ParticleSystem):
    """ A group of similar particles with basic mechanical properties.
    """
    def __init__(self,n,
            d=3,
            maxn=125,
            controllers=[],
            xmax=XMAX,
            ymax=YMAX,
            zmax=ZMAX,
            vmax=VMAX,
            simbox=None,
            mass=1.0,
            rinit=None,
            source=None,
            side=(5,5,5),
            integrator='rk4',
            twophase=None,
            spacing=0.1):
        """
        DIMENSIONS
        n -- initial number of particles.
        maxn -- maximum number of particles.
        dim -- number of spatial dimensions (1-3).
        steps -- number of steps taken
        nlists -- neighbour lists associated with this particle system.
        colour -- a 3 tuple giving the particles' RBG color.
        box -- the simulation box. Should replace this with constraint forces.
        nlists -- neighbour lists associated with this system.
        forces -- internal forces associated with this system.
        configuration -- initial positions.
        side -- side lengths (only for grid positions)
        rinit -- string code for initialisation strategy
               grid
               fcc
        """
        self.n = n
        self.dim = d
        self.maxn = maxn
        self.dt = 0.0
        self.steps = 0

        # Basic mechanical properties
        if simbox is None:
            print 'No box, setting mirror bounds'
            self.box = box.MirrorBox(p=self,xmax=xmax,ymax=ymax,zmax=zmax)
        else:
            self.box = simbox
            self.box.p = self
       
        # Select integrator
        integrator_mapping = {'euler':euler,
                              'ieuler':imp_euler,
                              'rk4':rk4}

        self.step = integrator_mapping[integrator]

        # Start with a random configuration
        self.r = self.box.xmax * np.random.random([self.maxn,self.dim])
        self.m = np.zeros(self.maxn,dtype=float)
        vrand = vmax * np.random.random([self.maxn,self.dim])
        vnorm = np.zeros(vrand.shape)
        vnorm[:,0] = vrand[:,0] - vrand[:,0].mean()
        vnorm[:,1] = vrand[:,1] - vrand[:,1].mean()
        vnorm[:,2] = vrand[:,2] - vrand[:,2].mean()
        self.v = vnorm
        self.rdot = np.zeros(self.r.shape)
        self.vdot = np.zeros(self.v.shape)
        self.mdot = np.zeros(self.m.shape)
    
        if rinit == 'random':
            self.r[0:n,:] = configuration.random3d(n,0,xmax,0,ymax,0,zmax)
            self.m[:] = mass 

        if rinit == 'grid':
            if self.dim == 3:
                self.r[0:n,:] = configuration.grid3d(
                    n,side,(xmax/2.,ymax/2.,zmax/2.),spacing=spacing)
                self.m[:] = mass 
            else:
                #def grid(n,xside,yside,origin,spacing=1.0):
                self.r[0:n,:] = configuration.grid(
                    n,side[0],side[1],(xmax/2.,ymax/2.),spacing=spacing)
                self.m[:] = mass 
        elif rinit == 'fcc':
            self.r[0:n,:] = configuration.fcc3d(n,side,
                #(xmax/2.,ymax/2.,zmax/2.)
                (0.0,0.0,0.0),
                spacing=spacing)
            self.m[:] = mass 

        elif rinit == 'twophase':
            # If twophase is selected there must be a dict
            # called twophase with the parameters
            # twophase3d(n,nh,side,sideh,origin,spacingh=0.7,spacing=2.0):
            nh = twophase['nh']
            sideh = twophase['sideh']
            spacingh = twophase['spacingh']
            self.r[0:n,:] = configuration.twophase3d(n,nh,side,sideh,
                (xmax/4.,ymax/4.,zmax/4.)
                ,spacing=spacing,spacingh=spacingh)
            self.m[:] = mass

        #elif rinit == 'disorder':
            # A relaxed disorder with no very close particles
        #    self.r[0:n,:] = configuration.twophase3d(n,nh,side,sideh,
        #        (xmax/4.,ymax/4.,zmax/4.)
        #        ,spacing=spacing,spacingh=spacingh)

        elif rinit == 'load':
            print 'LOADING',n
            # Load the configuration from a target.
            # Today the target is hard-coded
            # Make some assumptions about the source file.
            # source = os.environ.get('SPDATA') + '/' + source 
            # read_step(source,self,step='last')
            # This was commented because the read step method
            # assumes a smooth particle system.
            pass
        else:
            self.m[:] = mass

        # Initialise values
        self.colour = 1.0,0.0,0.0 

        # State vectors to pass to numerical integrators
        n_variables = 7
        self.x = np.zeros([n_variables,self.maxn])
        self.xdot = np.zeros([n_variables,self.maxn])

        self.nlists = []
        self.forces = []

        self.controllers = controllers
        for controller in self.controllers:
            controller.bind_particles(self)

        """ Variables for measuring performance. """
        self.timing = {}
        self.timing['force time'] = -1
        self.timing['deriv time'] = -1
        self.timing['pairsep time'] = -1
        self.timing['update time'] = -1
        self.timing['integrate time'] = -1

    def create_particle(self,r,v=(0.0,0.0,0.0)):
        """Adds a new particle to the system.
        """
        self.r[self.n] = r
        self.m[self.n] = self.m[self.n-1] 
        self.v[self.n] = v
        self.n = self.n+1
        self.rebuild_lists()

    def rebuild_lists(self):
        """ rebuilds all nlists """
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build()

    def update(self,dt):
        """ Update the particle system, using the
            neighbour list supplied.

        # To test
        # r1 = self.r[:,:].copy()
        # print self.r - r1
        """
        self.rebuild_lists()
        t = time()
        self.step(self.gather_state,self.derivatives, \
                       self.gather_derivatives,self.scatter_state,dt)
        self.timing['integrate time'] = time() - t
        self.box.apply(self)
        self.steps += 1

# Right now these are hard coded to 3d. I am still mulling over the
# best approach the being able to do 2 and 1d if I want to.

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.r[0:self.n,2]
        self.x[4,0:self.n] = self.v[0:self.n,0]
        self.x[5,0:self.n] = self.v[0:self.n,1]
        self.x[6,0:self.n] = self.v[0:self.n,2]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.r[0:self.n,2] = x[3,0:self.n]
        self.v[0:self.n:,0] = x[4,0:self.n]
        self.v[0:self.n:,1] = x[5,0:self.n]
        self.v[0:self.n:,2] = x[6,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = self.mdot[0:self.n] 
        self.xdot[1,0:self.n] = self.rdot[0:self.n,0]
        self.xdot[2,0:self.n] = self.rdot[0:self.n,1]
        self.xdot[3,0:self.n] = self.rdot[0:self.n,2]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[5,0:self.n] = self.vdot[0:self.n,1]
        self.xdot[6,0:self.n] = self.vdot[0:self.n,2]
        return self.xdot


class SmoothParticleSystem2d(ParticleSystem):
    """A particle system with additional properties to solve smooth
    particle equations of motion.

    TODO: Cleanup the sph properties. Currently it is calling
    3D only code.  Need to either write Fortran 2D versions or
    see what is already implemented in Cython.

    """

    def __init__(self,n,
            maxn=100,
            controllers=[],
            xmax=XMAX,
            ymax=YMAX,
            vmax=VMAX,
            rinit=None,
            side=None,
            mass=1.0,
            spacing=None,
            temperature=1.0,
            thermostat_temp=1.0,
            thermostat=False,
            hshort=1.0,
            hlong=3.0,
            source=None,
            integrator='ieuler',
            set_temperature=False,
            twophase=None,
            thermalk=5.0,
            simbox=None):
        ParticleSystem.__init__(self,n,
            xmax=xmax,
            ymax=ymax,
            rinit=rinit,
            side=side,
            mass=mass,
            source=source,
            spacing=spacing,
            twophase=twophase,
            integrator=integrator,
            vmax=vmax,
            simbox=simbox,
            maxn=maxn)

        self.dim = 2
        self.V_split = VSPLIT
        self.r_amalg = RAMAL

        """
        Number of Dynamical Variables


        SPH Properties
        --------------
        rho -- mass density
        rhodot -- time rate of change of mass density
        gradv --  spatial gradient of velocity
        t -- temperature
        h -- smoothing length
        p -- isotropic pressure (repulsive)
        pco -- isotropic pressure (cohesive)
        jq -- heat flux

        P -- pressure tensor

        timing -- a dictionary of average execution times for
                  particular subroutines.

        """

        self.rho = np.zeros(self.maxn)
        self.rho_lr = np.zeros(self.maxn)
        self.rhodot = np.zeros(self.rho.shape)
        self.gradv = np.zeros([self.maxn,self.dim,self.dim])
        self.grad_rho = np.zeros([self.maxn,self.dim])
        self.grad_rho_lr = np.zeros([self.maxn,self.dim])
        self.jq = np.zeros([self.maxn,self.dim])
        #thermal properties
        self.t = np.ones(self.maxn)
        self.u = np.ones(self.maxn,dtype=float)
        
        # Set temperature
        self.t[:] = temperature

        if (rinit == 'load'):
            source = os.environ.get('SPDATA') + '/' + source
            ncfile = netCDF4.Dataset(source,'r')
            time = ncfile.dimensions['timestep']
            u = ncfile.variables['internal_energy']
            t = ncfile.variables['temperature']
            rho = ncfile.variables['density']
            m = ncfile.variables['mass']
            i = len(time) - 1
            self.u[0:n] = u[i,0:n]
            self.t[0:n] = t[i,0:n]
            self.m[0:n] = m[i,0:n]
            # Check energy and temperature
            # These should be the same
            eos.adash=ncfile.ascl
            eos.bdash=ncfile.bscl
            eos.kbdash=ncfile.kbscl
            utest = eos.get_vdw_u2d(self.t[0:n],self.rho[0:n])
            ttest = eos.vdw_temp2d(self.rho[0:n],self.u[0:n])
            ttest2 = eos.vdw_temp2d(self.rho[0:n],utest)
            #print 'UTEST',utest
            #print 'TTEST',ttest
            #print 'TTEST',ttest2
            ncfile.close()
            read_step(source,self,step='last')
            # Ideally, we would
            # Update the temp using this internal energy
            # ideally the particle system has a reference
            # to an equation of state method
            # assuming that the temp will be updated here correctly
            # but this gets complicated because we need to compute
            # the density.
            #self.t = self.calc_temperature(u,self.t,rho)

        self.thermostat = bool(thermostat)
        print 'THERMOSTAT',thermostat

        self.udot = np.zeros(self.maxn,dtype=float)
        self.h = np.zeros(self.maxn)
        self.hlr = np.zeros(self.maxn)
        self.h[:] = hshort
        self.hlr[:] = hlong
        self.p = np.zeros([self.maxn])
        self.pco = np.zeros([self.maxn])
        self.P = np.zeros([self.maxn,self.dim,self.dim])
        self.thermostat_temp = thermostat_temp

        n_variables = 9
        self.x = np.zeros([n_variables,self.maxn])
        self.xdot = np.zeros([n_variables,self.maxn])
        self.timing['SPAM time'] = -1


    def split(self,i):
        """ Splits particle i into a ragtag collection of four
            particles. Why four? Because calculating the new
            positions is straightforward.

            The daughter particle is in the centre. The son particles
            are displaced slightly.

            See Feldman and Bonet, Dynamic refinement and boundary contact
            forces in SPH with applications in fluid flow problems.
            International Journal for Numerical Methods in Engineering.
            2007

            (2D only!)

        """
        alpha = 0.6
        eps = 2.6

        if self.n > self.maxn-3:
            print "cannot refine any further"
            return False
       
        # The son 
        self.m[i] = self.m[i] / 4.0
        #self.h[i] = self.h[i] * alpha

        # Daughter 1
        self.r[self.n] = self.r[i] + eps*np.array([0,1])
        self.m[self.n] = self.m[i] 
        self.v[self.n] = self.v[i]
        
        # Daughter 2
        self.r[self.n+1] = self.r[i] + eps*np.array([0.866025,-0.5])
        self.m[self.n+1] = self.m[i] 
        self.v[self.n+1] = self.v[i]
 
        # Daughter 3
        self.r[self.n+2] = self.r[i] + eps*np.array([-0.866025,-0.5])
        self.m[self.n+2] = self.m[i] 
        self.v[self.n+2] = self.v[i]
        
        self.n = self.n+3
        #print "There are now ",self.n,"particles"
        return True


    def check_refine(self):
        split = False
        for i in range(self.n):
            #V = self.m[i]/self.rho[i]
            #if V > self.V_split:
            if self.rho[i] < rhosplit:
                if verbose: print "V ",i," is ",V," - splitting"
                split=True
                self.split(i)        
            if split:
                for nl in self.nlists: 
                    nl.rebuild_list = True
              
    def amalgamate(self,i,j):
        """ Amalgamates particles i and j, merging them together to become
            one awesome robot with supernatural powers.
        """
        # conserve momentum
        self.v[i] = (self.v[i]*self.m[i]+self.v[j]*self.m[j])/ \
                    (self.m[i]+self.m[j])
        self.r[i] = (self.r[j] - self.r[i])/2 + self.r[j] 
        self.m[i] = self.m[i] + self.m[j]
        self.r[j] = self.r[self.n-1]
        self.v[j] = self.v[self.n-1]
        self.m[j] = self.m[self.n-1]
        self.n = self.n - 1

    def check_amalg(self,nl):
        for k in range(nl.nip):
            i = nl.iap[k,0]
            j = nl.iap[k,1]
            if nl.rij[k] < self.r_amalg:
                #print nl.rij[k]
                #print "amalgamating",i,j
                self.amalgamate(i,j)

    def rebuild_lists(self):
        """ rebuilds all nlists """
        
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build()
                nl.compress()

    def apply_thermostat(self,target_temp):
        """ Apply a scaling thermostat. """
        u_old =  self.u
        tav = self.t.mean()
        scale_factor = target_temp / tav
        self.t = self.t * scale_factor
        self.u =  eos.get_vdw_u(self.t,self.rho)
        uenv = self.u.sum() - u_old.sum()

    def update(self,dt):
        """ Update the particle system.
            The time() command seems to take about
            10e-6 seconds to 
            print is 10 microseconds
        """
        t1 = time()

        if SPLIT:
            self.check_refine()
        if AMALGAMATE:
            self.check_amalg(self.nl_default)

        t = time()
        self.rebuild_lists()
        self.timing['nlist rebuild time'] = time() - t

        # Is this derivative step required?
        #t = time()
        #self.derivatives()
        #self.timing['deriv time'] = time() - t
       
        t = time()
        self.step(self.gather_state,self.derivatives, \
            self.gather_derivatives,self.scatter_state,dt)
        self.timing['integrate time'] = time() - t
        
        self.box.apply(self)
        self.t = eos.vdw_temp(self.rho,self.u)

        if self.thermostat:
            self.apply_thermostat(self.thermostat_temp)
        
        self.timing['update time'] = time() - t1
        self.steps += 1

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.v[0:self.n,0]
        self.x[4,0:self.n] = self.v[0:self.n,1]
        self.x[5,0:self.n] = self.rho[0:self.n]
        self.x[6,0:self.n] = self.p[0:self.n]
        # added second component of pressure
        self.x[7,0:self.n] = self.pco[0:self.n]
        self.x[8,0:self.n] = self.u[0:self.n]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.v[0:self.n:,0] = x[3,0:self.n]
        self.v[0:self.n:,1] = x[4,0:self.n]
        self.rho[0:self.n] = x[5,0:self.n]
        self.p[0:self.n] = x[6,0:self.n]
        self.pco[0:self.n] = x[7,0:self.n]
        self.u[0:self.n] = x[8,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = self.mdot[0:self.n] 
        self.xdot[1,0:self.n] = self.rdot[0:self.n,0]
        self.xdot[2,0:self.n] = self.rdot[0:self.n,1]
        self.xdot[3,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,1]
        self.xdot[5,0:self.n] = self.rhodot[0:self.n] 
        self.xdot[6,0:self.n] = 0
        self.xdot[7,0:self.n] = 0
        self.xdot[8,0:self.n] = self.udot[0:self.n]
        return self.xdot

    def set_vdw_properties(self,abk):
        #properties.ADASH = abk[0]
        #properties.BDASH = abk[1]
        #properties.KBDASH = abk[2]
        #properties.adash = abk[0]
        #properties.bdash = abk[1]
        #properties.kbdash = abk[2]
        #properties.set_vdw_props(abk[0],abk[1],abk[2])
        eos.adash=abk[0]
        eos.bdash=abk[1]
        eos.kbdash=abk[2]
        eos.ADASH=abk[0]
        eos.BDASH=abk[1]
        eos.KBDASH=abk[2]

    def derivatives(self):
        """ get the rate of change of each variable 
            for every particle 
        """
        self.rdot = self.v
        self.vdot[:,:] = 0.0
        self.udot[:] = 0.0

        t = time()
        for nl in self.nlists:
            nl.compress()
            nl.separations()
        self.timing['pairsep time'] = (time() - t)

        t = time()
        properties.spam_properties2d(self,self.nl_default)
        self.timing['SPAM time'] = time() - t
        
        t = time()
        for force in self.forces:
            force.apply()
        self.timing['force time'] = time() - t
        
        if ADVECTIVE:
            self.rdot[:,:] = 0.0


class SmoothParticleSystem(ParticleSystem3D):
    """A particle system with additional properties to solve smooth
    particle equations of motion.

    A long-range/short range smooth particle system.
    I'd like to make the way this is structured more elegant - would
    like to have a simple smooth particle system, and have the long-short
    one inherit from it. Another solution is to have a smooth properties
    class and let the smooth particle system have multiple instances.

    As I don't know the best solution I am going to think on it.

    """

    def __init__(self,n,
            d=3,
            maxn=100,
            controllers=[],
            xmax=XMAX,
            ymax=YMAX,
            zmax=ZMAX,
            vmax=VMAX,
            rinit=None,
            side=None,
            mass=1.0,
            spacing=None,
            temperature=1.0,
            thermostat_temp=1.0,
            thermostat=False,
            hshort=1.0,
            hlong=3.0,
            source=None,
            integrator='ieuler',
            set_temperature=False,
            twophase=None,
            thermalk=5.0,
            simbox=None):
        ParticleSystem3D.__init__(self,n,d=d,
            xmax=xmax,
            ymax=ymax,
            zmax=zmax,
            rinit=rinit,
            side=side,
            mass=mass,
            source=source,
            spacing=spacing,
            twophase=twophase,
            integrator=integrator,
            vmax=vmax,
            simbox=simbox,
            maxn=maxn)
        self.V_split = VSPLIT
        self.r_amalg = RAMAL

        """
        SPH Properties
        --------------
        rho -- mass density
        rhodot -- time rate of change of mass density
        gradv --  spatial gradient of velocity
        t -- temperature
        h -- smoothing length
        p -- isotropic pressure (repulsive)
        pco -- isotropic pressure (cohesive)
        jq -- heat flux

        P -- pressure tensor

        timing -- a dictionary of average execution times for
                  particular subroutines.

        """

        self.rho = np.zeros(self.maxn)
        self.rho_lr = np.zeros(self.maxn)
        self.rhodot = np.zeros(self.rho.shape)
        self.gradv = np.zeros([self.maxn,self.dim,self.dim])
        self.grad_rho = np.zeros([self.maxn,self.dim])
        self.grad_rho_lr = np.zeros([self.maxn,self.dim])
        self.jq = np.zeros([self.maxn,self.dim])
        #thermal properties
        self.t = np.ones(self.maxn)
        self.u = np.ones(self.maxn,dtype=float)
        
        # Set temperature
        self.t[:] = temperature

        if (rinit == 'load'):
            source = os.environ.get('SPDATA') + '/' + source
            ncfile = netCDF4.Dataset(source,'r')
            time = ncfile.dimensions['timestep']
            u = ncfile.variables['internal_energy']
            t = ncfile.variables['temperature']
            rho = ncfile.variables['density']
            m = ncfile.variables['mass']
            i = len(time) - 1
            self.u[0:n] = u[i,0:n]
            self.t[0:n] = t[i,0:n]
            self.m[0:n] = m[i,0:n]
            # Check energy and temperature
            # These should be the same
            eos.adash=ncfile.ascl
            eos.bdash=ncfile.bscl
            eos.kbdash=ncfile.kbscl
            utest = eos.get_vdw_u(self.t[0:n],self.rho[0:n])
            ttest = eos.vdw_temp(self.rho[0:n],self.u[0:n])
            ttest2 = eos.vdw_temp(self.rho[0:n],utest)
            #print 'UTEST',utest
            #print 'TTEST',ttest
            #print 'TTEST',ttest2
            ncfile.close()
            read_step(source,self,step='last')
            # Ideally, we would
            # Update the temp using this internal energy
            # ideally the particle system has a reference
            # to an equation of state method
            # assuming that the temp will be updated here correctly
            # but this gets complicated because we need to compute
            # the density.
            #self.t = self.calc_temperature(u,self.t,rho)
        elif (rinit == 'hotspot'):
            self.r[0:n,:],self.t[0:n] = \
            configuration.hotspotgrid3d(n,side,
                (xmax/2.,ymax/2.,zmax/2.),spacing=spacing)

        self.thermostat = bool(thermostat)
        print 'THERMOSTAT',thermostat

        self.udot = np.zeros(self.maxn,dtype=float)

        self.h = np.zeros(self.maxn)
        self.hlr = np.zeros(self.maxn)
        self.h[:] = hshort
        self.hlr[:] = hlong
        self.p = np.zeros([self.maxn])
        self.pco = np.zeros([self.maxn])
        self.P = np.zeros([self.maxn,self.dim,self.dim])
        self.thermostat_temp = thermostat_temp

        #for nl in self.nlists: 
        #    nl.build()
        #    properties.spam_properties(self,nl,nl.cutoff_radius)

        #properties.spam_properties(self,self.nlists[0],self.nlists[1] \
        #    ,self.h[0:self.n],self.hlr[0:self.n])

        n_variables = 11
        self.x = np.zeros([n_variables,self.maxn])
        self.xdot = np.zeros([n_variables,self.maxn])
        self.timing['SPAM time'] = -1

    def rebuild_lists(self):
        """ rebuilds all nlists """
        
        for nl in self.nlists: 
            if nl.rebuild_list:
                nl.build()
                nl.compress()

    def apply_thermostat(self,target_temp):
        """ Apply a scaling thermostat. """
        u_old =  self.u
        tav = self.t.mean()
        scale_factor = target_temp / tav
        self.t = self.t * scale_factor
        self.u =  eos.get_vdw_u(self.t,self.rho)
        uenv = self.u.sum() - u_old.sum()

    def update(self,dt):
        """ Update the particle system.
            The time() command seems to take about
            10e-6 seconds to 
            print is 10 microseconds
        """
        t1 = time()

        if SPLIT:
            self.check_refine()
        if AMALGAMATE:
            self.check_amalg(self.nl_default)

        t = time()
        self.rebuild_lists()
        self.timing['nlist rebuild time'] = time() - t

        # Is this derivative step required?
        #t = time()
        #self.derivatives()
        #self.timing['deriv time'] = time() - t
       
        t = time()
        self.step(self.gather_state,self.derivatives, \
            self.gather_derivatives,self.scatter_state,dt)
        self.timing['integrate time'] = time() - t
        
        self.box.apply(self)
        self.t = eos.vdw_temp(self.rho,self.u)

        if self.thermostat:
            self.apply_thermostat(self.thermostat_temp)
        
        self.timing['update time'] = time() - t1
        self.steps += 1

    def gather_state(self):
        """ Maps the particle system to a state vector for integration
        """
        self.x[0,0:self.n] = self.m[0:self.n]
        self.x[1,0:self.n] = self.r[0:self.n,0]
        self.x[2,0:self.n] = self.r[0:self.n,1]
        self.x[3,0:self.n] = self.r[0:self.n,2]
        self.x[4,0:self.n] = self.v[0:self.n,0]
        self.x[5,0:self.n] = self.v[0:self.n,1]
        self.x[6,0:self.n] = self.v[0:self.n,2]
        self.x[7,0:self.n] = self.rho[0:self.n]
        self.x[8,0:self.n] = self.p[0:self.n]
        # added second component of pressure
        self.x[9,0:self.n] = self.pco[0:self.n]
        self.x[10,0:self.n] = self.u[0:self.n]
        return(self.x)

    def scatter_state(self,x):
        """ Maps the state vector to a particle system
        """
        self.m[0:self.n] = x[0,0:self.n] 
        self.r[0:self.n,0] = x[1,0:self.n]
        self.r[0:self.n,1] = x[2,0:self.n]
        self.r[0:self.n,2] = x[3,0:self.n]
        self.v[0:self.n:,0] = x[4,0:self.n]
        self.v[0:self.n:,1] = x[5,0:self.n]
        self.v[0:self.n:,2] = x[6,0:self.n]
        self.rho[0:self.n] = x[7,0:self.n]
        self.p[0:self.n] = x[8,0:self.n]
        self.pco[0:self.n] = x[9,0:self.n]
        self.u[0:self.n] = x[10,0:self.n]

    def gather_derivatives(self):
        """ Maps particle system's derivatives to a state vector
        """
        self.xdot[0,0:self.n] = self.mdot[0:self.n] 
        self.xdot[1,0:self.n] = self.rdot[0:self.n,0]
        self.xdot[2,0:self.n] = self.rdot[0:self.n,1]
        self.xdot[3,0:self.n] = self.rdot[0:self.n,2]
        self.xdot[4,0:self.n] = self.vdot[0:self.n,0]
        self.xdot[5,0:self.n] = self.vdot[0:self.n,1]
        self.xdot[6,0:self.n] = self.vdot[0:self.n,2]
        self.xdot[7,0:self.n] = self.rhodot[0:self.n] 
        self.xdot[8,0:self.n] = 0
        self.xdot[9,0:self.n] = 0
        self.xdot[10,0:self.n] = self.udot[0:self.n]
        return self.xdot

    def set_vdw_properties(self,abk):
        fproperties.ADASH = abk[0]
        fproperties.BDASH = abk[1]
        fproperties.KBDASH = abk[2]
        fproperties.adash = abk[0]
        fproperties.bdash = abk[1]
        fproperties.kbdash = abk[2]
        fproperties.set_vdw_props(abk[0],abk[1],abk[2])
        eos.adash=abk[0]
        eos.bdash=abk[1]
        eos.kbdash=abk[2]
        eos.ADASH=abk[0]
        eos.BDASH=abk[1]
        eos.KBDASH=abk[2]
        #eos_set_vdw_params(abk[0],abk[1],abk[2])

    def derivatives(self):
        """ get the rate of change of each variable 
            for every particle 
        """
        self.rdot = self.v
        self.vdot[:,:] = 0.0
        self.udot[:] = 0.0

        t = time()
        for nl in self.nlists:
            nl.compress()
            nl.separations()
            #nl.apply_minimum_image()
        self.timing['pairsep time'] = (time() - t)

        t = time()
        # The 'spam_complete_force' implementation computes
        # all properties
        if SPROPS:
            fproperties.spam_properties(self,self.nl_default \
                ,self.h[0:self.n],self.hlr[0:self.n])
        self.timing['SPAM time'] = time() - t
        
        t = time()
        for force in self.forces:
            force.apply()
        self.timing['force time'] = time() - t
        
        if ADVECTIVE:
            self.rdot[:,:] = 0.0


#class IsokineticSmoothParticleSystem(SmoothParticleSystem):
#    """ In this system the particles don't move.
#        Another
#    pass
