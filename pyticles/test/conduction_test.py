""" 
    Test that heat is conducted.
    hotspotgrid3d(n,side,origin,spacing=1.0,temp=(0.2,1.9)): 
"""

import sys
from time import time
import particles
import spam_complete_force
import neighbour_list
import f_properties as properties
from spam_nc import create_sph_ncfile, write_step, write_box
import numpy as np
import spkernel, spdensity
import forces
import c_forces
from box import PeriodicBox
from vsplib.c_interpolate import splash_map_3d
import feos
import matplotlib.pyplot as plt

def spawn_system(
        max_steps = 500,
        XMAX = 10,
        YMAX = 10,
        ZMAX = 10 ,
        NDIM = 3,
        NP = 125,
        SIDE = (5,5,5),
        VMAX = 0.0,
        dt = 0.05,
        SPACING = 0.8,
        TEMPERATURE = 0.95,
        HLONG = 4.0,
        HSHORT = 2.0,
        thermalk = 9.2e06,
        RINIT = 'hotspot',
        sfname = 'None',
        ascl = 7.45e+04,
        bscl = 5.84e-01,
        kbscl = 3.29e+04,
        pmass = 1.386e-01,
        cgrad = 8.35e02,
        thermostat = False,
        sigma = 0.0,
        rcoef = 0.0,
        eta = 0.0,
        zeta = 0.0,
        ofname = 'data/toybox.nc',
        write_frequency = 1,
        config = None,
        gravity = 0,
        gravk = 1.0
        ):
        """ Create a smooth particle system and integrate it forward in time.
            Returns the path of the netcdf file containing the results.
            config is a config object this will replace the other 
            variables soon.

        """

        #NP = SIDE[0]*SIDE[1]*SIDE[2]
        cnt = 0
        fps = 0
        print max_steps
        print "Initialising"
        box = PeriodicBox(xmax=XMAX,ymax=YMAX,zmax=ZMAX)
        p = particles.SmoothParticleSystem(
                NP,maxn=NP,
                d=3,
                rinit=RINIT,
                vmax=VMAX,
                side=SIDE,
                spacing=SPACING,
                xmax=XMAX,
                ymax=YMAX,
                zmax=ZMAX,
                source=sfname,
                temperature=TEMPERATURE,
                hlong=HLONG,
                hshort=HSHORT,
                thermostat_temp=TEMPERATURE,
                thermostat=False,
                mass=pmass,
                simbox=box
            )
        nl = neighbour_list.FastVerletList(p,cutoff=HLONG)
        nl.max_neighbours = 20
        p.nlists.append(nl)
        p.nl_default = nl
        p.forces.append(spam_complete_force.SpamHeat(p,nl))
        tstart = time() 

        nl.build()
        nl.separations()
        properties.ADASH = ascl
        properties.BDASH = bscl
        properties.KBDASH = kbscl
        properties.set_vdw_props(ascl,bscl,kbscl)
        p.set_vdw_properties((ascl,bscl,kbscl))
        properties.energy_from_temp(p,nl,p.h,p.hlr)
        properties.spam_properties(p,nl,p.h,p.hlr)

        # This is a crap way to set these parameters
        #p.set_vdw_properties((ascl,bscl,kbscl))
        print 'Built list and calc properties',time()-tstart
        cnt = 0
        attribs = {'creator':'Andrew', 'log':'functional test' , 'vdwa':ascl, 'vdwb':bscl}
        create_sph_ncfile(ofname,attribs,NP,NDIM,config=config)
        write_box(ofname,p)
        print "--------------------------------------------"
        print "Running"
        tstartrun = time()
        write_step(ofname,p)
        ui = p.u.copy()
        ti = p.t.copy()
        for i in range(max_steps):
            p.update(dt)
            #print time() - tstart, 'updating'
            
            #if i % write_frequency == 0:
            #    write_step(ofname,p)
            #    print time() - tstart,'nip',nl.nip,'writing step', i
            #print 'Step',i,'took',time()-tstart
            # ->-> What takes 118 seconds with this print statement
            
        print 'Completed',i,'steps, in',time()-tstartrun

        # We took the temperature and internal energy
        # of all particles earlier.
        # They are expected to approach a uniform temperature.
        deltau = p.u - ui
        deltat = p.t - ti
        print deltau
        print deltat
        return ofname,p

if __name__ == '__main__':
    ofname,p = spawn_system()



