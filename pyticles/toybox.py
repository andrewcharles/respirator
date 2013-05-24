#! /usr/local/bin/python

""" 
    A very small system used for functional testing.

    ALSO: the main driver program. It accepts a config structure.

    A 10x10x10 scaling length box
    Scaling length 2.81nm
    Scaling time 1ns
    Scaling mass 7520 au
    a = 7.45e+04
    b = 5.84e-01
    kb = 3.29e+04
    Number of particles = 6x6x6
    Target density = M / V
    V = 1000

    Particle mass 0.1386
    Temperature 0.95

    ## Copyright Andrew Charles, RMIT 2010                ###
    ## All rights reserved.                               ###

    todo: split into spawn system and run system.

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

def spawn_system(
        max_steps = 5,
        XMAX = 10,
        YMAX = 10,
        ZMAX = 10 ,
        NDIM = 3,
        NP = 125,
        SIDE = (5,5,5),
        VMAX = 0.0,
        dt = 0.001,
        SPACING = 1.0,
        TEMPERATURE = 0.95,
        HLONG = 4.0,
        HSHORT = 2.0,
        thermalk = 1.0,
        RINIT = 'grid',
        sfname = 'None',
        ascl = 7.45e+04,
        bscl = 5.84e-01,
        kbscl = 3.29e+04,
        pmass = 1.386e-01,
        cgrad = 0.0,
        thermostat = False,
        set_temperature = False,
        sigma = 0.0,
        rcoef = 0.0,
        avisc = True,
        eta = 0.0,
        zeta = 0.0,
        ofname = 'data/toybox.nc',
        write_frequency = 1,
        config = None,
        gravity = 0,
        gravk = 1.0,
        collide = 1.0,
        collide_dist = 0.6,
        terminate_dist = 0.2,
        ):
        """ Create a smooth particle system and integrate it forward in time.
            Returns the path of the netcdf file containing the results.
            config is a config object this will replace the other variables soon.
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
                set_temperature=set_temperature,
                thermostat=thermostat,
                mass=pmass,
                simbox=box
                #simbox=None
            )
        #nl = neighbour_list.VerletList(p,cutoff=HLONG)
        nl = neighbour_list.FastVerletList(p,cutoff=HLONG,tolerance=0.1)
        #nl = neighbour_list.KDList(p,cutoff=HLONG,tolerance=0.1)
        #nl.max_neighbours = 20
        p.nlists.append(nl)
        p.nl_default = nl
        p.forces.append(spam_complete_force.SpamComplete(
            p,nl,adash=ascl,bdash=bscl,kbdash=kbscl,cgrad=cgrad,eta=eta,
            zeta=zeta,sigma=sigma,rcoef=rcoef,thermalk=thermalk,
            art_viscosity=avisc))
        if collide:
            p.forces.append(forces.FortranCollisionForce(p,nl,
                cutoff=collide_dist))

        if gravity == 1:
            print 'GRAVITY'
            p.forces.append(forces.BodyForce3d(p,k=gravk,direction=[0,-1.0,0]))
        tstart = time() 

        nl.build()
        nl.separations()
        properties.set_vdw_props(ascl,bscl,kbscl)
        p.set_vdw_properties((ascl,bscl,kbscl))
        # TODO: This is not ideal
        properties.density(p,nl,p.h,p.hlr)
        properties.energy_from_temp(p,nl,p.h,p.hlr)
        properties.spam_properties(p,nl,p.h,p.hlr)
        # This is a crap way to set these parameters
        print 'Built list and calc properties',time()-tstart
        cnt = 0
        attribs = { 'creator':'Andrew', 
                    'log':'functional test',
                    'vdwa':ascl,
                    'vdwb':bscl}

        create_sph_ncfile(ofname,attribs,NP,NDIM,config=config)
        write_box(ofname,p)
        print "--------------------------------------------"
        print "STEP   INT  DERIV =  PAIR + SPAM +  FORCE   "
        tstartrun = time()
        print p.r
        write_step(ofname,p)
        for i in range(max_steps):
            
            # MAIN LOOP #

            tstart = time()
            p.update(dt)
            #print time() - tstart, 'updating'
            if np.isnan(p.r).any():
                print 'stopping due to nan'
                break
            if i % write_frequency == 0:
                write_step(ofname,p)
                print time() - tstart,'nip',nl.nip,'writing step', i

            if (nl.rij[0:nl.nip] < terminate_dist).any():
                print 'isclose'
                mind = nl.rij[0:nl.nip].min()
                print mind
                idx = np.where(nl.rij == mind)
                print 'pair',nl.iap[idx,:]
                pid1,pid2 = nl.iap[idx,0],nl.iap[idx,1]
                # print some more diagnostic information
                # velocity
                print p.rdot.shape
                print 'Velocity 1:',p.rdot[pid1,:]
                print 'Velocity 2:',p.rdot[pid2,:]
                # temperature
                # pressure
                # particle ids

                #break
            
            # Detailed profiling information
            g = p.timing.keys()
            g.sort()
            for k in g:
                # This overwrites the same line
                output = "\b" + k + " %4.3f" %(p.timing[k]) + ", "
                sys.stdout.write(output)
            #    sys.stdout.flush()
            # -- What takes 118 seconds with this print statement
            print 'Step',i,'took',time()-tstart

        print 'Completed',i,'steps, in',time()-tstartrun
        # -- Takes only 95 seconds without it
        return ofname,p

if __name__ == '__main__':
    spawn_system()


