#!/usr/bin/env python

""" Create a pyticles system, run it, then create a page with the output.
    Edit the spconfig parameters in this file to vary the system setup,
    run length, etc.

    We have a library of configs. Configs are specific to a model
    specification. These are just python class definitions.
    We also have a library of initial states.

    model spec / run scripts:
        - toybox.py
        - spbox.py
        - corebox.py

    Environment Variables 

"""

import sys
import os

from toybox import spawn_system
#from report3d import report
import os.path
import config_library

config_names = dir(config_library)

print 'Welcome to Smooth Particle Run Control'
#print config_names

if len(sys.argv) == 1:
    print 'Usage: rsprun.py run_name'
    print 'Usage: rsprun.py run_name config_name '
    print 'Doing a basic run '
    config_class = getattr(config_library,'basic')
    run_name = 'default_run'
    #print 'System Path: ',sys.path
    #print os.system('hostname')

if len(sys.argv) == 2:
    run_name = sys.argv[1]
    if sys.argv[1] in config_names:
        config_class = getattr(config_library,sys.argv[1])
    else:
        config_class = getattr(config_library,'basic')

if len(sys.argv) == 3:
    run_name = sys.argv[2]
    if sys.argv[1] in config_names:
        config_class = getattr(config_library,sys.argv[1])
    else:
        config_class = getattr(config_library,'basic')

spconfig = config_class()

# Make the data directory
data_base = 'data/runs'
if not os.path.exists(data_base):
    os.mkdir(data_base)

# The spawn_system function returns the name of the file with the data in it
toyrun = spawn_system(
        max_steps = spconfig.max_steps,
        XMAX = spconfig.XMAX,
        YMAX = spconfig.YMAX,
        ZMAX = spconfig.ZMAX,
        NDIM = spconfig.NDIM,
        NP = spconfig.NP,
        SIDE = spconfig.SIDE,
        VMAX = spconfig.VMAX,
        dt = spconfig.dt,
        SPACING = spconfig.SPACING,
        TEMPERATURE = spconfig.TEMPERATURE,
        HLONG = spconfig.HLONG,
        HSHORT = spconfig.HSHORT,
        RINIT = spconfig.RINIT,
        sfname = spconfig.sfname,
        ascl = spconfig.ascl,
        bscl = spconfig.bscl,
        kbscl = spconfig.kbscl,
        pmass = spconfig.pmass,
        cgrad = spconfig.cgrad,
        sigma = spconfig.sigma, 
        rcoef = spconfig.rcoef,
        terminate_dist = spconfig.terminate_dist,
        eta = spconfig.eta,
        zeta = spconfig.zeta,
        collide = bool(spconfig.collide),
        collide_dist = spconfig.collide_dist,
        gravity = spconfig.gravity,
        thermalk = spconfig.thermalk,
        avisc = spconfig.avisc,
        gravk = spconfig.gravk,
        set_temperature = bool(spconfig.set_temperature),
        thermostat = bool(spconfig.thermostat),
        write_frequency = spconfig.write_frequency,
        ofname = 'data/runs/'+run_name+'.nc',
        config = spconfig
        )

#toyrun = run_name+".nc"
#print 'Making report based on',toyrun
#report_base = 'reports'
#report_folder = report_base + '/' + run_name
#if not os.path.exists(report_folder):
#    os.mkdir(report_folder)
#report(rname=run_name,datafile='data/runs/'+toyrun,folder=report_folder)
os.system('./reportgen.py ' + run_name)
