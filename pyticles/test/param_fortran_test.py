#! /usr/local/bin/python

""" 
    Testing that the parameters passed at this level are communicated to the
    fortran code.

    Problems with inconsistent parameters (usually visible as
    temperature/energy) are traced to a failure to set required
    variables such as density.

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
import eos

# Here are the parameters for a particle system
# This is so ugly, there is a better Python pattern for this.
# TODO: get this from config_library
# TODO: use the config dictionary

max_steps = 5
XMAX = 10
YMAX = 10
ZMAX = 10
NDIM = 3
NP = 125
SIDE = (5,5,5)
VMAX = 0.0
dt = 0.05
SPACING = 1.0
TEMPERATURE = 0.95
HLONG = 4.0
HSHORT = 2.0
thermalk = 1.0
RINIT = 'grid'
sfname = 'None'
ascl = 7.45e+04
bscl = 5.84e-01
kbscl = 3.29e+04
pmass = 1.386e-01
cgrad = 0.0
thermostat = True
set_temperature = True
sigma = 0.0
rcoef = 0.0
eta = 0.0
zeta = 0.0
ofname = 'data/toybox.nc'
write_frequency = 1
config = None
gravity = 0
gravk = 1.0


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
        thermostat=False,
        mass=pmass,
        simbox=box
    )

nl = neighbour_list.FastVerletList(p,cutoff=HLONG)
nl.max_neighbours = 20
p.nlists.append(nl)
p.nl_default = nl

spamforce = spam_complete_force.SpamComplete(
    p,nl,adash=ascl,bdash=bscl,kbdash=kbscl,cgrad=cgrad,eta=eta,
    zeta=zeta,sigma=sigma,rcoef=rcoef)

p.forces.append(spamforce)

nl.build()
nl.separations()
p.set_vdw_properties((ascl,bscl,kbscl))

properties.density(p,nl,p.h,p.hlr)

# Need to get it started with the energy
import feos
feos.eos.adash = ascl
feos.eos.bdash = bscl
feos.eos.kbdash = kbscl
feos.eos.eos_set_vdw_params(ascl,bscl,kbscl)
feos.eos.calc_vdw_energy(p.u,p.t,p.rho)

u = {}
t = {}
u['start'] = p.u
t['start'] = p.t
feos.eos.calc_vdw_temp(p.u,t['start'],p.rho)

# Now the properties object
properties.ADASH = ascl
properties.BDASH = bscl
properties.KBDASH = kbscl
properties.set_vdw_props(ascl,bscl,kbscl)
properties.print_vdw_props()
properties.spam_properties(p,nl,p.h,p.hlr)

u['prop'] = p.u
t['prop'] = p.t

p.update(0.000001)

u['update'] = p.u
t['update'] = p.t

# Now that we have a particle system and have set some parameters,
# perform some calculations to determine if the parameters are being
# set correctly

# Compute energy and temperature directly

u['fortran'] = np.zeros(u['update'].shape)
t['fortran'] = np.zeros(t['update'].shape)
feos.eos.calc_vdw_temp(u['update'],t['fortran'],p.rho)
feos.eos.calc_vdw_energy(u['fortran'],t['fortran'],p.rho)

p.t[:] = TEMPERATURE
eos.adash = ascl
eos.bdash = bscl
eos.kbdash = kbscl
u['python'] = eos.get_vdw_u(p.t,p.rho)
t['python'] = eos.vdw_temp(p.rho,u['python'])
#feos.eos.calc_vdw_energy(p.u,p.t,rhoin)
#feos.eos.calc_vdw_temp(p.u,t,rhoin)
#p.t[:] = t

