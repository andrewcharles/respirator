#! /usr/local/bin/python

""" 
    Testing fortran and python energy implementations (again?)

"""

import f_properties as properties
from spam_nc import create_sph_ncfile, write_step, write_box
import numpy as np
import spkernel, spdensity
import forces
import c_forces
from box import PeriodicBox
import eos

TEMPERATURE = 0.95
ascl = 7.45e+04
bscl = 5.84e-01
kbscl = 3.29e+04
N = 5

u = {}
t = {}

us = np.zeros([N])
ts = np.zeros([N])
rho = np.zeros([N])

ts[:] = TEMPERATURE
rho[:] = 1.0

# Inline calculation
energy = (3./2.)*TEMPERATURE*kbscl - ascl*1.0
print energy

energy2 = (3./2.)*TEMPERATURE*1.0 - 2.0*1.0
print energy2

# Need to get it started with the energy
import feos
feos.eos.adash = ascl
feos.eos.bdash = bscl
feos.eos.kbdash = kbscl
feos.eos.eos_set_vdw_params(ascl,bscl,kbscl)
feos.eos.calc_vdw_energy(us,ts,rho)

print us

uf = np.zeros([N])
feos.eos.calc_vdw_energy(uf,ts,rho)
print uf

eos.adash = ascl
eos.bdash = bscl
eos.kbdash = kbscl
up = eos.get_vdw_u(ts,rho)
print up
# Fail if up and us are different


# Temperature...
#tf = np.zeros([N])
#feos.eos.calc_vdw_temp(us,tf,rho)
#feos.eos.calc_vdw_temp(us,tf,rho)
#print tf
## Fail if ts and tf are not pretty close to each other
#tf[:] = TEMPERATURE
# Fail if tp and ts are different
#tp = eos.vdw_temp(rho,up)



