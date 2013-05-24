
""" A set of configurations for rsprun.py
    These configurations determine the number of particles,
    steps, thermostat, initial configuration and all other
    run parameters of the smoooth particle system.
    todo: no config should have negative start temperature

    See config_lib for the imported configurations

    No True or False: use 1 or 0
"""

from config_lib.basic import *
from config_lib.droplet_500 import *
from config_lib.quench_500 import *
from config_lib.gas_c import *
from config_lib.grand_quench import *
from config_lib.grand_droplet import *
from config_lib.grand_expand import *
from config_lib.grand_long import *
from config_lib.grand_coexist_cube import *
from config_lib.mini_phase import *
from config_lib.twog import *
from config_lib.threeg import *
from config_lib.testing import *

class tiny_fcc_relax(object):
    """ A low compute cost basic functional testing configuration. """
    max_steps = 100
    XMAX = 4
    YMAX = 3
    ZMAX = 3
    NDIM = 3
    SIDE = (4,3,3)
    NP = 4*3*3*4
    VMAX = 0.0
    dt = 0.0001
    SPACING = 1.0 
    TEMPERATURE = 0.5
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'fcc'
    sfname = 'None',
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3 #1.386e-01
    cgrad = 8.35e02
    write_frequency = 5
    eta = 0.0#100.0
    zeta = 0.0#10.0
    gravity = 0
    gravk = 0.0
    sigma = 1.0
    rcoef = 1.0#e+04


class liquid_box(basic):
    """ Equilibrate a 512 particle liquid in a confined box. """
    max_steps = 500
    XMAX = 7
    YMAX = 7
    ZMAX = 7
    NDIM = 3
    SIDE = (8,8,8)
    NP = 8*8*8
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.8
    TEMPERATURE = 0.7
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3
    cgrad = 0
    write_frequency = 10
    eta = 0.0
    zeta = 0.0

class tiny_grav(basic):
    """ Gravitational Stratification at macro scale"""
    gravity = 1.0
    max_steps = 1000
    write_frequency = 10
    XMAX = 4
    YMAX = 4
    ZMAX = 4
    NDIM = 3
    SIDE = (4,4,4)
    VMAX = 0.0
    SPACING = 1.0
    dt = 0.001
    ascl = 1.045e+03
    bscl = 1.038e-03
    kbscl = 4.615e+02
    gravk = 9.807
    cgrad = 9.244e-17
    eta = 5.0e-04
    zeta = 5.0e-05

class grav_strat(basic):
    """ Gravitational Stratification """
    max_steps = 10
    XMAX = 7
    YMAX = 7
    ZMAX = 7
    NDIM = 3
    SIDE = (7,7,7)
    NP = 7*7*7
    VMAX = 0.0
    dt = 0.0001
    SPACING = 1.0
    TEMPERATURE = 1.0
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3
    cgrad = 0
    write_frequency = 10
    eta = 0.0
    zeta = 0.0
    gravity = 1
    gravk = 10.0 #2.794e+14 

class tiny_droplet(tiny_gas):
    XMAX = 10.25
    YMAX = 10.25
    ZMAX = 10.25
    dt = 0.0001
    max_steps = 1000
    set_temperature = 1
    TEMPERATURE = 0.2

# alex
# ALEX

class bigpac(basic):
    """ A big one to see what vpac can do """
    max_steps = 3000
    XMAX = 20
    YMAX = 10
    ZMAX = 10
    NDIM = 3
    SIDE = (20,10,10)
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.7
    TEMPERATURE = 0.5
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 1.386e-01
    cgrad = 1.46e-40
    write_frequency = 10
    eta = 0.0
    zeta = 0.0

class middling(basic):
    """ A more resource hungy small job than basic """
    max_steps = 500
    XMAX = 5
    YMAX = 5
    ZMAX = 5
    NDIM = 3
    #SIDE = (20,10,10)
    SIDE = (4,4,4)
    VMAX = 0.0
    dt = 0.001
    SPACING = 1.0
    TEMPERATURE = 0.5
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 1.386e-01
    cgrad = 1.46e-40
    write_frequency = 1
    eta = 0.0
    zeta = 0.0




class bedaux_box_2_liquid(basic):
    """ One thousand particles.
        Equilibrate a liquid
    """
    max_steps = 4000
    XMAX = 15
    YMAX = 10
    ZMAX = 10
    NDIM = 3
    NP = 1500
    SIDE = (15,10,10)
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.9
    TEMPERATURE = 0.95
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 1.0
    write_frequency = 5 
    eta = 100.0
    zeta = 0.0

class bedaux_box_2_gas(basic):
    """ Equilibrate a gas
    """
    max_steps = 4000
    XMAX = 15
    YMAX = 10
    ZMAX = 10
    NDIM = 3
    NP = 1500
    SIDE = (15,10,10)
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.9
    TEMPERATURE = 0.95
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.2
    write_frequency = 5 
    eta = 100.0
    zeta = 0.0

class bedaux_combined(basic):
    max_steps = 3000
    write_frequency = 10
    NP = 3000
    XMAX = 30
    YMAX = 10
    ZMAX = 10
    NDIM = 3
    VMAX = 0.0
    RINIT = 'load'
    sfname = 'runs/bedaux_box_combine.nc'
    dt = 0.0001
    eta = 10.0
    zeta = 0.0
    TEMPERATURE = 0.95
    thermal = 1e06

class mini_box_liquid(basic):
    """ 
        Equilibrate a liquid
    """
    max_steps = 4000
    XMAX = 5
    YMAX = 5
    ZMAX = 5
    NP = 125
    NDIM = 3
    SIDE = (5,5,5)
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.9
    TEMPERATURE = 0.95
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 1.0
    write_frequency = 5 
    eta = 100.0
    zeta = 0.0

class mini_box_gas(basic):
    """ Equilibrate a gas
    """
    max_steps = 4000
    XMAX = 5
    YMAX = 5
    ZMAX = 5
    NDIM = 3
    NP = 125
    SIDE = (5,5,5)
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.9
    TEMPERATURE = 0.95
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.2
    write_frequency = 5 
    eta = 100.0
    zeta = 0.0

class mini_box_smash(basic):
    """ Test the collisions!
    """
    max_steps = 1000
    XMAX = 10
    YMAX = 10
    ZMAX = 10 
    NDIM = 3
    NP = 125
    SIDE = (5,5,5)
    VMAX = 0.0
    dt = 0.001
    SPACING = 1.3
    TEMPERATURE = 0.2
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3
    write_frequency = 5 
    eta = 100.0
    zeta = 0.0
    rcoef = -10.0

class mini_box_combined(basic):
    max_steps = 5000
    write_frequency = 10
    NP = 250
    XMAX = 10
    YMAX = 5
    ZMAX = 5
    NDIM = 3
    VMAX = 0.0
    RINIT = 'load'
    sfname = 'runs/mini_box_combine.nc'
    dt = 0.0001
    eta = 10.0
    zeta = 0.0
    TEMPERATURE = 0.95

class groovy(basic):
    """ A moderate compute cost, visually interesting setup. """
    max_steps = 1000
    XMAX = 5
    YMAX = 5
    ZMAX = 5
    NDIM = 3
    #SIDE = (20,10,10)
    SIDE = (6,6,6)
    VMAX = 0.0
    dt = 0.00005
    SPACING = 0.75
    TEMPERATURE = 0.2
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.5#1.386e-01
    cgrad = 8.35e02
    write_frequency = 5
    eta = 1.0
    zeta = 0.1

class gravity(object):
    """ A low compute cost basic functional testing configuration. """
    max_steps = 100
    XMAX = 4
    YMAX = 4
    ZMAX = 4
    NDIM = 3
    SIDE = (3,3,3)
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.75
    TEMPERATURE = 0.5
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    sfname = 'None',
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3 #1.386e-01
    cgrad = 8.35e02
    write_frequency = 5
    eta = 100.0
    zeta = 100.0

#class
#class threeg_fill_container(basic):
