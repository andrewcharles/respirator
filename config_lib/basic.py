""" Simple basic configurations. """

class basic(object):
    """ A low compute cost basic functional testing configuration. """
    max_steps = 100
    XMAX = 3
    YMAX = 3
    ZMAX = 3
    NDIM = 3
    SIDE = (3,3,3)
    NP = 27
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.66
    TEMPERATURE = 0.1
    set_temperature = 1
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    sfname = 'None',
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3
    cgrad = 0.0#8.35e02
    write_frequency = 5
    eta = 100.0
    zeta = 10.0
    thermalk = 9.2e06
    gravity = 0
    gravk = 100.0
    sigma = 0.7
    rcoef = -1.0e+04
    thermostat = 1
    collide = 1
    collide_dist = 0.4
    avisc = 1
    terminate_dist = 0.2


class tiny_gas(basic):
    """ Equilibrate a 125 particle gas in a confined box. """
    max_steps = 4000
    XMAX = 4.25
    YMAX = 4.25
    ZMAX = 4.25
    NDIM = 3
    SIDE = (5,5,5)
    NP = 125
    VMAX = 0.0
    dt = 0.0001
    SPACING = 0.8
    TEMPERATURE = 1.5
    HLONG = 4.0
    HSHORT = 2.0
    RINIT = 'grid'
    ascl = 7.45e+04
    bscl = 5.84e-01
    kbscl = 3.29e+04
    pmass = 0.3
    cgrad = 8.35e02
    write_frequency = 1
    eta = 100.0
    zeta = 10.0
    sigma = 1.5
    rcoef = -1000.0
    thermalk = 9.2e06
    collide = 1
    collide_dist = 0.4
    avisc = 1
    gravity = 0
    gravk = 0.0
