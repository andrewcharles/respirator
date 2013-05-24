""" Some equations of state for use with smooth particle simulations. """

adash=2
bdash=0.5
kbdash=1

adash = 7.446e04
bdash = 5.842e-01   
kbdash = 3.285e04

# Provided for backwards compatibility
#ADASH = adash #2.0
#BDASH = bdash #0.5
#KBDASH = kbdash #1.0
RHONAUGHT = 1.0

def ideal_isothermal(rho,t):
    """ Calculates the pressure from the kinetic
        equation of state. 
        Isothermal equation of state
    """
    return (rho * kbdash)

def art_water(rho,t):
    """ Equation of state. Isothermal, with a reference density.
        Models a compressible liquid.
        Artificial liquid state water.
    """
    return ((rho - RHONAUGHT)*kbdash)

def vdw_eos(rho,t):
    """ Given a density and an internal energy returns a
        temperature and a pressure.
    """
    return (rho*kbdash*t)/(1.-rho*bdash) - adash*rho*rho

def vdw(rho,t):
    """ Van der Waals repulsive, cohesive in a tuple.
    """
    #return vdw_eos(rho,t)
    return (rho*kbdash*t)/(1-rho*bdash), - adash*rho*rho

def get_vdw_u(T,rho):
   #return (3./2.)*T*kbdash - adash*rho
   return T*kbdash - adash*rho

def get_vdw_u2d(T,rho):
   return T*kbdash - adash*rho

def vdw_energy(rho,t):
    """ Returns van der waals energy. """
    return get_vdw_u(t,rho)
    #return t * KBDASH - ADASH * rho

def vdw_energy2d(rho,t):
    """ Returns van der waals energy. """
    return get_vdw_u2d(t,rho)

def vdw_temp(rho,u):
    return (u+adash*rho)/kbdash
    #return (2./3.)*(u+adash*rho)/kbdash

def vdw_temp2d(rho,u):
    return (u+adash*rho)/kbdash

#def vdw_temp(rho,u):
#    return (u + ADASH * rho)/KBDASH
                
# Separated van der waals equation of state
def vdw_hc(rho,t):
    return (rho*kbdash*t)/(1-rho*bdash) 

def vdw_co(rho,t):
    return - adash*rho*rho

def bha_eos(rho,t):
    """ Untested.  Barker-henderson-abraham equation of state. """
    f1 = -10.5045
    f2 =  20.7126
    f3 = -13.0723
    f4 =  2.50803
    g1 = -2.8466
    g2 =  10.9930
    g3 = -11.1515

    return rho * t \
        + rho * t * \
        ( (1 + f1 * rho + f2 * rho**2 + f3 * rho**3 + f4 * rho**4)
         /(1 + g1 * rho + g2 * rho**2 + g3 * rho**3) )

