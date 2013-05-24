""" This module contains generic numerical integration methods.
    Copyright (C) Andrew Charles 2008
    BSD License

    Most time is taken up in calc_derivs
    
"""
# scipy integrators have the interface
# X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
# where t is a sequence of time points
# 


def euler(get_state,calc_derivs,get_derivs,set_state,dt):
    """ Simple euler integration. 

        state -- vector describing the system's state.
        get_state -- function to obtain the system state
	    calc_derivs -- calculates the derivatives for the system
	    get_derivs -- gathers the derivates from the system into a vector
        scatter_state -- a function to scatter the state vector to a
                more descriptive objects (e.g. an n-body system organised in
                particles)
        dt -- timestep size

        These are usually functions of an object, e.g.

        euler(particle.get_state,particle.get_derivs,particle.set_state,0.1)
        
        typical application is to a particle system:
	
        def get_state():
          x[0,:] = p.r[:,1]
          x[1,:] = p.r[:,2]
          x[2,:] = p.v[:,1]
          ...etc
    """
    calc_derivs()
    x = get_state()
    xdot = get_derivs()
    x = x+xdot*dt
    set_state(x)

from time import time
def imp_euler(get_state,calc_derivs,get_derivs,set_state,dt):
    """ Improved euler integration (two step predictor-corrector
    """
    # First guess
    calc_derivs()
    x_start = get_state()
    xdot = get_derivs()
    c1 = xdot*dt
    # Take Euler step and get the second term
    x = x_start + c1
    set_state(x)
    calc_derivs()
    xdot = get_derivs()
    c2 = xdot * dt
    x = x_start + (c1+c2)/2.
    set_state(x)

def rk4(get_state,calc_derivs,get_derivs,set_state,dt):
    """ runge kutta fourth order integrator
    """

    # First guess 
    calc_derivs()
    x_start = get_state()
    xdot = get_derivs()
    c1 = xdot*dt

    # Take half step and get the second term
    x = x_start + c1 / 2.0
    set_state(x)
    calc_derivs()
    xdot = get_derivs()
    c2 = xdot * dt

    # Take another half step with c2 to get c3
    x = x_start + c2 / 2.0
    set_state(x)
    calc_derivs()
    xdot = get_derivs()
    c3 = xdot * dt

    # c4 step
    x = x_start + c3
    set_state(x)
    calc_derivs()
    xdot = get_derivs()
    c4 = xdot * dt

    # Final step
    x = x_start + (1.0/6.0) * (c1 + 2.*c2 + 2.*c3 + c4) 
    set_state(x) 
    
    del x_start
    del xdot
    del x
    del c1
    del c2
    del c3
    del c4












