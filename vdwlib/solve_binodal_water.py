#! /usr/local/bin/python
""" Numerically solve for the binodal points of the 
    van der Waals equation of state.

    OUTPUT
    ------
    binodes = np.load('binodal_data_water.npy')
        temperature,pressure,binodal1,homogeneous,binodal2

    Given the van der Waals equation in the volume, with pressure
    as the dependent variable:

        1. Compute the area under the curve between two candiate points (these
            are the binodals).
        2. Compute the area in the constant pressure rectangle bounded by
            those points. (for what pressure?)

    The next stage for this file is to plot the binodal points for
    a number of temperatures.

    [Brent1973]	Brent, R. P., Algorithms for Minimization Without Derivatives. 
    Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

    We use scale values of the van der Waals parameters for water, based on
    Bedaux [ref]

    Scaling variables are:

        xscl = 2.81nm
        tscl = 1.00ns
        mscl = 7500 atomic units (1.247e-23 kilogram)

        a = 7.446e04 (scaled)
        b = 5.842e-01 (scaled)
        kb = 3.285e04 (scaled)
        kb = 7.869e01 (scaled)

    Alternative values for the van der Waals parameters are:
        a) Nugent and Posch vdw values a=2.0,b=0.5,kb=1.
        b) Usual dimensionless vdw values,a=3.0,b=1./3.,kb=8./3.

    This program has been cloned to plot_water_binodal_t98 so that I can always
    make a simple plot.

    matplotlib 1.0.1 required

"""

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import scipy
from scipy import optimize
from scipy import integrate,interpolate
from numpy import cos, arange, ones, asarray, zeros, mat, array

# Plotting parameters (matplotlib)
params = {'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12}

#rcParams.update(params)
plt.ion()
fig = plt.gcf()
plt.clf()
fig.set_dpi(100)
fig.set_size_inches(6,4)

# Scaled values for water
a = 7.446e04
b = 5.842e-01   
kb = 3.285e04

# Singularity at v = b
singular_volume = b
vmin = singular_volume + 0.1
vmax = 100.0
# T -- temperature at which the Maxwell construction is plotted
T = 1.0
V = np.arange(vmin,vmax,vmax/1000.,dtype="float64")
import sys
# Temperatures at which binodal points are plotted
#temperatures = [0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.05,1.1,1.12]
temperatures = [0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.05,1.1,1.12]
#temperatures = [1.0]

def vdweos(v,T):
    """ Here v is the volume of one scaling mass. """
    return( ((T*kb)/(v - b))  - (a/(v**2)) )

def diff(f,x,y,dx):
    """ Calculates a finite difference derivative
        using dx as the interval.
        F is a function of two variables because that's what
        I'm working with at the moment.
    """
    d = (f(x+dx,y) - f(x-dx,y) ) / (2.*dx) 
    return d

def calc_critical_temperature():
    return ( (8*a) / (27*b) ) / kb

def calc_spinodal_points(v,t):
    """ Returns the spinodal points as a 2 tuple.
        Find the roots of the gradient.
        This will fail for the long ends of the tail where the gradient looks
        like zero everywhere.
    """
    d = np.zeros(v.shape)
    for i in range(v.size):
        d[i] = diff(vdweos,v[i],t,0.001)
    spline = interpolate.splrep(v,d,k=3,s=2)
    roots = interpolate.sproot(spline,mest=3)
    if len(roots) != 2:
        #print 'Invalid roots'
        #print 'v=',v,'t=',t
        #plt.plot(v,d,'r')
        #plt.plot(v,vdweos(v,t),'b')
        #print t
        #pmax = p.max()
        #plot([roots[0],roots[1]],[pmax,pmax])
        #sys.exit()
        #p = vdweos(v,t)
        #plt.plot(v,p)
        return None
    else:
        a,b = roots
        return a,b

def area_function(a,b,t):
    area = integrate.quad(lambda v: vdweos(v,t),a,b)
    return area[0]

def const_pressure_area(a,b,p):
    """ Returns the area of the constant pressure
        rectangle between a and b."""
    return (b-a)*p

def solve_vdw_volume(p,t,V):
    """ Given a pressure, solve the vdw equation for the volumes """
    """ Assuming the temperature is sub critical.
    """
    #global vmax
    sa,sb = calc_spinodal_points(V,t)
    vdw_ct = lambda v: vdweos(v,t) - p
    # Find roots seperately
    pvals = vdw_ct(V)
    start = np.min(V[ (pvals > 0) & (V < sa)])
    result1 = optimize.brentq(vdw_ct,start,sa)
    result2 = optimize.brentq(vdw_ct,sa,sb)
    try:
        result3 = optimize.brentq(vdw_ct,sb,vmax)
    except ValueError:
        print 'Solve for volume failed...'
        print 'candidate pressure',p
        print 'interval',sb,vmax
        print 'psb',vdweos(sb,t)
        print 'pvmax',vdweos(vmax,t)
        sys.exit()
    return result1,result2,result3

def binodal_area_test(p,v,t):
    """ Compute the area given p,v,t
    """
    # First find the spinodal points
    sa,sb = calc_spinodal_points(v,t)
    # Now find the coexisting volumes for this pressure
    b1,b2,b3 = solve_vdw_volume(p,t,v)
    ba = min(b1,b2,b3)
    bb = max(b1,b2,b3)
    # Now compute the area
    area = integrate.quad(lambda v: vdweos(v,t),ba,bb)
    return area[0] - const_pressure_area(ba,bb,p)

def solve_binodes(V,T):
    """ Solve the van der waals equations for the binodal points at a
        given temperature.
    """
    sa,sb = calc_spinodal_points(V,T)
    psa = vdweos(sa,T)
    if psa < vdweos(vmax,T):
        psa = vdweos(vmax,T)
        psa = vdweos(vmax,T)
    psb = vdweos(sb,T)
    areafunc = lambda p: binodal_area_test(p,V,T)
    coex_pres = optimize.brentq(areafunc,psa,psb)
    bv1, bv2, bv3 = solve_vdw_volume(coex_pres,T,V)
    return coex_pres,bv1,bv2,bv3

def compute_binodes(v,T):
    """ Compute and save as temperature,pressure,binodal1,homogeneous,binodal2 """
    nt =len(T) 
    binodes = np.zeros([nt,5])
    for i in range(nt):
        t = T[i]
        print 'finding binodes for T=',t
        p, bv1, bv2, bv3 = solve_binodes(v,t)
        binodes[i,0:5] = (t,p,bv1,bv2,bv3) 
    np.save('data/binodal_data_water.npy',binodes)

if __name__ == '__main__':
    # 1. Solve for binodal points for the set of temperatures
    from os.path import exists
    if not exists('data/binodal_data_water.npy'):
        compute_binodes(V,temperatures)

    # Solve for the temperature at which the maxwell construction is plotted
    sa,sb = calc_spinodal_points(V,T)
    # The pressure at the spinodal points
    psa = vdweos(sa,T)
    psb = vdweos(sb,T)
    print 'Spinodal points',sa,sb
    print 'Spinodal pressures',psa,psb

    # Find the coexisting pressure by a brentq optimisation between these two
    # pressures
    areafunc = lambda p: binodal_area_test(p,V,T)
    if psa < vdweos(vmax,T):
        psa = vdweos(vmax,T)
    coex_pres = optimize.brentq(areafunc,psa,psb)
    print 'Coexisting pressures=',coex_pres
    print 'Area under curve - constant pressure area=',areafunc(coex_pres)

    # Now solve for the volume given the pressure
    bv1, bv2, bv3 = solve_vdw_volume(coex_pres,T,V)

    # Generate the pressure/volume curve
    p = vdweos(V,T)

    # Main plotting section
    # ---------------------

    # Set axis limits
    pmax = 2 * np.max(p[ (V>sa) & (V<sb)] )
    plt.axis([vmin,vmax,0,pmax])
    #ax = plt.axes([vmin,vmax,0,pmax])
    ax = plt.gca()
    ax.set_autoscale_on(False)
    ax.set_xscale('log')
    plt.ylim(0,pmax)
    plt.xlim(vmin,vmax)
    ax.grid(True)

    # Plot binodal points
    binodes = np.load('data/binodal_data_water.npy')
    nt = binodes.shape[0]
    for i in range(nt):
        pi = binodes[i,1]
        ti = binodes[i,0]
        #ax2.plot([binodes[i,2],binodes[i,3]],[ti,ti],'go')
        ax.plot([binodes[i,2],binodes[i,4]],[pi,pi],'go')

    #ax.yaxis.set_ticks_position('left')
    #ax.yaxis.set_label_position('left')

    # Plot the vanderwaals equation of state
    ax.plot(V,p,label="vdw equation of state")
    ax.plot(V,np.zeros(V.size),'k')
    plt.xlabel("Volume/22.2 nm3")
    plt.ylabel("Pressure/0.44 Bar")
    plt.title("Binodal line, van der Waals Equation, Water")

    # Plot spinodal points
    ax.plot([sa,sb],[psa,psb],'ko')

    # Plot critical point
    tcrit = calc_critical_temperature()
    pcrit = vdweos(3*b,tcrit)
    ax.plot([3*b,3*b],[pcrit,pcrit],'ko')

    region = np.zeros(V.shape,dtype='b')
    # Plot binodal points
    # bv1,2,3 are binodal volumes
    ax.plot([bv1,bv2,bv3],[coex_pres,coex_pres,coex_pres],'ro')
    coexistence_pressure = np.repeat(coex_pres,V.shape)
    region[:]=True
    region[ (V<bv1) | (V>bv3) ] = False
    ax.fill_between(V,coexistence_pressure,p,where=region)

    #ax.plot(binodes[:,3],binodes[:,1],'ko')
    #plt.text(0.72,0.8,'xscl=2.81nm',transform=ax.transAxes)
    #plt.text(0.72,0.75,'tscl=1.00ns',transform=ax.transAxes)
    #plt.text(0.72,0.7,'mscl=7500u',transform=ax.transAxes)
    #plt.text(0.72,0.65,'Tscl=562K',transform=ax.transAxes)

    ax2 = plt.twinx()
    plt.ylim(0,pmax)
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    bix = (binodes[:,0] == 0.6) | (binodes[:,0] == 0.7) | (binodes[:,0] == 0.8) | ( binodes[:,0] == 1.1) |(binodes[:,0] == 0.9) 
    ax2.set_yticks(binodes[bix,1])
    ax2.set_yticklabels(map(str,binodes[bix,0]))
    plt.ylabel("Temperature/562K")

    plt.subplots_adjust(left=0.16,bottom=0.13,right=0.9,top=0.87)
    plt.savefig("img/vdw_water.png",format='png')


