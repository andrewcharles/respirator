#! /usr/local/bin/python
""" This script generates a plot of the van der waals phase diagram.

    The binodal curve plot is a vexing problem and has its own
    file (binodal.py, solve_binodal.py) dedicated to it. 
   
    Values of X originally from ~/masters/active/vasp/pred_densities.txt
    These were computed somewhat laboriously using Maple to solve the 
    minimisation problem - laborious because the trial densities needed
    to be manually input for each temperature. There is of course a better
    way, I just don't know what it is.

    To plot the often used dimensionless vdw equation use
    kb = 8./3.
    b = 1./3.
    a = 2.0

    To plot the Nugent and Posch version use
    kb = 1.0
    b = 0.5
    a = 2.0

"""

import matplotlib
from pylab import *
import numpy as np
import scipy
import os
from vdwlib.binodal_spline import binodal_spline_coex_density
from vdwlib.solve_binodal import diff
from scipy import optimize
from scipy import integrate,interpolate
from numpy import cos, arange, ones, asarray, zeros, mat, array
a = 9
params = {'axes.labelsize': a,
          'text.fontsize': a,
          'legend.fontsize': a,
          'xtick.labelsize': a,
          'ytick.labelsize': a}

rcParams.update(params)


kb = 1.     
b = 0.5     
a =  2.0 

# b is the singular volume
# rho * b < 1
rhomax = (1.0 / b) - 0.01

rho = np.arange(0.0,rhomax,0.005,dtype="float")
V = 1.0 / rho[::-1] #np.arange(0.01,10.995,0.05,dtype="float")
T = np.arange(0,1.8,0.005,dtype="float")
BASE_DIR = os.environ.get("VSP_BASE")

def vdweos(rho,T):
    return( (rho*T*kb)/(1-(rho*b)) - (a*(rho**2)) )

def d_vdweos_drho(rho,T):
    return( -2*a*rho + (kb*T*(1-rho*b) - rho*kb*T*b)/(1 - rho*b)**2 )

def dvdw(rho,t):
    """ Computes a finite difference derivative for
        the vdw_eos function
    """
    return(diff(vdweos,rho,t,0.01) )

def dvdv(v,t):
    """ Computes a finite difference derivative.
        for two volumes
    """
    rho1 = 1./v[0]
    rho2 = 1./v[1]
    return (diff(vdweos,rho1,t,0.01),diff(vdweos,rho2,t,0.01)) 

def plot_vdw(rho,T):
    # Plot the vanderwaals equation of state
    plot(rho,vdweos(rho,T),label="vdw equation of state")
    d = np.zeros(rho.shape)
    for i in range(rho.size):
        d[i] = diff(vdweos,rho[i],T,0.005)
    plot(rho[1:rho.size-1],d[1:d.size-1],'r',label="dp/drho")
    plot(rho,np.zeros(rho.size),'k')
    ylim(-2.0,2.0)
    xlim(0,2.0)
    xlabel("density")
    ylabel("pressure")
    title("p(rho)")
    legend()

def plot_vdw_vol(v,T):
    # Plot the vanderwaals equation of state by volume
    #rho = np.arange(0.01,1.995,0.005,dtype="float")
    ax = gca()
    ax.set_xscale('log')
    p = vdweos(1./v,T)
    print p
    print p.max()
    plot(v,p,'b',label="vdw equation of state")
    d = np.zeros(v.shape)
    f = np.zeros(v.shape)
    for i in range(v.size):
        d[i] = -diff(vdweos,1./v[i],T,0.005)
    plot(v[1:],d[1:],'r',label="dp/dv")
    plot([0,100],[0,0],'k')
    ylim(-2.0,2.0)
    ax = gca()
    ax.set_xscale('log')
    xlim(0,100.0) 
    xlabel("Volume")
    title("p(V)")
    ylabel("Pressure")
    legend()

def spline_vdw(T):
    pressure = lambda T: vdweos(rho,T)
    spline = interpolate.splrep(rho,pressure(rho),s=0)
    splinex = rho
    spliney = interpolate.splev(splinex,spline,der=0)
    plot(splinex,spliney)
    print interpolate.sproot(spline)
    rho1 = 0.1
    rho2 = 1.8
    v1 = 1./rho1
    v2 = 1./rho2
    area = interpolate.splint(rho1,rho2,spline)
    print area
    print area - vdweos(rho2,T)*(v2 - v1) 


def rend_pressure():
    p = np.zeros([rho.size,T.size])
    for i in range(rho.size):
        for j in range(T.size):
            p[i,j] = vdweos(rho[i],T[j]) 
       
    p = p.transpose()
    cols = np.arange(-3,10,0.2)
    contourf(rho,T,p,cols,extend='both')
    xlabel('Density')
    ylabel('Temperature')
    #imshow(p.transpose())

def plot_spinodal():
    T = np.arange(0.1,1.2,0.1,dtype="float")
    rho = np.arange(0.0,1.99,0.001,dtype="float")
    d = np.zeros(rho.shape)

    n = size(T)
    ps = np.zeros([2*n,2])
    j=0

    for t in T:
        for i in range(rho.size):
            d[i] = diff(vdweos,rho[i],t,0.001)
            #d[i] = d_vdweos_drho(rho[i],t)
        spline = interpolate.splrep(rho,d,s=0)
        splinex = rho
        spliney = interpolate.splev(splinex,spline,der=0)
        axis([0,2,0,2])
        pt = interpolate.sproot(spline)

        if len(pt) == 2:
            #scatter((pt[0],),(t,))
            #scatter((pt[1],),(t,))
            ps[j,0] = pt[0]
            ps[j+1,0] = pt[1]
            ps[j,1] = t
            ps[j+1,1] = t
        j+=2

    psnz = ps[ps.nonzero()]
    num = size(psnz)
    psnz = psnz.reshape(num/2,2)
    psnz = psnz[psnz[:,0].argsort(),]
    #psnz = vstack((np.array([0.0,0.0]),psnz))
    #psnz = vstack((psnz,[2.0,0.0]))
    splineodal = interpolate.splrep(psnz[:,0],psnz[:,1])
    splinexodal = rho
    splineyodal = interpolate.splev(splinexodal,splineodal,der=0)
    axis([0,2,0,2])
    return plot(splinexodal,splineyodal,'r--')

def render_deriv():
    # plot derivative wrt rho
    dp = np.zeros([rho.size,T.size])
    for i in range(rho.size):
        for j in range(T.size):
            dp[i,j] = diff(vdweos,rho[i],T[j],0.001) 
       
    dp = dp.transpose()
    #imshow(p.transpose())
    plot_spinodal()
    contourf(rho,T,dp,[-420,0,4200])
    #colorbar()

def render_deriv_vol():
    # plot derivative wrt rho
    rho = np.arange(0.01,1.995,0.005,dtype="float")
    v = np.arange(0.01,10.995,0.05,dtype="float")
    dp = np.zeros([v.size,T.size])
    for i in range(v.size):
        for j in range(T.size):
            dp[i,j] = -diff(vdweos,1./v[i],T[j],0.01) 
    dp = dp.transpose()
    #imshow(p.transpose())
    contourf(v,T,dp,[-42,0,42])
    #colorbar()

def plot_state_point(V,M,T):
    #render_deriv()
    return plot_predicted(M/V,T)

def plot_predicted(x,y):
    """
        Plot the binodal curve based on analysis points
        with a dot where specified by x and y.
        Plot a line so we can read off the densities.

    """
    A,B = binodal_spline_coex_density(y)
    ax = gca()
    binodal_line = ax.plot(A,B,'b-')
    point = ax.scatter([x],[y],marker='o',c='b')
    temperature_line = ax.plot([min(rho),2.0],[y,y],'k-'
        ,label="Thermostat temperature") 
    axis([0,2,0,1.6])
    plot_spinodal()
    return binodal_line,temperature_line,point
    #ax.legend()

def calc_spinodal_points(t):
    """ Returns the spinodal points as a 2 tuple.
        Find the roots of the gradient.
    """
    # def diffpt(f,x,y,i):
    # print optimize.fsolve(lambda rho: diffpt(vdweos(rho,t),0.7 ))
    rho = np.arange(0.01,1.9,0.001,dtype="float")
    d = np.zeros(rho.shape)
    for i in range(rho.size):
        d[i] = diff(vdweos,rho[i],t,0.001) 
    spline = interpolate.splrep(rho,d,s=0)
    print "spinodals at",t,interpolate.sproot(spline)
    a,b =  interpolate.sproot(spline)
    return a,b 
     
def test():
    # These subroutines are working
    t = 0.5
    clf() 
    # The PV equation of state
    subplot(3,2,1)
    plot_vdw_vol(V[0:V.size],t)

    # The P-rho equation of state
    subplot(3,2,2)
    plot_vdw(rho,t)
   
    # Area under volume curve
    subplot(3,2,3)
    render_deriv_vol()
   
    subplot(3,2,4)
    render_deriv()

    subplot(3,2,5)
    
    subplot(3,2,6)
    rend_pressure()
    
    #plot_state_point(1.0,1.0,1.0)


    subplots_adjust(wspace=0.33,hspace=0.33)


if __name__ == "__main__":
    
    ion()
    test()
    #show()

