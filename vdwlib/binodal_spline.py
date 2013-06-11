""" Using the water binodal data, get the coexisting
    densities for any temperature by fitting a spline
    and rooting it.

    Shows how to use scipy's spline capabilities to
    fit to data.

    Depends on data stored in

    binodal_data_water.npy
    0 Temperature
    1 Pressure
    2 Binodal densities for this pressure

    Change this data to change the output ...

    Andrew Charles 2012

"""

import numpy as np
from scipy import interpolate
from solve_binodal_water import vdweos
from scipy import optimize
import matplotlib.pyplot as plt
import os
import sys

DATA_DIR = os.environ.get('VSP_BASE') + '/data'
TCRIT = 1.15

def coex_pressure_spline(T):
    """ Well it looks like a spline of pressure against temperature
        maybe for analytical water?
    """
    binodal_points = np.load(DATA_DIR + '/binodal_data_water.npy')
    temperature = binodal_points[:,0]
    pressure = binodal_points[:,1]
    rdx = np.argsort(temperature)
    temperature = temperature[rdx]
    pressure = pressure[rdx]
    spline = interpolate.splrep(temperature,pressure,s=0.001,k=3)
    return interpolate.splev(T,spline)

def binodal_spline_coex_density(T):
    """ A spline fitted to analytical water data...
        This is not actually a spline, it returns the
        coexistence densities...
        (Formerly called binodal_spline_density
    """
    binodal_points = np.load(DATA_DIR + '/binodal_data_water.npy')
    binodal_densities = 1. / np.concatenate(
        (binodal_points[:,2],binodal_points[:,4]))
    binodal_temps = np.concatenate(
        (binodal_points[:,0],binodal_points[:,0]))
    binodal_temps = binodal_temps
    rdx = np.argsort(binodal_densities)
    binodal_densities = binodal_densities[rdx]
    binodal_temps = binodal_temps[rdx]

    bin_spline = interpolate.splrep(
        binodal_densities,binodal_temps,s=0.001,k=3)
    splinex = np.arange(-0.3,2.0,0.01)
    spliney = interpolate.splev(splinex,bin_spline,der=0)
    #plt.plot(splinex,spliney)
    #plt.plot([-0.3,2.0],[0.,0.0])

    spline = lambda v: interpolate.splev(v,bin_spline) - T

    # Find the critical temperature and return the vapour
    # density

    result1 = optimize.brentq(spline,-0.3,0.5)
    result2 = optimize.brentq(spline,0.5,2.0)
    roots = result2,result1
    rholiq,rhogas = roots
    return rholiq,rhogas

if __name__ == '__main__':
    print coex_pressure_spline(0.5)
    print binodal_spline_coex_density(float(sys.argv[1]))
