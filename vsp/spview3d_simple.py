#! /usr/local/bin/python

""" Load a SPAM output file and take a look. Uses the axes3d module
    of matplotlib. (3D)
"""

import netCDF4 
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from pylab import figure, show, ion, clf, gcf
from splib.renderspam import get_grid_map3d
from splib.interpolate import splashmap_density_3dvol
import os
#ion()

if 'SPDATA' in os.environ.keys():
    BASEDIR = os.environ.get('SPDATA') + '/runs'
else:
    BASEDIR = '.'


class spview3d:
    """ Loads a smooth particle netcdf file.

        Initialisation
        rname is the name of the ouput file #without the extension
    """
    def __init__(self,rname
                ,basedir=BASEDIR
                ,fig='None'):
        fname = basedir + '/' + rname #+ '.nc'
        print 'Opening',fname
        self.fname = fname
        if fig == 'None':
            self.fig = figure()
            self.ax = axes3d.Axes3D(self.fig)
        else:
            self.fig = fig
            self.ax = axes3d.Axes3D(fig)

    def openfile(self):
        self.ncfile = netCDF4.Dataset(self.fname,'r')
        self.timev = self.ncfile.variables['timestep']
        self.n = self.ncfile.variables['particle']

    def closefile(self):
        self.ncfile.close()

    def save(self,name='spimg.png',basedir=BASEDIR):
        plt.savefig(basedir + '/' + name)

    def get_ntime(self):
        ncfile = netCDF4.Dataset(self.fname,'r')
        runtime = ncfile.dimensions['timestep']
        nt = runtime.__len__()
        ncfile.close()
        return nt

    def scatterstep(self,step):
        self.openfile()
        self.r = self.ncfile.variables['position']
        #self.v = self.ncfile.variables['velocity']
        #self.a = self.ncfile.variables['acceleration']
        x = self.r[step,:,:]
        n = self.n[:].max()
        ascl = 0.1
        vscl = 0.1
        if np.isnan(x).any():
            print 'NaN present'
        else:
            x = x.squeeze()
            ax = self.ax
        
        ax.plot(x[:,0],x[:,1],x[:,2],'o' )

        spat = self.ncfile.variables['spatial']
        xmax = spat.xmax
        ymax = spat.ymax
        zmax = spat.zmax
        ax.set_xlim3d(0,xmax)
        ax.set_ylim3d(0,ymax)
        ax.set_zlim3d(0,zmax)

        self.closefile()


    def clear(self):
        clf()
        self.ax = axes3d.Axes3D(gcf())

    def plotall(self):

        for i in range(10):
            x = r[i,:,:]
            if np.isnan(x).any():
                print 'NaN present'
                break
            else:
                x = x.squeeze()
                ax.scatter3D(x[:,0],x[:,1],x[:,2])
                sleep(2.0)

if __name__ == '__main__':
        print 'Opening file',sys.argv[1]
        infile = sys.argv[1]
        infile = 'grand_gas_more.nc'
        sp = spview3d(infile)
        ncfile = netCDF4.Dataset(BASEDIR + '/' + infile,'r')
        time = ncfile.variables['timestep']
        nt = np.max(time[:])
        Z = sp.scatterstep(0)

