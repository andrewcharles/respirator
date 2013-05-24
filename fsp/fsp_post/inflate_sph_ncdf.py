#!/usr/bin/python

""" Inflates the SPH netCDF file into the text output.
"""

import numpy
import netCDF4
import sphstate
import os
import glob
import sys


def main():
    if( len(sys.argv) <= 1):
        print "using current path"
        path = os.path.abspath('.')
        name = os.path.basename(path)
        filename = name + ".nc"
    else:
        filename= sys.argv[1]+".nc"
    inflate_sph_ncfile(filename)

def output_step(filename,t,name):
    """ 
        Writes one step from a SPAC netcdf file to
        the ASCII format that SPAC uses.

        f: netcdf file
        t: step to output
            if t=-1 get the last step
        name: filename to use
    """
    f = netCDF4.Dataset(filename,"r")
    time = f.dimensions['timestep']
    print time
    r = f.variables["position"]
    timesteps = r.shape[0]
    if t==-1:
        t = timesteps-1
    pop = r.shape[1]
    v = f.variables["velocity"]
    a = f.variables["acceleration"]
    temp = f.variables["temperature"]
    energy = f.variables["internal_energy"]
    mass = f.variables["mass"]
    rho = f.variables["density"]
    press = f.variables["pressure"]
    ss = f.variables["sound_speed"]
    visc = f.variables["viscosity"]
    h = f.variables["smoothing_length"]
    q = f.variables["heat_flux"]
    vsm = f.variables["smoothed_velocity"]
    psm = f.variables["smoothed_pressure"]
    tmpsm = f.variables["smoothed_temperature"]
    grad_rho = f.variables["density_gradient"]
    ofile = open(name,"w")
    for i in range(pop):
        ofile.write("%1.12e %1.12e " %(r[t,i,0],r[t,i,1]) )
	#print r[t,i,:]
        ofile.write("%1.12e %1.12e " %(v[t,i,0],v[t,i,1]) )
        ofile.write("%1.12e %1.12e " %(a[t,i,0],a[t,i,1]) )
        ofile.write("%1.12e " %temp[t,i])    
        ofile.write("%1.12e " %energy[t,i])
        ofile.write("%1.12e " %mass[t,i])
        ofile.write("%1.12e " %rho[t,i])
        ofile.write("%1.12e " %press[t,i])
        ofile.write("%1.12e " %ss[t,i])
        ofile.write("%1.12e " %visc[t,i])
        ofile.write("%1.12e " %h[t,i])
        ofile.write("%1.12e %1.12e " %(q[t,i,0],q[t,i,1]))
        ofile.write("%1.12e " %vsm[t,i,0])
        ofile.write("%1.12e " %vsm[t,i,1])
        ofile.write("%1.12e " %psm[t,i])
        ofile.write("%1.12e " %tmpsm[t,i])
        ofile.write("%1.12e " %grad_rho[t,i,0])
        ofile.write("%1.12e\n" %grad_rho[t,i,1] )
    f.close()


def inflate_sph_ncfile(filename):
    #this one takes the netcdf file and converts it into a bunch
    #of text documents
    #load the netcdf file
    f = netCDF4.Dataset(filename,"r")
    print f.dimensions.keys()
    print f.variables.keys()
    #Question:
    #check the order of the variables - is it the same
    #as when I wrote the file out?
    #Answer:
    #No, it seems to be random

    time = f.dimensions['timestep']
    print time.shape

    #Get all the variables
    r = f.variables["position"]
    print r.shape
    timesteps = r.shape[0]
    pop = r.shape[1]
    print r[1,2,:]

    v = f.variables["velocity"]
    a = f.variables["acceleration"]
    temp = f.variables["temperature"]
    energy = f.variables["internal_energy"]
    mass = f.variables["mass"]
    rho = f.variables["density"]
    press = f.variables["pressure"]
    ss = f.variables["sound_speed"]
    visc = f.variables["viscosity"]
    h = f.variables["smoothing_length"]
    q = f.variables["heat_flux"]
    vsm = f.variables["smoothed_velocity"]
    psm = f.variables["smoothed_pressure"]
    tmpsm = f.variables["smoothed_temperature"]
    grad_rho = f.variables["density_gradient"]

    # dummy
    timestep_size = 1

    # for each time
    for t in range(timesteps):
        #create a new sphstate file
        filename = "sphstate.%08d" %(timestep_size*t)
        print filename
        ofile = open(filename,"w")
        #write out this time's data to it
        for i in range(pop):
            ofile.write("%1.12e %1.12e " %(r[t,i,0],r[t,i,1]) )
            ofile.write("%1.12e %1.12e " %(v[t,i,0],v[t,i,1]) )
            ofile.write("%1.12e %1.12e " %(a[t,i,0],a[t,i,1]) )
            ofile.write("%1.12e " %temp[t,i])    
            ofile.write("%1.12e " %energy[t,i])
            ofile.write("%1.12e " %mass[t,i])
            ofile.write("%1.12e " %rho[t,i])
            ofile.write("%1.12e " %press[t,i])
            ofile.write("%1.12e " %ss[t,i])
            ofile.write("%1.12e " %visc[t,i])
            ofile.write("%1.12e " %h[t,i])
            ofile.write("%1.12e %1.12e " %(q[t,i,0],q[t,i,1]))
            ofile.write("%1.12e " %vsm[t,i,0])
            ofile.write("%1.12e " %vsm[t,i,1])
            ofile.write("%1.12e " %psm[t,i])
            ofile.write("%1.12e " %tmpsm[t,i])
            ofile.write("%1.12e " %grad_rho[t,i,0])
            ofile.write("%1.12e\n" %grad_rho[t,i,1] )
            



    f.close()

if __name__ == "__main__":
    main()
