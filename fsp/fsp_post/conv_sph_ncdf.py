#!/usr/bin/env python

""" 
    Converts the ASCII output from the SPAC sph code into a nice netCDF format.

"""

import numpy
import netCDF4 
import sphstate
import os
import glob
import sys
import load_sphvars

def main():
    #if( len(sys.argv) < 1):
    #   print "using default filename"
        #filename = "default.nc"
    if( len(sys.argv) < 2):
        print "using current path"
        path = os.path.abspath('.')
        name = os.path.basename(path)
        filename = name + ".nc"
    else:
        filename= sys.argv[1]+".nc"
    print "creating ",filename
    create_sph_ncfile(filename)

def create_sph_ncfile(filename):
    print filename
    nc_file = netCDF4.Dataset(filename,'w')

    #Set up attributes
    #attributes should cover most of the options in
    #the configuration file
   
    #old ofile option deprecated 
    #ofile = open('sphvars.var','r' )
    ifile = open('spinput.txt','r' )
    vars = load_sphvars.load_sphvars_new(ifile)

    # Miscellaneous attributes
    setattr(nc_file,'Date',1)
    setattr(nc_file,'Creator','ac')

    # sphvars file attributes
    setattr(nc_file,'Run',filename)
    setattr(nc_file,'Number',int(vars['NP']))
    setattr(nc_file,'Dimensions',int(vars['DIM']))
    setattr(nc_file,'Maxinteractions',int(vars['MAXIP']))
    setattr(nc_file,'nsteps',int(vars['TSTEPS']))
    setattr(nc_file,'Timestep_size',float(vars['DT']))
    setattr(nc_file,'xmax',float(vars['XMAX']))
    setattr(nc_file,'ymax',float(vars['YMAX']))
    setattr(nc_file,'snapshot_frequency',int(vars['SNAPFREQ']))
    setattr(nc_file,'cutoff_radius',float(vars['CUTOFF']))
    setattr(nc_file,'spacing',float(vars['SPACING']))
    setattr(nc_file,'box_side',float(vars['SIDE']))
    setattr(nc_file,'integrator',int(vars['INTEGRATOR']))
    
    setattr(nc_file,'art_visc',int(vars['ART_VISC']))
    setattr(nc_file,'entropy_type',int(vars['ENTROPY_TYPE']))
    setattr(nc_file,'conduction',int(vars['CONDUCTION']))
    setattr(nc_file,'core_size',float(vars['CORE_SIZE']))
    setattr(nc_file,'repulsion',float(vars['REPULSION']))
    setattr(nc_file,'restrict_v',int(vars['RESTRICT_V']))
    setattr(nc_file,'rho_eq',int(vars['RHO_EQ']))
    setattr(nc_file,'unstable_rho',int(vars['UNSTABLE_RHO']))
    
    setattr(nc_file,'eta',float(vars['ETA']))
    setattr(nc_file,'zeta',float(vars['ZETA']))
    setattr(nc_file,'sound_speed',float(vars['SOUND_SPEED']))
    
    setattr(nc_file,'thermostat_type',int(vars['THERMOSTAT_TYPE']))
    setattr(nc_file,'thermostat_temp',float(vars['THERMOSTAT_TEMP']))

    setattr(nc_file,'eos',int(vars['EOS']))
   
    setattr(nc_file,'adash',float(vars['ADASH']))
    setattr(nc_file,'bdash',float(vars['BDASH']))
    setattr(nc_file,'kbdash',float(vars['KBDASH']))
    setattr(nc_file,'cgrad',float(vars['CGRAD']))

    setattr(nc_file,'ndim',float(vars['NDIM']))
    setattr(nc_file,'origin',float(vars['ORIGIN']))

    setattr(nc_file,'tolerance',float(vars['TOLERANCE']))
    setattr(nc_file,'n_boundary_1',float(vars['NBP_1']))
    setattr(nc_file,'n_boundary_2',float(vars['NBP_2']))
    setattr(nc_file,'n_boundaries',float(vars['NBOUNDS']))
 
    # .. and so on. Perhaps there is a quick way to do
    # this automagically by parsing the sphvars.var file?

    # Also need to figure out what type each variable is.

    # get the number of particles
    pop = int(vars['NP']) + int(vars['NBP_1']) + int(vars['NBP_2'])
    
    #pop = load_sphvars.getPop(ofile)

    # Create netcdf dimensions
    # number of particles
    # spatial dimensions
    # timestep number
    nc_file.createDimension('timestep',None)
    nc_file.createDimension('particle',pop)
    nc_file.createDimension('spatial',2)

    # Create variables for the dimensions, and populate them
    tstep = nc_file.createVariable('timestep','d',('timestep',))
    part = nc_file.createVariable('particle','i',('particle',))
    space = nc_file.createVariable('spatial','i',('spatial',))
   
    part[:] = numpy.array(range(pop))
    space[:] = numpy.array([0,1])

    dimnames = nc_file.dimensions.keys()

    #Set up variables
    #every particle property has a variable
    #and there are also variables for the box size
    # (and later for the box dimensions)
    #a variable for 'time elapsed' at each step (for variable stepping)
    # see http://amber.scripps.edu/netcdf/nctraj.html  for inspiration
    #and some of the constants(?)

    #each variable needs a "units" attribute

    #vector variables
    v_dims =('timestep','particle','spatial')

    #scalar variables
    sc_dims = ('timestep','particle')
    
    #histogram variables
    hist_dims = ('timestep')

    #total and average variables
    tot_dims = ('timestep')

    r = nc_file.createVariable('position','d',v_dims)
    v = nc_file.createVariable('velocity','d',v_dims)
    a = nc_file.createVariable('acceleration','d',v_dims)
    temp = nc_file.createVariable('temperature','d',sc_dims)
    energy = nc_file.createVariable('internal_energy','d',sc_dims)
    mass = nc_file.createVariable('mass','d',sc_dims)
    rho = nc_file.createVariable('density','d',sc_dims)
    press = nc_file.createVariable('pressure','d',sc_dims)
    ss =nc_file.createVariable('sound_speed','d',sc_dims)
    visc =nc_file.createVariable('viscosity','d',sc_dims)
    h = nc_file.createVariable('smoothing_length','d',sc_dims)
    hl = nc_file.createVariable('long_smoothing_length','d',sc_dims)
    q = nc_file.createVariable('heat_flux','d',v_dims)
    vsm= nc_file.createVariable('smoothed_velocity','d',v_dims)
    psm =nc_file.createVariable('smoothed_pressure','d',sc_dims)
    tmpsm =nc_file.createVariable('smoothed_temperature','d',sc_dims)
    grad_rho = nc_file.createVariable('density_gradient','d',v_dims)
    ptype = nc_file.createVariable('particle_type','u1',sc_dims)

    #now set up the non-particle averaged or total system variables
    # kinetic energy, internal energy, isolated Hamiltonian

    V = nc_file.createVariable('total_kinetic_energy','d',tot_dims) 
    T = nc_file.createVariable('total_internal_energy','d',tot_dims)
    tav = nc_file.createVariable('average_temp','d',tot_dims)
    rhoav = nc_file.createVariable('rho_average','d',tot_dims)
    tstat_energy = nc_file.createVariable('thermostat_energy','d',tot_dims)
    TV =  nc_file.createVariable('hamiltonian','d',tot_dims)
    dti = nc_file.createVariable('dt','d',tot_dims)
    sys_dt = nc_file.createVariable('systime','d',tot_dims)
     
    #the path thing is an issue
    path = os.path.dirname(filename)
    print path
    #state_files = glob.glob(path + '/sphstate.????????')
    state_files = glob.glob('sphstate.????????')
    state_files.sort()
    if len(state_files) == 0:
        print 'no state files - no changes made'
        nc_file.close()
        exit(0)
    i = 0

    for file in state_files:
        tstep[i]=i
        print file
        data = sphstate.load_sphstate(file)
        print "loading "+ file    
        if r.shape[1] != data.shape[0]:
            print r.shape
            print data.shape
            print "oops! shape mismatch"
            exit(0)
        
        r[i,:,:] = data[:,0:2]
        v[i,:,:] = data[:,2:4]
        a[i,:,:] = data[:,4:6]
        temp[i,:] = data[:,6]
        energy[i,:] = data[:,7]
        mass[i,:] = data[:,8]
        rho[i,:] = data[:,9]
        press[i,:] = data[:,10]
        ss[i,:] = data[:,11]
        visc[i,:] = data[:,12]
        h[i,:] = data[:,13]
        q[i,:,:] = data[:,14:16]
        vsm[i,:,:] = data[:,16:18]
        psm[i,:] = data[:,18]
        tmpsm[i,:] = data[:,19]
        grad_rho[i,:,:] = data[:,20:22]
        ptype[i,:] = data[:,22]
        hl[i,:] = data[:,23]

        if data[:,9].shape != rho[i,:].shape:
            print "shape error"
        i=i+1
        
    nc_file.sync()
    
    tot_data = numpy.loadtxt("properties.output")
    indices = numpy.array(range(0,tot_data.shape[0],nc_file.snapshot_frequency))
    print 'number of indices',indices.shape
    print 'tot data shape',tot_data.shape
    print V.shape
    print indices
    if indices.shape[0] > V.shape[0]:
        print indices.shape[0],'gt',V.shape[0]
        indices = indices[0,0:V.shape[0]]
        print indices.shape
    V[:] = tot_data[indices,0]
    T[:] = tot_data[indices,1] 
    tav[:] = tot_data[indices,2]
    rhoav[:] = tot_data[indices,3]
    tstat_energy[:] = tot_data[indices,4]  
    TV[:] =  tot_data[indices,5]
    dti[:] = tot_data[indices,6]
    sys_dt[:] = tot_data[indices,7]
    nc_file.sync()
    nc_file.close()

if __name__ == "__main__":
    main()
