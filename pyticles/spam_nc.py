#! /usr/local/bin/python

""" 
    Create and write a netcdf file of SPAM output.

    todo - would be nice to have a dimension for every timestep, and a
    dimension just for the frames we've snapshotted. should read up
    on how to dothis beofre making changes that might cause regression.
    It would have to be a separate time dimension.

    dimension time_fine
    dimension time_coarse

    position[time_coarse,spatial]
    system_energy[time_fine]

"""

import numpy
import netCDF4
import os
import glob
import sys

def create_sph_ncfile(filename,attribs,n,dim,config=None):
    """ Create a netCDF file for SPAM output data, roughly
        following the AMBER specification.
        see http://amber.scripps.edu/netcdf/nctraj.html  for inspiration
        and some of the constants(?)

        attribs -- a dictionary of attributes for the nc file name, value
        n -- number of particles
        dim -- spatial dimensions

    """
    print 'Creating ',filename
    if not os.path.exists('data/runs'):
        os.makedirs('data/runs')
    nc_file = netCDF4.Dataset(filename,'w')

    # Miscellaneous attributes
    #getattr and setattr are built in python functions
    setattr(nc_file,'Date',1)
    setattr(nc_file,'Creator','ac')

    # sphvars file attributes
    for name,val in attribs.iteritems():
        setattr(nc_file,name,val)

    if config is not None:
        if type(config) == dict:
            for name,val in config.iteritems():
                setattr(nc_file,name,val)
        else:
            attdic = dict( (name, getattr(config,name)) 
                for name in dir(config) if not name.startswith('__'))      
            for name,val in attdic.iteritems():
                setattr(nc_file,name,val)
    else:
        print 'No way'

    # Create netcdf dimensions
    # number of particles
    # spatial dimensions
    # timestep number
    nc_file.createDimension('timestep',None)
    nc_file.createDimension('particle',n)
    nc_file.createDimension('spatial',dim)

    # Create variables for the dimensions, and populate them
    tstep = nc_file.createVariable('timestep','d',('timestep',))
    part = nc_file.createVariable('particle','i',('particle',))
    space = nc_file.createVariable('spatial','i',('spatial',))
   
    part[:] = numpy.array(range(n))
    space[:] = numpy.array([0,1,2])

    dimnames = nc_file.dimensions.keys()

    # Set up variables
    # every particle property has a variable
    # and there are also variables for the box size
    # and the box dimensions.
    # a variable for 'time elapsed' at each step (for variable stepping)

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
    #temp = nc_file.createVariable('temperature','d',sc_dims)
    energy = nc_file.createVariable('internal_energy','d',sc_dims)
    mass = nc_file.createVariable('mass','d',sc_dims)
    rho = nc_file.createVariable('density','d',sc_dims)
    t = nc_file.createVariable('temperature','d',sc_dims)
    #vdwb = nc_file.createVariable('vdwa','d',tot_dims)
    #vdwa = nc_file.createVariable('vdwb','d',tot_dims)
    bounds = nc_file.createVariable('box','d','spatial')
    
    #rho_lr = nc_file.createVariable('density_lr','d',sc_dims)
    #jq = nc_file.createVariable('density_gradient','d',v_dims)
    #thermal properties
    #t = np.ones(self.maxn)
    #thermostat=thermostat

    #u = np.ones(self.maxn,dtype=float)
    #udot = np.zeros(self.maxn,dtype=float)

    #h = np.zeros(self.maxn)
    #hlr = np.zeros(self.maxn)
    #p = np.zeros([self.maxn])
    #pco = np.zeros([self.maxn])
    #P = np.zeros([self.maxn,self.dim,self.dim])
    #thermostat_temp = thermostat_temp

    #rho = 
    press = nc_file.createVariable('repulsive_pressure','d',sc_dims)
    press = nc_file.createVariable('cohesive_pressure','d',sc_dims)
    #ss =nc_file.createVariable('sound_speed','d',sc_dims)
    #visc =nc_file.createVariable('viscosity','d',sc_dims)
    #h = nc_file.createVariable('smoothing_length','d',sc_dims)
    #hl = nc_file.createVariable('long_smoothing_length','d',sc_dims)
    #q = nc_file.createVariable('heat_flux','d',v_dims)
    #vsm= nc_file.createVariable('smoothed_velocity','d',v_dims)
    #psm =nc_file.createVariable('smoothed_pressure','d',sc_dims)
    #tmpsm =nc_file.createVariable('smoothed_temperature','d',sc_dims)
    #grad_rho = nc_file.createVariable('density_gradient','d',v_dims)
    #ptype = nc_file.createVariable('particle_type','u1',sc_dims)

    #now set up the non-particle averaged or total system variables
    # kinetic energy, internal energy, isolated Hamiltonian

    #V = nc_file.createVariable('total_kinetic_energy','d',tot_dims) 
    #T = nc_file.createVariable('total_internal_energy','d',tot_dims)
    #tav = nc_file.createVariable('average_temp','d',tot_dims)
    #rhoav = nc_file.createVariable('rho_average','d',tot_dims)
    #tstat_energy = nc_file.createVariable('thermostat_energy','d',tot_dims)
    #TV =  nc_file.createVariable('hamiltonian','d',tot_dims)
    #dti = nc_file.createVariable('dt','d',tot_dims)
    #sys_dt = nc_file.createVariable('systime','d',tot_dims)
     
    nc_file.sync()
    nc_file.close()

def write_box(filename,p):
    b = p.box
    ncfile = netCDF4.Dataset(filename,'r+')
    box = ncfile.variables['box']
    box[0],box[1],box[2] = b.xmax,b.ymax,b.zmax
    spat = ncfile.variables['spatial']
    spat.xmax = b.xmax
    spat.ymax = b.ymax
    spat.zmax = b.zmax
    ncfile.sync()
    ncfile.close()

def write_step(filename,p):
    """ Write out this step to the output file. 
    
        p -- particle object
    
    """
    ncfile = netCDF4.Dataset(filename,'r+')
    time = ncfile.dimensions['timestep']
    i = len(time)
    timev = ncfile.variables['timestep']
    timev[i] = i+1

    r = ncfile.variables['position']
    v = ncfile.variables['velocity']
    a = ncfile.variables['acceleration']
    u = ncfile.variables['internal_energy']
    m = ncfile.variables['mass']
    t = ncfile.variables['temperature']
    rho = ncfile.variables['density']
    pre = ncfile.variables['repulsive_pressure']
    pco = ncfile.variables['cohesive_pressure']
    n = p.n
        
    r[i,:,:] = p.r[0:n,:]
    v[i,:,:] = p.v[0:n,:]
    a[i,:,:] = p.vdot[0:n,:]
    u[i,:] = p.u[0:n]
    m[i,:] = p.m[0:n]
    t[i,:] = p.t[0:n]
    rho[i,:] = p.rho[0:n]
    pre[i,:] = p.p[0:n]
    pco[i,:] = p.pco[0:n]

    ncfile.sync()
    ncfile.close()

def particle_from_file(filename,step='last'):
    """ Return a particle object from a file step. """

    ncfile = netCDF4.Dataset(filename,'r')
    # Create the particle object
    npart = ncfile.dimensions['particle'].size
    p = particles.SmoothParticleSystem(
            npart,maxn=npart,
            d=3,
            rinit='load',
            xmax=ncfile.xmax,
            ymax=ncfile.YMAX,
            zmax=ncfile.ZMAX,
            source=filename,
            temperature=ncfile.TEMPERATURE,
            hlong=ncfile.HLONG,
            hshort=ncfile.HSHORT,
            thermostat_temp=ncfile.TEMPERATURE,
            thermostat=ncfile.THERMOSTAT,
            mass=ncfile.pmass
        )
    # WORK IN PROGRESS

def read_step(filename,p,step='last'):
    """ Read a step i from a file and make this the state of a particle object. 
        TODO: check that the number of particles is the same!
    """
    try:
        ncfile = netCDF4.Dataset(filename,'r')
    except RuntimeError:
        print 'File not present?',filename
        sys.exit()
    time = ncfile.dimensions['timestep']
    timev = ncfile.variables['timestep']

    r = ncfile.variables['position']
    v = ncfile.variables['velocity']
    #u = ncfile.variables['internal_energy']
    m = ncfile.variables['mass']
    n = p.n

    if step == 'last':
        i = len(time) - 1
    else:
        # test for integer
        i = step
        
    p.r[0:n,:] = r[i,0:n,:]
    p.v[0:n,:] = v[i,0:n,:]
    #p.u[0:n] = u[i,0:n]
    p.m[0:n] = m[i,0:n]
    ncfile.close()






if __name__ == "__main__":
    print 'Testing'
    import particles

    # TEST NC WRITE
    # 1. Create a netcdf file
    attribs = {'name':'Andrew', 'age':33}
    filename = 'test.nc'
    filename2 = 'test2.nc'
    filename3 = 'test3.nc'

    # 2. Create a particle
    SIDE = (3,3,3)
    p = particles.SmoothParticleSystem(27,maxn=30,d=3,rinit='grid',vmax=1.0
        ,side=SIDE,spacing=1.0,xmax=10,ymax=10,zmax=10)
    
    # 3. Write the particle data out
    create_sph_ncfile(filename,attribs,27,3)
    write_step(filename,p)

    # This section is dependent on the VSP modules
    from vsp.spview3d import spview3d

    # Plot the positions
    viewer = spview3d(filename,basedir='.')
    viewer.scatterstep(0)

    # 4. Create another particle structure with a different configuration
    p2 = particles.SmoothParticleSystem(27,maxn=30,d=3,rinit='grid',vmax=1.0
        ,side=SIDE,spacing=0.5,xmax=10,ymax=10,zmax=10)

    # Plot the positions
    create_sph_ncfile(filename2,attribs,27,3)
    write_step(filename2,p2)
    viewer2 = spview3d(filename2,basedir='.')
    viewer2.scatterstep(0)
    
    # 5. Use the previously saved file as the template to inititialise this
    read_step(filename2,p)
    create_sph_ncfile(filename3,attribs,27,3)
    write_step(filename3,p)

    # Plot the positions
    viewer3 = spview3d(filename3,basedir='.')
    viewer3.scatterstep(0)




