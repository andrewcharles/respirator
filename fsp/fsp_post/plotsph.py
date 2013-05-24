#!/usr/bin/env python
"""
:Info: Part of the VSP set of smooth particle plotting modules
:Author: Andrew Charles <ac1201@gmail.com>
:Date: $Date:2008-10-29$

Run as a script, expects a command line argument, which is the name of
the netcdf file to plot. It must be the name, not the
directory. 
Generates plots for sph netCDF files and sphstate.* files

    Command line arguments
    ----------------------
    -f       : Name of sph netcdf file to plot
    --render : Creates a render animation of density
    --trender: Creates a render animation of temperature
    --energy : Creates a line plot of energies
    --state  : Animates the system state on the vdw phase diagram
    --temps  : Plot temperatures (line plot)
    --step   : Plot timestep

"""

import numpy
import scipy
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
import pylab as plt
import matplotlib
import netCDF4
import getopt

import mpegify
import glob
import os
import sys
import sphstate
import renderspam
import file_ops
import shutil
import vsp.vdw.plot_vdw
import time

#Globals
xmin = 0
xmax = 12
ymin = 0
ymax = 12
res = 1.0
filename = "default"
PLOT_ENERGY = False 
RENDER = False 
STATE = False
HIST = False
TRENDER = False
STEP = False
TEMP_SERIES = False

# Process cl arguments
def clargs():
    global filename, PLOT_ENERGY, STATE, HIST, RENDER, TRENDER, STEP, TEMP_SERIES
    try:
        opts,args = getopt.getopt(sys.argv[1:],"f:",["file=","energy",
        "render","state","hist","trender","step=","temps"])
    except getopt.getopterror:
        print "command line argument error"
        sys.exit(2)
    for opt,arg in opts:
        if opt in ("-f","--file"):
            print "plotting " + str(arg)
            filename = str(arg)
        if opt in ("--energy"):
            print "Plotting energies"
            PLOT_ENERGY = True
        if opt in ("--render"):
            print "Rendering"
            RENDER = True
        if opt in ("--state"):
            print "Plotting State"
            STATE = True
        if opt in ("--hist"):
            print "Plotting histograms"
            HIST = True
        if opt in ("--trender"):
            print "Rendering Temperature"
            TRENDER = True
        if opt in ("--temps"):
            print "Temperature time series"
            TEMP_SERIES = True
        if opt in ("--step"):
            print "Plot one step"
            STEP = True
            step = int(arg)

def get_sph_file_list(path):
    """ Creates a list of all sphstate files in a directory. """
    state_files=glob.glob(path+'/sphstate.????????')
    return state_files 

def plot_positions(r,style='ro'):
    plt.plot(r[:,0],r[:,1],style)

def plot_bounds(r):
    plt.plot(r[:,0],r[:,1],'ro',markersize=5.0)

def plot_density_profile(r,rho):
    """ For starters just use the hack's method of plotting
        the particle densities on a line
    """
    plt.plot(r[:,0],rho[:],'bx')

def plot_energies_raw(ncfile):
    """ Plots the internl, kinetic and isolated hamiltonian
        energies as recorded in the sphstate file.
    """
    fig = plt.figure()
    axen = fig.add_axes([0.1,0.1,0.8,0.7])
    k = ncfile.variables["total_kinetic_energy"][:]
    u = ncfile.variables["total_internal_energy"][:]
    l = ncfile.variables["hamiltonian"][:]
    ts = ncfile.variables["thermostat_energy"][:]
    up, = axen.plot(u,label="internal energy")
    kp, = axen.plot(k,label="kinetic energy")
    lp, = axen.plot(l,label="isolated hamiltonian")
    tsp, = axen.plot(ts,label="thermostat_energy")
    #tup, = axen.plot(u+ts,label="total internal energy")
    #tvp, = axen.plot(k+u,label="particle total energy")
    #ptp, = axen.plot(k+u+ts,label="k+u+thermostat")
    #pyp, = axen.plot(k+u-ts,label="k+u-thermostat")
    axen.set_ylabel("Energy")
    axen.set_xlabel("Timesteps")
    axen.set_title("Energies")
    plt.figlegend((up,kp,lp,tsp),(up.get_label(),kp.get_label(),lp.get_label(),tsp.get_label()),(0.7,0.8),markerscale=0.6,prop=FontProperties(size='smaller'))

def plot_energies(u,v,m):
    """ Computes and plots the kinetic and internal energy evolution curves
        takes u,v,m as numpy arrays.
    """
    utot = u.sum(axis=1)
    k = zeros(u.shape)
    for i in range( v.shape[0]):
        for j in range( v.shape[1]):
            k[i,j] = (0.5) * m[i,j]*dot(v[i,j,:],v[i,j,:])
    ktot = k.sum(axis=1)
    internal = plt.plot(utot,label="internal energy")
    kinetic = plt.plot(ktot,label="kinetic energy")
    plt.ylabel("Energy")
    plt.xlabel("Timesteps")
    plt.title("Energies")
    plt.figlegend()

def plot_temperatures(t,n):
    """ Plots the time series of temperatures of normal particles
        and boundary particles.
        n is the number of normal particles.
    """
    for i in range(n):
        plt.plot(t[:,i],'k-')

    for i in range (n,t.shape[1]):
        plt.plot(t[:,i],'r-')
    
    plt.ylabel("Temperature")
    plt.xlabel("Timesteps")
    plt.title("Temperatures")

   
def plot_density_histogram(rho,n):
    """ Plots a density histogram
    """
    nbins = n/10.0 
    plt.hist(rho,bins=nbins,normed=True)
    #xlim(0,2.0)


def plot_parrows(r,a):
    """ Given arrays containing positions and
        acceleration vectors respecitively,
        plots these quantities
    """
    plt.plot(r[:,0],r[:,1],'k.')
    # this next condition is a bit of a hack - for some reason
    # quiver quivers when a is all zero.
    if (not (a.all() == 0.0)):
        plt.quiver(r[:,0],r[:,1],a[:,0],a[:,1],units='dots',width=1, \
                     color='grey')

def renderplot(r,vdot,m,rho,a,h,cutoff,res,bounds):
    """ Uses renderspam to make a nice coloured plot
        bounds: (xmin,xmax,ymin,ymax)
    """
    #global xmax
    #global ymax
    xmin = bounds[0]
    xmax = bounds[1]
    ymin = bounds[2]
    ymax = bounds[3]
    #ax = plt.axes([0.3,0.1,,1],frameon=True)
    plt.axis([xmin,xmax,ymin,ymax])
    #ax.set_axis_off() 
    #ax = plt.axes([xmin,xmax,ymin,ymax],frameon=False)
    gridx,gridy = renderspam.get_grid_map(xmin,xmax,ymin,ymax,res)
    plot_parrows(r,vdot)
    renderspam.renderspam2d(gridx,gridy,r,m,rho,a,h,xmin,xmax,ymin,ymax,cutoff)
    #plt.colorbar(ax=ax)
    #plt.axis([xmin,xmax,ymin,ymax])

def plot_step(ncfile,i,style):
    """ Given a netcdf file, plots the given timestep in the given
        style.
        
        For now this is designed to use interactively with matplotlib
        
        Use a dictionary to map style names to styles.
        ncfile - file object
        i - step to plot
        style - style to plot in

        nb to convert netcdf variable type to  numpy array just take a 
        slice e.g.g v_arr = v[:,:]
    """

    global xmax
    global ymax

    xmax = getattr(ncfile,"xmax")
    ymax = getattr(ncfile,"ymax")
    r = ncfile.variables["position"][i,:,:]
    vdot = ncfile.variables["acceleration"][i,:,:]
    m = ncfile.variables["mass"][i,:]
    a = ncfile.variables["density"][i,:]
    rho = ncfile.variables["density"][i,:]
    h = ncfile.variables["smoothing_length"][i,:]
    cutoff = getattr(ncfile,"cutoff_radius")
    r = numpy.atleast_2d(r)
    plt.plot(r[:,0],r[:,1],'k.')
    #ax = axes([xmin,xmax,ymin,ymax])
    #ax = get_axes()
    #renderplot(r[:,:],vdot[:,:],m[:],rho[:],a[:],h[:],cutoff)
    #title(build_title(ncfile,i))
    #figtext(0.1,0.1,r'$\Gamma$',fontsize=36)
    #ct = ['wtf','wtf','bbq','hax','max','spax']
    #rl = ['r1','r2','r3']
    #cl = ['c1','c2','c3']
    #tbl = plt.table(cellText=[["wtf"],["bbq"]])
    #,rowLabels=rl,colLabels=cl,loc='bottom')
    #ax.add_table(tbl)

def build_title(ncfile,i):
    """ Constructs a title for a spam netcdf plot.
        should check for existence of variables in file.
    """
    title = "" 
    #if 'Run' in ncfile.ncattrs():
    #    title = title  + str(ncfile.Run)
    if 'thermostat_temp' in ncfile.ncattrs():
        title = title + "Temperature "+str(ncfile.thermostat_temp) + " "
    #if 'eos' in ncfile.ncattrs():
    #    title = title + "EOS="+str(ncfile.eos)
    #if 'cgrad' in ncfile.ncattrs():
    #    title = title + "c="+str(ncfile.cgrad) + " "
    
    frame = ncfile.variables["timestep"][i]
    step = frame
    title = title + "frame %04d" %(step) + " "
    if 'dt' in ncfile.ncattrs():
        time = step * ncfile.variables["dt"][0] 
        title = title + "time=" + str(time)

    return title

def plotstate(ncfile):
    """ Plots a nice picture of the state point. 
    """
    global xmax
    global ymax
    m = ncfile.variables["mass"]
    xmax = getattr(ncfile,"xmax")
    ymax = getattr(ncfile,"ymax")
    big_v = xmax*ymax
    big_m = sum(m[0,:]) 
    return plot_vdw.plot_state_point(big_v,big_m,ncfile.thermostat_temp)

def plot_energies_boundary(ncfile):
    """ Plots the boundary particles' energies seperately.
    """
    n = getattr(ncfile,"Number")
    v = ncfile.variables["velocity"][:,:,:]

    if (n < v.shape[1]):
    
        u = ncfile.variables["internal_energy"][:,0:n]
        bu = ncfile.variables["internal_energy"][:,n:]
        #bv = ncfile.variables["velocity"][:,n:,:]
        m = ncfile.variables["mass"]

        utot = u.sum(axis=1)
        butot = bu.sum(axis=1)
        k = zeros(u.shape)
        bk = zeros(bu.shape)
        #for i in range( v.shape[0]):
        #    for j in range(n):
        #        k[i,j] = (0.5) * m[i,j]*dot(v[i,j,:],v[i,j,:])
      
        nb = v.shape[1] - n 
        for i in range(v.shape[0]):
            for j in range(nb):
                bk[i,j] = (0.5) * m[i,n+j]*numpy.dot(v[i,n+j,:],v[i,n+j,:])

        ktot = k.sum(axis=1)
        bktot = bk.sum(axis=1)
        internal = plt.plot(utot,label="internal energy")
        #kinetic = plot(ktot,label="kinetic energy")
        b_internal = plt.plot(butot,label="b internal energy")
        b_kinetic = plt.plot(bktot,label="b kinetic energy")
        plt.ylabel("Energy")
        plt.xlabel("Timesteps")
        plt.title("Energies")
        plt.legend()
        plt.savefig("b_energy.png",format='png')
        plt.clf()


def render_density_all(ncfile):
    """ Plots a set of png renderings of the density for a given
        netcdf file. 
        
        Arguments
        ---------

    """
    #global xmax
    #global ymax
    #global res

    # q = ncfile.variables["heat_flux"]
    t = ncfile.variables["temperature"]
    # x and y are backwards?   

#    res=xmax/50.0
    timesteps = t.shape[0]
    
    #plot_temperatures(t[:,:],n)

    for i in range(timesteps):
        render_density(ncfile,i)
        plt.clf()



def modsim_plot(ncfile,i,indices=[1],imgname='sequence.png'):
    """ Modsim nice plots for one timestep
    """
    res = 1.0
    num_plots = len(indices)
    # Assume 6 frames

    n = ncfile.Number
    cutoff = getattr(ncfile,"cutoff_radius")
    m = ncfile.variables["mass"]

    fig = plt.gcf()
    fig.set_dpi(100)
    fig.set_size_inches(8,12.8)

    # 1 1 3 4 7
    
    xmax = getattr(ncfile,"xmax")
    ymax = getattr(ncfile,"ymax")
    xmin = 0
    ymin = 0
    gridx,gridy = renderspam.get_grid_map(xmin,xmax,ymin,ymax,res)
    plt.clf() 
  
    font = {'family' : 'sans-serif',
            'color'  : 'k',
            'weight' : 'normal',
            'size'   : 10,
            }

    #legend_axes= fig.add_axes([0.1,0.1,1.0,0.2])
    cbar_axes= fig.add_axes([0.88,0.25,0.02,0.7])
    for j in range(num_plots):
        spatial_bounds = [0.3,1.0-(j+1)*0.13,0.6,0.1]
        ax = fig.add_axes(spatial_bounds)
        ax.set_xticks([0,20.0,40.0,60.0])
        if j==5:
            ax.set_xlabel('Density field',font)
            ax.set_xticklabels(['0','20','40','60'],font)
        else:
            ax.set_xticklabels([])
        ax.set_yticks([0,10,20])
        ax.set_yticklabels(['0','','20'],font)
        i = indices[j]
        r = ncfile.variables["position"][i,:,:]
        vdot = ncfile.variables["acceleration"][i,:,:]
        m = ncfile.variables["mass"][i,:]
        a = ncfile.variables["density"][i,:]
        rho = ncfile.variables["density"][i,:]
        h = ncfile.variables["smoothing_length"][i,:]
        img = renderspam.renderspam2d(gridx,gridy,r,m,rho,a,h,xmin,xmax,ymin,ymax,cutoff)
        plot_parrows(r[0:n,:],vdot[0:n,:])
        #plt.axis([xmin,xmax,ymin,ymax])
        title_string = "" 
        title_string = title_string + "T="+str(ncfile.thermostat_temp) + " "
        title_string = title_string + "c="+str(ncfile.cgrad) + " "
        frame = ncfile.variables["timestep"][i]
        step = frame * 50 
        title_string = title_string + "step=" + str(step) + " "
        time = step * ncfile.variables["dt"][0] 
        title_string = title_string + "time=" + str(time)
        ax.text(0.01,0.9,'t=%d' %(time),transform=ax.transAxes,weight='heavy',color=(0.8,0.8,0.8),backgroundcolor=(0.1,0.1,0.1))
    
    plt.colorbar(img,cbar_axes)

    for j in range(num_plots):
        histo_bounds =   [0.1,1.0-(j+1)*0.13,0.2,0.1]
        ax = fig.add_axes(histo_bounds)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([0,0.5,1.0,1.5,2.0])
        lines = plotstate(ncfile)
        ax2 = plt.twinx()
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        
        if j==5:
            ax.set_xticklabels(['0','','1.0','','2.0'],font)
            ax.set_xlabel(r'$\rho$',font)
            ax.set_ylabel('Temperature',font)
        else:
            ax.set_xticklabels([])

        i = indices[j]
        rho = ncfile.variables["density"][i,:]
        plot_density_histogram(rho[:],n)
        plt.xlim(0,2.0)
        ax.set_ylim(0,1.5)
        #plt.text(0,0.95,"Phase diagram")

    #legend_axes.
    lns = lines[0][0],lines[1][0]
    #pt = lines[2][0]
    #print pt
    #print dir(pt)
    #print pt.get_linestyle()
    #print pt.get_marker()
    #leg = plt.figlegend(lns,('Coexistence line','Thermostat temperature')
    #                       ,(0.05,0.12),numpoints=1)
    #leg2 = plt.figlegend(lines[2],('Quenched state'),(0.55,0.05),numpoints=1,markerscale=1)
    #leg.get_lines()[2].set_marker('o')
    #ax.text(0.6,0.1,'State',transform=ax.transAxes)

    #legend_axes.plot([0.6],[0.1],'bo')
    #for t in leg.get_texts():
    #    t.set_fontsize('small')
    plt.draw()
    plt.savefig(imgname,format='png')

def render_density(ncfile,i):
    """ Renders the density at one timestep
    """ 
    # note that by passing slices we are converting from netCDF
    # variables to numpy arrays
    print "Rendering density for step",i
    n = ncfile.Number
    cutoff = getattr(ncfile,"cutoff_radius")
    m = ncfile.variables["mass"]
    r = ncfile.variables["position"][i,:,:]
    vdot = ncfile.variables["acceleration"][i,:,:]
    m = ncfile.variables["mass"][i,:]
    a = ncfile.variables["density"][i,:]
    rho = ncfile.variables["density"][i,:]
    h = ncfile.variables["smoothing_length"][i,:]
    xmax = getattr(ncfile,"xmax")
    ymax = getattr(ncfile,"ymax")
    xmin = 0
    ymin = 0

    if type(cutoff) == type('wtf'):
        cutoff = ymax / 5.0
    renderplot(r[:,:],vdot[:,:],m[:],rho[:],a[:],h[:],cutoff,res
              ,(xmin,xmax,ymin,ymax))
    #plt.title(build_title(ncfile,i))
    #plot_parrows(r[0:n,:],vdot[0:n,:])

    if r.shape[0] > n:
        plot_bounds(r[n:r.shape[0],:])

    num = "%05d" %i
    plt.savefig("ra"+num+"I.png",format='png')

    #plot_density_profile(r[i-1,:,:],rho[i,:])
    #savefig("final_profile.png",format='png')

def render_temperature(ncfile):
    """ Plots a set of png renderings of the density for a given
        netcdf file. Takes the file pointer as its argument
    """
    global xmax
    global ymax
    global res
    r = ncfile.variables["position"]
    q = ncfile.variables["heat_flux"]
    m = ncfile.variables["mass"]
    rho = ncfile.variables["density"]
    h = ncfile.variables["smoothing_length"]
    t = ncfile.variables["temperature"]
    xmax = getattr(ncfile,"xmax")
    ymax = getattr(ncfile,"ymax")
    cutoff = getattr(ncfile,"cutoff_radius")
    timesteps = r.shape[0]

    u = ncfile.variables["internal_energy"]
    m = ncfile.variables["mass"]

    res=1.0
    
    for i in range(timesteps):
        # note that by passing slices we are converting from netCDF
        # variables to numpy arrays
        print "Plotting step",i
        plt.clf()
        plt.title("Temperature Plot")
        renderplot(r[i,:,:],q[i,:,:],m[i,:],rho[i,:],t[i,:],h[i,:],cutoff,res,
            (xmin,xmax,ymin,ymax))
        plot_parrows(r[i,:,:],q[i,:,:])
        num = "%05d" %i
        plt.savefig("temp"+num+"I.png",format='png')
        plt.clf()


def plot_state_hist(ncfile):
    """ 
        Make the state-histogram plot
    """
    global xmax
    global ymax
    global res
    n = ncfile.Number
    rho = ncfile.variables["density"]
    print rho.shape
    timesteps = rho.shape[0]
    res = 1.0
    
    for i in range(timesteps):
        # note that by passing slices we are converting from netCDF
        # variables to numpy arrays
        print "state histogram plot, step",i
        #plt.title(build_title(ncfile,i))
        num = "%05d" %i
        ax1 = plt.subplot(111)
        plotstate(ncfile)
        ax2 = plt.twinx()
        plot_density_histogram(rho[i,:],n)
        plt.xlim(0,2.0)
        plt.savefig("rho"+num+"I.png",format='png')
        plt.clf()

def save_images(ncfile,prefix):
    """ Selects five images to keep and changes their filenames
    """
    r = ncfile.variables["position"]
    range = r.shape[0]
    inc = range/5
    figs = [0,inc,2*inc,3*inc,\
           4*inc,range-1]
    for f in figs:
        num = "%05d" %f
        shutil.copyfile(prefix+num+"I.png",prefix+num+".png.keep")

def rename_kept_images():
    """ this does some renaming
    """
    
    keepers = glob.glob("*.png.keep")
    for keep in keepers:
        newname = keep[:len(keep)-5]
        os.rename(keep,newname)

def getfile(path='default'):
    """ Opens a netcdf file given a path
        Most of the plotting methods expect an open file
        If no path is supplied, defaults to the directory
        name.
    """
    if path == 'default':
        path = os.path.abspath('.')
        name = os.path.basename(path)
        filename = name + ".nc"
        print "Using default filename:",filename
    else:
        filename = path
    if os.path.isfile(filename):
        f = netCDF4.Dataset(filename,"r")
        return f
    else:
        print 'no file'
        return False

def main():
    """ Main function for plotting SPH data.
    """
    global filename
    plt.ioff()
    print "Plotting smooth particle data"
    ta = time.time()
    f = getfile(path = filename)

    if not(f):
        print 'No file'
        return

    if RENDER:
        render_density_all(f)
        mpegify.process(["ra"])
        save_images(f,"ra")
        file_ops.trash_pngs("path")

    if TRENDER:
        render_temperature(f)
        mpegify.process(["temp"])
        save_images(f,"temp")
        file_ops.trash_pngs("path")

    if PLOT_ENERGY:    
        plt.clf()
        plot_energies_raw(f)
        plt.savefig("energy2.png",format='png')

    if STATE: 
        if 'thermostat_temp' in f.ncattrs():
            plt.clf()
            plot_state_hist(f)
            mpegify.process(["rho"])
            save_images(f,"rho")
            file_ops.trash_pngs("path")
            
            #plt.savefig("state.png",format='png')

    if TEMP_SERIES:
        print 'Temperature series'
        t = f.variables["temperature"]
        n = getattr(f,"Number")
        plot_temperatures(t,n)
        plt.savefig("t_av.png",format='png')
        plt.clf()

    if STEP:
        render_density(f,step)
        plt.clf()
        plot_energies_raw(f)
        plt.savefig("energy2.png",format='png')

    f.close()
    
    rename_kept_images()
    #plot_energies_boundary(f)
    print "Finished, took",time.time()-ta,"seconds"

if __name__ == "__main__":
    clargs()
    main()


# DETRITUS
""" def plot_old_style(filename):
    # code fragment for ploting old style (raw) files 
     inflate_sph_ncfile(filename)
     state_files=get_sph_file_list('.')
     for infile in state_files:
         matplotlib
         matplot(infile) """
