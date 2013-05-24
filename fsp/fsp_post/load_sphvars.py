#! /usr/local/bin/python
"""  This loads the sphvars.var file into
     a python dictionary
     A python dictionary is an unordered set of
     key-value pairs.
    
     If called from the command line it spits out a bit of tex
     describing the config file. (not yet implemented)

"""


import sys
import os
import string
import generate_html
import pickle

#vardic = {'timestep': 0.001, 'temperature':1.05}
#vardic['size'] = 3
#print vardic['timestep']
#print vardic

#NB note that python will display the binary rounded
#value (e.g. 0.7 is actually 0.6999999000011 or so
# There is a decimal module that might be worth looking at
# to resolve this.

#Step 1

#load the sphvars file
#create a dictionary entry with the appropriate label

vardic = {}

def main():
    load_sphvars()
    outputstuff =  'number of particles:%d <br />' %(vardic['n'])
    outputstuff += 'dimension:%d <br />' %(vardic['dim'])
    outputstuff += 'timestep:%d <br />' %(vardic['timestep'])
    outputstuff += 'x length:%d <br />' %(vardic['x_len'])
    outputstuff += 'y length:%d <br />' %(vardic['y_len'])
    outputstuff += 'starting temperature:%d <br />' %(vardic['start_temp'])

    generate_html.generate_html(outputstuff)

def load_mapping():
    """ loads a mapping of an fsph input file
        to python variable dictionary names """
    ifile = open('~/masters/active/vasp/sphvars_mapping.txt','r' )


def load_sphvars_new(ifile):
    """ Loads sph initialisation variables from the new style
        of key-value configuration file.
    """
    for line in ifile:
        lsp = line.split()
        if len(lsp) > 1:
            if (lsp[0][0] != '#') :
                vardic[lsp[0]] = lsp[1]
    return vardic


def write_sphvars(vars,oname):
    """ Takes an sphvars dictionary and writes the sph config
        file out. It does not preserve the order. Oh well.
    """
    print vars
    ofile = open(oname,'w')
    for key,val in vars.items():
       ofile.write(key + " " + str(val) + "\n") 
    
 

def load_sphvars(ofile):
    """ Loads old line-dependent config file 
    """
#   ofile = open(parent+'/sphvars.var','r' )
    fclist = ofile.readlines()
    j=0
    for i in fclist:
        if i == "Input_data\n":
            # 1  number of particles (must be square)900,2485
            line = fclist[j+1]
            linelist = line.split()
            vardic['n'] = int(linelist[0])
            #print linelist, vardic['n']

            # 2                   dimension of system
            line = fclist[j+2]
            linelist = line.split()
            vardic['dim'] = int(linelist[0])
            #print linelist

            # 3     max interacting pairs404550,3123750
            line = fclist[j+3]
            linelist = line.split()
            vardic['max_inter'] = int(linelist[0])
            #print linelist

            # 4      number of steps to take
            line = fclist[j+4]
            linelist = line.split()
            vardic['max_steps'] = int(linelist[0])
            #print linelist, vardic['max_steps']

            # 5       timestep size
            line = fclist[j+5]
            linelist = line.split()
            vardic['timestep'] = float(linelist[0])
            #print linelist

            # 6   snapshot frequency
            line = fclist[j+6]
            linelist = line.split()
            vardic['snap_freq'] = int(linelist[0])
            #print linelist

            # 7       cutoff radius
            line = fclist[j+7]
            linelist = line.split()
            vardic['cutoff'] = int(linelist[0])
            #print linelist

            # 8.     debug on/off
            line = fclist[j+8]
            linelist = line.split()
            vardic['debug'] = bool(linelist[0])
            #print linelist

            # 9      debug level 2 on/off (verbose)
            line = fclist[j+9]
            linelist = line.split()
            vardic['verbose'] = bool(linelist[0])
            #print linelist

            # 10     confirm step?
            line = fclist[j+10]
            linelist = line.split()
            vardic['confirm'] = bool(linelist[0])
            #print linelist

            # 11   x dimension of box
            line = fclist[j+11]
            linelist = line.split()
            vardic['x_len'] = int(linelist[0])
            #print linelist

            # 12       y dimension of box
            line = fclist[j+12]
            linelist = line.split()
            vardic['y_len'] = int(linelist[0])

            # 13     side of starting box10.8,21.6
            line = fclist[j+13]
            linelist = line.split()
            vardic['box_side'] = float(linelist[0])

            # 14      spacing of particles
            line = fclist[j+14]
            linelist = line.split()
            vardic['spacing'] = float(linelist[0])

            # 15     starting temperature
            line = fclist[j+15]
            linelist = line.split()
            vardic['start_temp'] = float(linelist[0])

            # 16     Viscous entropy type (1= Liu, 2=Sigalotti ("physical viscosity")
            line = fclist[j+16]
            linelist = line.split()
            vardic['visc_type'] = int(linelist[0])

            # 17     Integration step type (1=Improved Euler, 2=Leapfrog (sig), 3=runge_kutta )
            line = fclist[j+17]
            linelist = line.split()
            vardic['integrator'] = int(linelist[0])

            # 18    Run type (0=normal, 1=load from file, 2=load from file but set temp)
            line = fclist[j+18]
            linelist = line.split()
            vardic['run_type'] = int(linelist[0])

            # 19    Artificial viscosity (0=none, 1=vnr, 2=monaghan beta term)
            line = fclist[j+19]
            linelist = line.split()
            vardic['art_visc'] = int(linelist[0])

            # 20    Conduction type (1= monaghan, 2=full heat flux)
            line = fclist[j+20]
            linelist = line.split()
            vardic['q_type'] = int(linelist[0])

            # 21    Profile (timing) on/off
            line = fclist[j+21]
            linelist = line.split()
            vardic['profile'] = bool(linelist[0])

            # 22    Sigma (repulsive core size)
            line = fclist[j+22]
            linelist = line.split()
            vardic['sigma'] = float(linelist[0])
            #print linelist, vardic['sigma']

            # 23    Repulsive coefficient (negative for repulsion)
            line = fclist[j+23]
            linelist = line.split()
            vardic['repulsive_coeff'] = float(linelist[0])

            # 24    Whether to restrict maximum velocity
            line = fclist[j+24]
            linelist = line.split()
            vardic['restrict_v'] = int(linelist[0])

            # 25    Density equation (0=summation, 1=continuity)
            line = fclist[j+25]
            linelist = line.split()
            vardic['density_type'] = int(linelist[0])

        #if i == "Particle_input_data\n":
        
        j=j+1
        # it would be neat to get the mapping from another file somewhere
        # create a dictionary (just for the first set for now)
#   ofile.close()
#   print vardic
    return vardic


def getPop(vf):
    """get the number of particles, for old style sphvars files"""
    #vf is a pointer to a file os sphvars.var type
    vf.seek(0)
    fclist = vf.readlines()
    j=0
    for i in fclist:
        if i == "Input_data\n":
            line = fclist[j+1]
            linelist = line.split()
            return int(linelist[0])
        j=j+1
    print "should not hit this point"


   
def getTimestep(vf):
    """extract the timestep from the sphvars.var file"""
    #vf is a pointer to a file os sphvars.var type
    vf.seek(0)
    fclist = vf.readlines()
    j=0
    for i in fclist:
        if i == "Input_data\n":
            line = fclist[j+5]
            linelist = line.split()
            return float(linelist[0])
        j=j+1
    print "should not hit this point"
    return none

#main()


# use the dict() constructor to build dictionaries
# from lists of key-value pairs stored as tuples
# vardic2 = dict( [('n',20),('eta',0.1)])
# print vardic2
