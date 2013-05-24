#!/usr/bin/env python

""" 
Generates a restructured text report in a spac output directory.

rundir - location of netcdf output file
reportdir - path of output directory

"""

# This script moves through all subdirectories
# of the specified directory, and looks for
# sphvars.var files. It returns a report with the name
# contents.txt, containing selected information from
# sphvars.var

# Externals
import sys
import os
import shutil
import glob

# Built externals
import netCDF4

# This package
import sph_results
import load_sphvars

def make_report(rundir,reportdir='None'):
    if not os.path.isdir(rundir):
        print "No run directory."
        return

    if reportdir == 'None':
        reportdir = rundir

    if not os.path.isdir(reportdir):
        print "Creating report directory"
        os.mkdir(reportdir)

    contents = os.listdir(rundir)
    ofile = open(reportdir + '/report.txt','w' )
    ofile.write('Run Name:' + rundir + '\n\n')

    # List the netcdf files in this directory
    nc_files = glob.glob(rundir + '/*.nc')
    print nc_files
        
    # Open the netcdf file
    f = netCDF4.Dataset(nc_files[0],"r")

    # Write out the important parts of the configuration
    for att in f.ncattrs():
        ofile.write('| ' + str(att) + ' ' + str(getattr(f,str(att))) + '\n')

    ofile.write('\n')
    # Assume plotsph has been run

    # List the images in this directory
    state_images = glob.glob(rundir + '/ra*.png')

    for img in state_images:
        ofile.write('.. image:: ' + os.path.basename(img)+'\n')
        ofile.write('    :width: 120px\n\n')
    
    # We need an optional report location
    # to which everything will be copied if it is not the
    # same location as the run directory

    if not(rundir == reportdir):
        print 'Copying image files...'
        for img in state_images:
            shutil.copy(img,reportdir)
        shutil.copy(rundir + '/energy2.png',reportdir)

    ofile.write('.. image:: energy2.png\n')
    ofile.write('    :width: 180px\n\n')

    ofile.close()
    os.system("rst2latex.py --documentclass=scrartcl --documentoptions=9pt,a4paper,twocolumn --stylesheet=/Users/acharles/local/rstlatex.sty "
    + reportdir + "/report.txt " + reportdir + "/report.tex")
    
    os.system("pdflatex report.tex")


def main():
    """main subroutine contains high level flow"""
        
    print 'Looking for smooth particle results in' + sys.argv[1]
    make_report(sys.argv[1])

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Expecting more arguments"
        #raise Exception
        sys.exit()
    main()
