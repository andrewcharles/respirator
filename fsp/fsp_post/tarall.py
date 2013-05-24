#! /usr/bin/python
""" searches the current directory
 for old style sph run folders
 makes a tar zip of all old state files
 makes a tar zip of all old gofr files
 deletes all pngs
 deletes all old state files
"""

import sys
import os
import shutil
import glob
import file_ops

# vasp modules
#import report

#make a list of all directories
contents = os.listdir('.')
contents = contents + ["."]
#print contents

for a in contents:
    if os.path.isdir(a) and os.path.isfile(a+"/sphvars.var"):
        #for each sph directory
        print a
        os.chdir(a)
       
        file_ops.compress_states(a)
        file_ops.compress_gofrs(a)        
        #file_ops.trash_pngs(a):
       
        os.chdir("..")




