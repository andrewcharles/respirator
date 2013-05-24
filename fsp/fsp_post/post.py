#! /usr/local/bin/python
"""
    Runs the post processing for an fsph run.
    Takes a command line argument, which is the name of
    the netcdf file to create (without the .nc)
    Compresses and gets rid of sphstate and gofr files.
    The function post_process takes the run name 
    (the filename without the .nc)
        
    Andrew Charles 2008
    
"""
import sys
import os
import file_ops
import glob
# calls the post-processing bits
DIR = os.environ.get('FSP_BASE')

def post_process(name):
    # if there are no sphstate files, don't do the
    # conversion
    if len(glob.glob("sphstate.????????")) > 1:
        os.system(DIR + "/fsp_post/conv_sph_ncdf.py "+ name)
    file_ops.compress_states(".")
    os.system(DIR + "/fsp_post/plotsph.py -f "+ name +".nc " + "--energy --render")
    file_ops.compress_gofrs(".")
    file_ops.trash_pngs(".")

if __name__ == "__main__":
    #test for the existence of the argument!
    post_process(sys.argv[1])
