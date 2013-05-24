#! /usr/bin/python

# various operations for handling sph files
# mostly to deal with old formats and handle
# compress, copying, etc
# This is very specific to fsp

import sys
import glob
import os

def compress_states(a):
    """compress and delete all sphstate files in path a"""
    state_files=glob.glob("sphstate.????????")
    state_files.append("properties.output")
    if len(state_files) > 0:
        os.system("tar -cvzf sphstate_all.tgz sphstate.*")
        if(os.path.isfile("properties.output")):
            os.system("tar -cvzf properties.tgz properties.output")
        print "removing sphstate files"
    for st in state_files:
        try:
            os.remove(st)
        except os.error:
            pass

def compress_gofrs(a):
    """compress and delete all gofr files in path a"""
    gofrs=glob.glob("gofr.????????")
    if len(gofrs) > 0:
        os.system("tar -cvzf gofr_all.tgz gofr.*")
        print "removing gofr files"
    for gr in gofrs:
         try:
             os.remove(gr)
         except os.error:
             pass

def trash_pngs(a):
    """ The trailing I before the png is a convention that means
        the image is only temporary and should be deleted.
    """
    pngs=glob.glob("*I.png")
    if len(pngs) > 0:
        print "removing pngs files"
    for pn in pngs:
         try:
             os.remove(pn)
         except os.error:
             pass

    
