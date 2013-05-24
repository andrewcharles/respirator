#! /usr/bin/python
import sys
import math
import string
import os
import glob

#this contains utility functions for making sense of my sph data

def find_largest(path):
   
   #count the number of sphstate files
   #print 'searching' + path + '/sphstate.*'
   state_files = glob.glob(path + '/sphstate.*')
   i=0
   number = []
   #for each sphstate file, extract the number (extension) into a list
   for name in state_files:
#      if os.path.basename(name) != 'sphstate.input':
      try:
         number.append( int( os.path.basename(name).strip('sphstate.') ))
         #print number[i], type(number[i])
         i = i+1
      except ValueError:
         pass
   highest = max(number)
   #print highest
   return highest     
      #print os.path.basename(i), type(os.path.basename(i))
      #if type(os.path.basename(state_files[4])) == str:
      #   print 'it is a string' 
#find_largest(sys.argv[1])
