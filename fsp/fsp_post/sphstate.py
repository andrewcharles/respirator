#! /usr/bin/python
""" this loads an sphstate file
    to a numpy array
    and loads the column definition
    to a list
"""

import sys
import os
import string
import numpy
from scipy import *

def load_sphstate(fname):
    """ Loads the data from the sphstate file (ascii) at the path fname 
        read the whole file into an array
    """
    ifile = open(fname)
    data = numpy.loadtxt(ifile)
    ifile.close()
    return data

   
def load_column_mapping(cname):
    """ Under construction!.
        loads the same file as splash to get the column-variable mapping
        Should just define a column mapping dictionary in this file itself.
    """
    columns=[]
    colfile = open(cname)
    cols = colfile.readlines()
    for col in cols:
        columns = columns + [(col.rstrip(),'f4')]
    columns = numpy.dtype(columns)
    colfile.close()
    return columns

def main():
    #need to read sphvars to get n
    #should output it to the netcdf file anyways
    print 'executing main'
    columns = load_column_mapping("columns")
    data = load_sphstate("sphstate.00000500")

if __name__ == "__main__":
    main()


# use the dict() constructor to build dictionaries
# from lists of key-value pairs stored as tuples
# vardic2 =:wdict( [('n',20),('eta',0.1)])
# print vardic2

#old code from trying to make record arrays
#    i=0
#    j=0
#    data.transpose()
#    for d in data:
#        for g in d:
#            data[i]=tuple(data[i])
#            j=j+1
#        i=i+1
#    data = tuple(data)
    # can't get the record array working at the moment
    # data = io.array_import.read_array(ifile)
    # read_array is deprecated
    # dt = array(data,dtype=columns)
    # data = tuple(data)
    # convert to a record array so that we can
    # access the data by named fields    
    # use the column mappings to name each column
    # print columns
    # data = array(data,columns) 
    # recarr = zeros( data.shape, columns )

