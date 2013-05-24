""" Tensor maths functions.
"""

import numpy as np

def symmetric_traceless(a):
    # should test for squareness
    aos = a.copy()
    sz = a.shape
    ndim = sz[0]
    for i in range(ndim):
        for j in range(ndim):
                aos[i,j] = (1./2) * (a[i,j] + a[j,i])
    for i in range(ndim):
        aos[i,i] = aos[i,i] - (1./ndim)*np.trace(a)
    return aos
