from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext# as build_pyx
import numpy as np
from socket import gethostname
import os

if __name__ == '__main__':
    include_dirs = [np.get_include(),'.','splib']
    setup(name = 'splib',
          packages=['splib'],
          package_dir={'splib':'splib'},
          ext_modules=[
             Extension('c_interpolate', ['interpolate.pyx'],include_dirs=include_dirs),
             Extension('c_neighbours2d', ['c_neighbours2d.pyx'],include_dirs=include_dirs),
             Extension('c_neighbours3d', ['c_neighbours3d.pyx'],include_dirs=include_dirs),
             Extension('pairsep', ['pairsep.pyx'],include_dirs=include_dirs),
             Extension('c_forces', ['c_forces.pyx'],include_dirs=include_dirs),
             Extension('c_properties', ['c_properties.pyx'],include_dirs=include_dirs)
             ],
          cmdclass = { 'build_ext': build_ext })

    #os.system('mv c_interpolate.so splib')
    #os.system('mv c_neighbours2d.so splib')
    #os.system('mv c_neighbours3d.so splib')


