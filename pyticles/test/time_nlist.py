""" Time the nlist creation and pair seperation.
    Create random particle positions for NP particles
    plot_nlist_time.py will plot the data generated by
    this script.

    Now recording time to build list and compute pair separation.

"""
import particles
import neighbour_list
from time import time
import numpy as np

NP = 1000

def make_list(code,p):
   if code == 'fortran':
        return neighbour_list.FortranVerletList(p,cutoff=10,tolerance=2)
   if code == 'cython':
        return neighbour_list.CythonVerletList(p,cutoff=10,tolerance=2)
   if code == 'fastcython':
        return neighbour_list.FastVerletList(p,cutoff=10,tolerance=2)
   if code == 'python':
        return neighbour_list.VerletList(p,cutoff=10,tolerance=2)
   if code == 'kd':
        return neighbour_list.KDList(p,cutoff=10,tolerance=2)
   if code == 'pdist':
        return neighbour_list.BruteScipy(p)

dt_tot = {}

def time_list(code):
    fname = code + 'pairtime'  + '.dat'
    dt_tot[code] = 0
    ofile = open(fname,'w')

    j = 1
    for n in range(100,NP,100):
        tbuild1 = time()
        jmax = 2 * len(range(100,NP,100))
        p = particles.ParticleSystem3D(n,maxn=n)
        nl = make_list(code,p)
        for i in range(n):
            p.r[i,:] = np.random.random(3) * 10
        nl.build()
        nl.compress()
        tsep1 = time()
        nl.separations()
        # execution time is in seconds
        dt_sep = (time() - tsep1)
        dt_build = (time() - tbuild1)
        ofile.write('%d %d %f %f\n' %(n,nl.nip,dt_sep,dt_build))
        print code,nl.nip,'pairs, %7.6fs, %d of %d' %(dt_sep,j,jmax)
        dt_tot[code] += dt_build
        j += 1
    ofile.close()

time_list('fortran')
time_list('cython')
time_list('fastcython')
time_list('kd')
#time_list('pdist')
#time_list('python')

print 'Comparing neighbour lists.'

for key in dt_tot.keys():
    key,'%7.6f' %(dt_tot[key])


