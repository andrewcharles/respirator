import sys
sys.path.append('\Users\acharles\masters\active\pyticles\trunk')
sys.path.append('\Users\acharles\masters\active\vasp')
import numpy as np
import pylab as plt
plt.ion()
import fkernel
import sphlib
from eos import vdw_co,vdw_hc
from tensor import symmetric_traceless

""" Test subroutines in sphlib related to velocity gradient 
    and pressure tensor for very small systems of particles.

    Uses the fortran version of code where available.

"""

# Data for a system of two particles
v = np.array([     [1.0,0.0],
                   [0.0,0.0]])
n = 2
ni = 1
d = 2
T = np.zeros([n])
T[:] = 1.0
nlist = np.array( [ [0,1] ],order='F' )
sml = np.array( [ [2.0], [2.0] ] )
w = np.array([[0.0],[0.0],[0.0],[0.0],[0.0]])
dwdx = np.array([[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]])
dv = np.zeros(dwdx.shape,order='F')
rij = np.array([ [1.0], [1.0], [1.41], [3.0], [3.0] ])
drij = np.array([ [1.0,0.0], [0.0,1.0], [1.0,1.0], [1.0,2.0], [2.0,1.0] ])
rho = np.zeros((n))
grad_rho = np.zeros((n,d),order='F')


# Velocity diff
sphlib.sphlib.calc_dv(dv,nlist+1,v)
print 'dv'
print dv

for i in range(ni):
    dv[i,:] = v[nlist[i,0],:] - v[nlist[i,1],:]

print dv

# Kernels and kernel gradients
w,dwdx = fkernel.kernel.smoothing_kernels(rij,drij,nlist,sml,1)
print 'w'
print w
print 'dwdx'
print dwdx

# Density summation
mass = np.ones((n))
fkernel.kernel.density_sum(rho,grad_rho,nlist,sml,mass,w,dwdx,1)
print 'rho'
print rho
print 'grad_rho'
print grad_rho

# Grad v
grad_v = np.zeros((n,d,d),order='F')
sphlib.sphlib.calc_grad_v(grad_v,nlist,dwdx,dv,mass,rho)
print 'grad v'
print grad_v

# Div v
div_v = np.zeros([n],order='F')
sphlib.sphlib.calc_div_v(grad_v,div_v)
print 'div v'
print div_v

#Symmetric traceless velocity gradient
grad_v_os = np.zeros((n,d,d),order='F')
sphlib.sphlib.calc_grad_v_os(grad_v_os,grad_v,div_v)

print 'grad_v_os'
print grad_v_os

aos = symmetric_traceless(grad_v[0,:,:])
print aos[:,:]

#Compute a symmetric traceless part of the pressure tensor
#    -- calc_pi_os(pi_os,grad_v_os,eta,n,ndim)
pi_os = np.zeros((n,d,d),order='F')
eta = np.zeros([n],order='F')
eta[:] = 1
sphlib.sphlib.calc_pi_os(pi_os,grad_v_os,eta)

print 'pi_os'
print pi_os

#Compute the isotropic irreversible pressure
#    -- calc_pi_one(pi_one,div_v,zeta,n)
pi_one = np.zeros([n],order='F')
zeta = np.zeros([n],order='F')
zeta[:] = 1
sphlib.sphlib.calc_pi_one(pi_one,div_v,zeta)

print 'pi_one'
print pi_one

# Isotropic equilibrium (reversible) pressure
# (python implementation)
phc = np.zeros([n],order='F')
pco = np.zeros([n],order='F')

phc = vdw_hc(rho,T)
pco = vdw_co(rho,T)

print 'phc'
print phc

print 'pco'
print pco

P = np.zeros([n,d,d],order='F')

for i in range(n):
    for j in range(d):
        for k in range(d):
            P[i,j,k] = phc[i] + pco[i] +  pi_one[i] + pi_os[i,j,k]

print 'pressure tensor'
np.set_printoptions(precision=5,suppress=True)
print P
