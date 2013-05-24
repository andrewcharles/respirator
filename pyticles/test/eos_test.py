import feos
import numpy as np

thermalk = 1.0
shape = [12]
U = np.zeros(shape)
T = np.zeros(shape)
rho = np.zeros(shape)

U[:] = 0.
T[:] = 1.
rho[:] = 1.

feos.eos.calc_vdw_energy(U,T,rho)
print U
print T

feos.eos.calc_vdw_temp(U,T,rho)
print U
print T
