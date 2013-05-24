import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import platform
plt.ion()

#lists = ['fortran','cython','fastcython','kd','pdist']
#colors = {  'fortran':'b',
#            'cython':'r',
#            'fastcython':'g',
#            'kd':'k',
#            'pdist':'c'
#         }

lists = ['fastcython','kd']
colors = {  
            'fastcython':'g',
            'kd':'k',
         }

data = {}
number = {}
pairs = {}
times = {}
nline = {}
pline = {}

for list_type in lists:
    fname = list_type + 'pairtime' + '.dat'
    ofile = open(fname,'r')
    data = np.loadtxt(ofile)
    n = data.shape[0]
    number[list_type] = data[0:n,0]
    pairs[list_type] = data[0:n,1]
    times[list_type] = data[0:n,2]

regs = {}
for list_type in lists:
    regs[list_type] =  linregress(pairs[list_type],times[list_type])
    print list_type,regs[list_type][0]

#print regs['fastcython'][0]/regs['cython'][0] 

fig = plt.gcf()
quad = fig.add_axes([0.14,0.1,0.75,0.75])
quad.set_ylabel('Number of pairs')
quad.set_xlabel('Number of particles')

for list_type in lists:
    pline[list_type], = quad.plot(number[list_type],pairs[list_type],
        colors[list_type]+'o')

leg = plt.legend(
    [pline[key] for key in sorted(pline.keys())],
    sorted(pline.keys()),
    loc='upper left')
plt.title('Pairs vs N ' + platform.platform())

plt.savefig('kd_vs_verlet.png')

