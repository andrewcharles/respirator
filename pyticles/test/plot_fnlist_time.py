import numpy as np
import matplotlib.pyplot as plt
plt.ion()

#fname = 'pairtime.dat'

for j in range(1):
        fname = 'pairtime' + str(j) + '.dat'
        fname2 = 'fpairtime' + str(j) + '.dat'
        data_in = np.loadtxt(fname)
        data_in2 = np.loadtxt(fname2)
        if j == 0:
            data = data_in
            data2 = data_in2
        else:
            data = data + data_in
            data2 = data2 + data_in2

number = data[0:150,0]
pairs = data[0:150,1]
times = data[0:150,2]
times2 = data2[0:150,2]
#plt.subplot(3,1,1)
#plt.title('Time against pairs')
#plt.plot(pairs,times)
fig = plt.gcf()
quad = fig.add_axes([0.1,0.1,0.75,0.75])
quad.set_ylabel('Execution time')
quad.set_xlabel('Number of particles')
quad.set_ylim([0,0.4])
pair_ax = quad.twinx()
pair_ax.set_ylabel('Number interacting pairs')
p1 = quad.plot(number,times,'b--')
p1 = quad.plot(number,times2,'g--')
p2 = pair_ax.plot(number,pairs,'r-.')
leg = plt.legend((p1,p2),('Execution Time','Interacting Pairs'),loc='left')
plt.title('Time to compute pair separation 3D (2Ghz 64 bit)')

