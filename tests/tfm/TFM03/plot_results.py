#!/usr/bin/python

import numpy as np
import pylab as plt

# load data files and ignore top three text lines
with open('u_profile.dat') as fu:
  data_lines = fu.readlines()[3:]

data_u = np.loadtxt(data_lines)

with open('v_profile.dat') as fv:
  data_lines = fv.readlines()[3:]

data_v = np.loadtxt(data_lines)

# plot vertical centerline comparisons in 1st sub plot
# plot horizontal centerline comparisons in 2nd sub plot
plt.figure(1, figsize=(7.0,3.0)) # size of the full figure in inches

plt.subplot(1,2,1)

plot1 = plt.plot(data_u[:,0], data_u[:,1], '-b+', linewidth=2)
plot2 = plt.plot(data_u[:,0], data_u[:,2], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('y (m)')
plt.ylabel('U_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ghia'), 'best')
plt.title('(a)',x=0.5,y=-0.3)

plt.subplot(1,2,2)
plot3 = plt.plot(data_v[:,0], data_v[:,1], '-b+', linewidth=2)
plot4 = plt.plot(data_v[:,0], data_v[:,2], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('V_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ghia'), 'best')
#plt.locator_params(axis='x',nbins=5)
plt.title('(b)',x=0.5,y=-0.3)

plt.subplots_adjust(wspace=0.3)

plt.savefig('tfm03_01.png', bbox_inches='tight')
plt.close(1)
