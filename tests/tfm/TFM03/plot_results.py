#!/usr/bin/python

import numpy as np
import pylab as plt

# load data files and ignore top three text lines
with open('u_profile.dat') as f:
  data_lines = f.readlines()[3:]

data_u = np.loadtxt(data_lines)

with open('v_profile.dat') as f:
  data_lines = f.readlines()[3:]

data_v = np.loadtxt(data_lines)

with open('ghia_results_u.dat') as f:
  data_lines = f.readlines()[3:]

data_ug = np.loadtxt(data_lines)

with open('ghia_results_v.dat') as f:
  data_lines = f.readlines()[3:]

data_vg = np.loadtxt(data_lines)

# plot vertical centerline comparisons in 1st sub plot
# plot horizontal centerline comparisons in 2nd sub plot
plt.figure(1, figsize=(9.0,4.0)) # size of the full figure in inches

plt.subplot(1,2,1)

plot1 = plt.plot(data_u[:,1], data_u[:,3], '-b', linewidth=2)
plot2 = plt.plot(data_ug[:,1], data_ug[:,2], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('y (m)')
plt.ylabel('U_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'best')
plt.title('(a) x=0.5, Re=100',x=0.5,y=-0.25)

plt.subplot(1,2,2)
plot3 = plt.plot(data_v[:,0], data_v[:,3], '-b', linewidth=2)
plot4 = plt.plot(data_vg[:,0], data_vg[:,2], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('V_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'best')
#plt.locator_params(axis='x',nbins=5)
plt.title('(b) y=0.5, Re=100',x=0.5,y=-0.25)

plt.subplots_adjust(wspace=0.3)

plt.savefig('tfm03_01.png', bbox_inches='tight')
plt.close(1)
