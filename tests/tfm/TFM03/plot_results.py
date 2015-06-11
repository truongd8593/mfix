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

# plot vertical centerline comparisons in 1st plot
# plot horizontal centerline comparisons in 2nd plot
plt.figure(1, figsize=(7.0,3.0))

plt.subplot(1,2,1)

plot1 = plt.plot(data[:,0], data[:,3], 'ro')
plot2 = plt.plot(data[:,0], data[:,4], '-b+', linewidth=2)
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('Temperature (K)')
plt.xlim(0.0, 0.2)
plt.ylim(300.0, 450.0)
plt.legend([plot1, plot2], ('MFIX','Exact'), 'best')
plt.title('(a)',x=0.5,y=-0.3)

plt.subplot(1,2,2)
plot3 = plt.plot(data[:,0], data[:,5], '-b+')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('Error (K)')
#plt.locator_params(axis='x',nbins=5)
plt.title('(b)',x=0.5,y=-0.3)

plt.subplots_adjust(wspace=0.3)

plt.savefig('tfm02_01.png', bbox_inches='tight')
plt.close(1)
