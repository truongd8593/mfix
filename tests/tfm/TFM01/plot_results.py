#!/usr/bin/python

import numpy as np
import pylab as pl

# load data file and ignore top three text lines
with open('solution_x_velocity_profile.dat') as f:
  data_lines = f.readlines()[3:]

data = np.loadtxt(data_lines)

# plot u_g (numerical) and u_g (exact) vs. y variable
pl.figure(1, figsize=(3.0,3.0))
plot1 = pl.plot(data[:,3], data[:,1], 'ro')
plot2 = pl.plot(data[:,4], data[:,1], '-b+', linewidth=2)
pl.xlabel('x-velocity (m/s)')
pl.ylabel('y (m)')
pl.xlim(0.0, 20.0)
pl.ylim(0.0, 0.01)
pl.legend([plot1, plot2], ('MFiX','Exact'), 'best')
pl.savefig('x_velocity_profile.png', bbox_inches='tight')
pl.close(1)

#pl.figure(2, figsize=(4.0,3.0))
#pl.plot(data[:,5], data[:,1], 'ro')
#pl.xlabel('Error (m/s)')
#pl.ylabel('y (m)')
##pl.xlim(0.0, 20.0)
##pl.ylim(0.0, 0.01)
#
#pl.savefig('x_velocity_error.png', bbox_inches='tight')
#pl.close(2)
