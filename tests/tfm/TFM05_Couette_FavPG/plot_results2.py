#!/usr/bin/python

import numpy as np
import pylab as plt

# load data file and ignore top three text lines
with open('N128_Uprofile.dat') as f:
  data_lines = f.readlines()[3:]
data1 = np.loadtxt(data_lines)

with open('N256_Uprofile.dat') as f:
  data_lines= f.readlines()[3:]
data2 = np.loadtxt(data_lines)

with open('N512_Uprofile.dat') as f:
  data_lines = f.readlines()[3:]
data3 = np.loadtxt(data_lines)

#with open('N1024_Uprofile.dat') as f:
with open('../test2/N1024_Uprofile.dat') as f:
  data_lines= f.readlines()[3:]
data4 = np.loadtxt(data_lines)

#with open('N8_denorm_Superbee_U10.dat') as f:
#  data_lines = f.readlines()[4:]
#error1 = np.loadtxt(data_lines)

#with open('N16_denorm_Superbee_U10.dat') as f:
#  data_lines= f.readlines()[4:]
#error2 = np.loadtxt(data_lines)

#with open('N32_denorm_Superbee_U10.dat') as f:
#  data_lines = f.readlines()[4:]
#error3 = np.loadtxt(data_lines)

#with open('N64_denorm_Superbee_U10.dat') as f:
#  data_lines= f.readlines()[4:]
#error4 = np.loadtxt(data_lines)


# plot numerical and exact solution in first subplot
# plot error in second subplot
plt.figure(1, figsize=(7.0,3.0))

plt.subplot(1,2,1)

#exacty = np.arange(0.0,0.101,0.0125)
#exactU = exacty*(10.0/0.1)

plot1 = plt.plot(data1[:,3], data1[:,1], '-k',linewidth=2)
plot2 = plt.plot(data2[:,3], data2[:,1], '-b', linewidth=2)
plot3 = plt.plot(data3[:,3], data3[:,1], '-r',linewidth=2)
plot4 = plt.plot(data4[:,3], data4[:,1], '-r',linewidth=2)
plot5 = plt.plot(data1[1::8,4], data1[1::8,1], '-ks')
plt.xlabel('x-velocity (m/s)',labelpad=10)
plt.ylabel('y (m)')
plt.xlim(0.0, 15.0)
plt.ylim(0.0, 0.01)
plt.legend([plot1, plot2, plot3, plot5], ('N=128','N=256','N=512','Exact'),loc=2,prop={'size':8})
#plt.title('(a)',x=0.5,y=-0.4)

plt.subplot(1,2,2)

#N = np.array([0.0125,0.00625,0.003125,0.0015625])
#error = np.array([error1[1,1],error2[1,1],error3[1,1],error4[1,1]])
#order1x = np.arange(0.001,0.05,0.002)
#order1y = np.arange(0.001,0.05,0.002)
#order1y = pow(order1x,1)
#order2x = np.arange(0.001,0.05,0.002)
#order2y = pow(order2x,2)

#plot1 = plt.loglog(N,error,'-k',linewidth=2)
#plot2 = plt.loglog(order1x,order1y,'-b',linewidth=2)
#plot3 = plt.loglog(order2x,order2y,'-r',linewidth=2)
#plt.xlabel('N',labelpad=20)
#plt.ylabel('L_2 norm')
#plt.legend([plot1, plot2, plot3],('MFIX','O(1)','O(2)'),loc=4, prop={'size':10})

# plt.xlim(0.0, 10.0)
#plt.ylim(0.004, 0.006)

plot1 = plt.plot(abs(data3[:,5]/data3[:,4]),data3[:,1],'-k',linewidth=2)
#plt.xlim(-1e-2, 0.0*1e-2)
plt.ylim(0.0, 0.01)
plt.xlabel('Relative error',labelpad=10)
plt.ylabel('y (m)')
plt.ticklabel_format(style='sci',axis='x',scilimits=(-1,-4))

plt.subplots_adjust(wspace=0.5)

plt.savefig('FavPG_U10.eps', bbox_inches='tight')
plt.savefig('FavPG_U10.png', bbox_inches='tight')
plt.close(1)

