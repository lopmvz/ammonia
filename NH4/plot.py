#!/usr/env python
from sys import argv, stdout
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

filelist = glob.glob('PE_NH4*.dat')
data = np.genfromtxt('PE_NH4_4H2O_step2220.dat')
YT = data[:,1]
YT *= 0
n = 0

for fname in filelist:
#   print fname
#   print n
    n+=1
    data = np.genfromtxt(fname)
    X = data[:,0]
    Y = data[:,1]
    plt.plot(X,Y,'cadetblue', linewidth=1)
    for i in range(0,len(Y),1):
        YT[i] += Y [i]


ip_T = 0
n_ip = 0

filelist_ip = glob.glob('*IP.dat')
for fname_ip in filelist_ip:
    n_ip += 1
    ip = np.genfromtxt(fname_ip)
    ip_T += ip
    plt.plot([ip,ip],[0,1],'skyblue',linewidth=1.0, dashes=[5, 2])


YT /= n
ip_T /= n_ip

#for i in range(0, len(X), 1):
#    print X[i], YT[i]
plt.plot(X,YT,'k-', linewidth=2)
plt.plot([ip_T,ip_T],[0,1],'k-', linewidth=2)
plt.ylim(0, 0.14)
plt.xlim(400, 420)
plt.title('$NH_{4}$ + $4H_{2}O}$ (6-311G**)')
plt.xlabel('$\omega$ (eV)')
plt.ylabel('Intensity (arb. units)')
plt.legend(loc='best',fontsize='small')
plt.savefig('PE_NH4_4H2O_6311Gss.pdf')
plt.show()
