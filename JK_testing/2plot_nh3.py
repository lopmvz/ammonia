import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

data0 = np.genfromtxt('PE_NH3_4H2O_step5000_631G_crit2.out.txt')
CCSD0_eV = data0[:,0] * 27.2114
fCCSD0 = (data0[:,1])

data0 = np.genfromtxt('PE_NH3_4H2O_step5000_631G_crit2.out.dat')
omega = data0[:,0] 
speCCSD0 = (data0[:,1])

data1 = np.genfromtxt('PE_nto_NH3_4H2O_step5000_nd2.out.txt')
CCSD1_eV = data1[:,0] * 27.2114
fCCSD1 = (data1[:,1])

data1 = np.genfromtxt('PE_nto_NH3_4H2O_step5000_nd2.out.dat')
omega = data1[:,0] 
speCCSD1 = (data1[:,1])

data11 = np.genfromtxt('PE_NH3_4H2O_step5000_631gss.out.txt')
CCSD11_eV = data11[:,0] * 27.2114
fCCSD11 = (data11[:,1])

data11 = np.genfromtxt('PE_NH3_4H2O_step5000_631gss.out.dat')
omega = data11[:,0] 
speCCSD11 = (data11[:,1])

data2 = np.genfromtxt('6-311Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.txt')
CCSD2_eV = data2[:,0] * 27.2114
fCCSD2 = (data2[:,1])

data2 = np.genfromtxt('6-311Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.dat')
omega = data2[:,0] 
speCCSD2 = (data2[:,1])

data3 = np.genfromtxt('6-311Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.txt')
CCSD3_eV = data3[:,0] * 27.3114
fCCSD3 = (data3[:,1])

data3 = np.genfromtxt('6-311Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.dat')
omega = data3[:,0] 
speCCSD3 = (data3[:,1])

data4 = np.genfromtxt('6-311++Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.txt')
CCSD4_eV = data4[:,0] * 27.2114
fCCSD4 = (data4[:,1])

data4 = np.genfromtxt('6-311++Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.dat')
omega = data4[:,0] 
speCCSD4 = (data4[:,1])

data5 = np.genfromtxt('6-311++Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.txt')
CCSD5_eV = data5[:,0] * 27.3114
fCCSD5 = (data5[:,1])

data5 = np.genfromtxt('6-311++Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.dat')
omega = data5[:,0] 
speCCSD5 = (data5[:,1])

data6 = np.genfromtxt('6-311++Gss_6-311Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.txt')
CCSD6_eV = data6[:,0] * 27.2114
fCCSD6 = (data6[:,1])

data6 = np.genfromtxt('6-311++Gss_6-311Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.dat')
omega = data6[:,0] 
speCCSD6 = (data6[:,1])

data7 = np.genfromtxt('6-311++Gss_6-311Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.txt')
CCSD7_eV = data7[:,0] * 27.3114
fCCSD7 = (data7[:,1])

data7 = np.genfromtxt('6-311++Gss_6-31Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.dat')
omega = data7[:,0] 
speCCSD7 = (data7[:,1])

data8 = np.genfromtxt('6-311++Gss_6-311Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.txt')
CCSD8_eV = data8[:,0] * 27.211
fCCSD8 = (data8[:,1])

data8 = np.genfromtxt('6-311++Gss_6-31Gss/PE_NH3_4H2O_step5000_nd_ecp_new_NH3_4H2O_step5000.out.dat')
omega = data8[:,0] 
speCCSD8 = (data8[:,1])

data9 = np.genfromtxt('6-311++Gss_6-31Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.txt')
CCSD9_eV = data9[:,0] * 27.3114
fCCSD9 = (data9[:,1])

data9 = np.genfromtxt('6-311++Gss_6-31Gss/PE_NH3_4H2O_step5000_nd_new_NH3_4H2O_step5000.out.dat')
omega = data9[:,0] 
speCCSD9 = (data9[:,1])

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
#subplot(3,1,1)
plt.suptitle('$NH_{3}$ + $4H_{2}O (l)$')
#xticks([]), yticks([])

plt.xlabel('Excitation energy (eV)')
plt.ylabel('                                        Intensity (arb. units)')

plt.xlim(400, 420)
plt.ylim(0, 0.075)
#ax1.plot(omega, speCCSD1,'k-', label='6-311G** (Marta)', linewidth=1)
#ax2.plot(omega, speCCSD2,'r-', label='6-311G** new ecp', linewidth=1)
#ax3.plot(omega, speCCSD3,'g-', label='6-311G** new', linewidth=1)
#ax4.plot(omega, speCCSD4,'y-', label='6-311++G** new ecp', linewidth=1)
#ax5.plot(omega, speCCSD5,'m-', label='6-311++G** new ', linewidth=1)
ax1.plot(omega, speCCSD1, '#000000', label='6-311G** (Marta)', linewidth=1)
ax1.plot(omega, speCCSD0, '#2f4f4f', label='6-311++G**/6-31G** (Marta)', linewidth=1)
ax1.plot(omega, speCCSD11,'#696969' , label='6-31G** (Marta)', linewidth=1)

ax2.plot(omega, speCCSD2, '#000080', label='6-311G** new ecp', linewidth=1)
ax2.plot(omega, speCCSD3, '#6495ed', label='6-311G** new', linewidth=1)

ax2.plot(omega, speCCSD6, '#b22222', label='6-311++G**/6-31G** new ecp', linewidth=1)
ax2.plot(omega, speCCSD7, '#f08080', label='6-311++G**/6-31G** new ', linewidth=1)

ax2.plot(omega, speCCSD4, '#228b22', label='6-311++G** new ecp', linewidth=1)
ax2.plot(omega, speCCSD5, '#9acd32', label='6-311++G** new ', linewidth=1)

ax2.plot(omega, speCCSD8, '#ff8c00', label='6-311++G**/6-311G** new ecp', linewidth=1)
ax2.plot(omega, speCCSD9, '#f4a460', label='6-311++G**/6-311G** new ', linewidth=1)

ax1.legend(loc='best',fontsize='small')
ax2.legend(loc='best',fontsize='small')
#n=len(CCSD1_eV)
#for i in range(n):
#    plt.plot([CCSD1_eV[i],CCSD1_eV[i]],[0,fCCSD1[i]],'r-',linewidth=1.0)

#subplot(3,1,2)
#xticks([]), yticks([])
#plt.xlim(400, 420)
#plt.ylim(0, 0.05)
#plt.plot(omega, speCCSD2, 'b-', label='6-311G**', linewidth=1)
plt.legend(loc='best',fontsize='small')
#n=len(CCSD2_eV)
#for i in range(n):
#    plt.plot([CCSD2_eV[i],CCSD2_eV[i]],[0,fCCSD2[i]],'b-',linewidth=1.0)
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

#plt.savefig('comp_ryd_gas_nh3_step5000.pdf')
plt.savefig('comp_nh3_step5000_2.pdf')
plt.show()
