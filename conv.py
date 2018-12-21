#!/usr/env python

from sys import argv, stdout
import sys
import numpy as np
from pylab import *

def lorentz(omega, exci, gamma,forza):
    grecopi=3.14159265358979323846264338327950288
    return (gamma*forza)/(grecopi*exci*((omega-exci)**2+gamma**2))

def gaussian(omega, exci, gamma,forza):
    grecopi=3.14159265358979323846264338327950288
    return (forza/(exci*gamma*sqrt(2*grecopi))*exp(-(omega-exci)**2/(2*gamma**2)))*0.02

def Mgaussian(omega, exci, sigma,weight):
    grecopi=3.14159265358979323846264338327950288
    fact=1.0/(sigma*sqrt(2*grecopi))
    func=fact*exp(-(omega-exci)**2/(2*sigma**2))
    return weight*func

def read_file( fd ):

        data = np.genfromtxt(fd)
        exci_auSD1 = data[:,0] #Reads excitation energies in a.u.
        fCCSD1 = (data[:,1]) #Reads oscillator strength

        sigma = 0.8/(2*sqrt(2*(log(2)))) # sigma = 0.3397 eV = 0.0125 a.u. 
        exci_SD1 = exci_auSD1*27.2114

        n=len(exci_auSD1)
        for i in range(0,n,1):
           fCCSD1[i] = fCCSD1[i]

        step = 0.01
        omega = np.arange(400,430,step)
        speCCSD1 = []
        speCCSDM = []

        for i in range(0, len(omega), 1):
#               speCCSD1.append(0.0)
                speCCSDM.append(0.0)
        for i in range(0, len(exci_auSD1), 1):
            for n in range(0, len(omega), 1):
#               speCCSD1[n]=speCCSD1[n]+omega[n]*gaussian(omega[n],exci_auSD1[i],hwhm,fCCSD1[i])
#               speCCSDM[n]=speCCSDM[n]+Mgaussian(omega[n],exci_SD1[i],sigma,fCCSD1[i])*omega[n]
                speCCSDM[n]=speCCSDM[n]+Mgaussian(omega[n],exci_SD1[i],sigma,fCCSD1[i])#*step

        for i in range(0, len(omega), 1):
            print omega[i], speCCSDM[i]


if __name__ == '__main__':
    with open(argv[1],'r') as fd:
        read_file(fd)

