#!/usr/env python

from sys import argv, stdout

def read_file( fd ):
        lineit = iter(fd)
        E = {}
        f = {}
        if1 = 0
        jf1 = 0
        sym = {}
        n = {}
        s = 0
        for line in lineit:
                if 'CCSD       Excitation energies' in line:
                        for i in range(3):
                            next(lineit)
                        ip = float(next(lineit).split()[7])
                        print(ip)

if __name__ == '__main__':
    with open(argv[1],'r') as fd:
        read_file(fd)
