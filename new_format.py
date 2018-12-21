#!/usr/bin/env python

import sys
import numpy as np

data = np.loadtxt(sys.argv[1], skiprows=2)
with open(sys.argv[1], 'r') as f:
    unit = f.readline().strip()

with open(f'new_{sys.argv[1]}', 'w') as f:
    f.write(f'!Converted with {sys.argv[0]} from {sys.argv[1]}\n')
    f.write('@COORDINATES\n')
    f.write(f'{data.shape[0]}\n')
    f.write(f'{unit}\n')
    for i in range(data.shape[0]):
        f.write(f'X     {data[i,1]:12.9f} {data[i,2]:12.9f} {data[i,3]:12.9f}    {i+1}\n')
    f.write('@MULTIPOLES\n')
    f.write('ORDER 0\n')
    f.write(f'{data.shape[0]}\n')
    for i in range(data.shape[0]):
        f.write(f'{i+1:<6d} {data[i,4]:12.9f}\n')
    f.write('@POLARIZABILITIES\n')
    f.write('ORDER 1 1\n')
    f.write(f'{data.shape[0]}\n')
    for i in range(data.shape[0]):
                             #xx               xy           xz          yy               yz         zz
        f.write(f'{i+1:<6d} {data[i,5]:12.9f} {0.0:12.9f} {0.0:12.9f} {data[i,5]:12.9f} {0.0:12.9f} {data[i,5]:12.9f}\n')
    f.write('EXCLISTS\n')
    f.write(f'{data.shape[0]} 3\n')
    for i in range(data.shape[0]//3):
        f.write(f'{i*3+1} {i*3+2} {i*3+3}\n')
        f.write(f'{i*3+2} {i*3+1} {i*3+3}\n')
        f.write(f'{i*3+3} {i*3+1} {i*3+2}\n')

    
