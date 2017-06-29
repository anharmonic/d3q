from __future__ import print_function
from numpy import loadtxt
import sys

M = loadtxt(sys.argv[1])

c1 = int(sys.argv[2])-1
c2 = int(sys.argv[3])-1
c3 = int(sys.argv[4])-1
c4 = int(sys.argv[5])
f=float(sys.argv[6])

print('set object {0} polygon '.format(c4), end='')
for i in range(0,len(M)):
    if (i == 0):
        print('from {0},{1} '.format(M[i][c1], M[i][c2]+f*M[i][c3]), end='')
    else:
        print('to {0},{1} '.format(M[i][c1], M[i][c2]+f*M[i][c3]), end='')
for i in range(len(M)-1,-1,-1):
    print('to {0},{1} '.format(M[i][c1], M[i][c2]-f*M[i][c3]), end='')

