#!/usr/bin/env python

print """

This is a very rough script that imports the FORCE_CONSTANT
file produced by phonopy (by A. Togo). It assumes that the
file has been produced by a Quantum-ESPRESSO calculation, 
(i.e. units Ry/bohr^2)

This script will try to read two files:
1. "pw.in" input file for the unit cell calculation
2. "FORCE_CONSTANTS" with the force costants from phonopy

In order to generate file 2, you need to run phonopy with a
config file containing the option
  FORCE_CONSTANTS = WRITE
I.e. if you put it in a file called mkcf.conf, you can then run:
phonopy --pwscf --dim="2 2 2" -c pw.in -p mkfc.conf

This code uses the ASE python library, you can usually install
it with:
 pip install ase
or by hand (https://wiki.fysik.dtu.dk/ase/).

I do not have time to add support for any other code,
you should be able to do it yourself. Please note that the
units of measure used inside FORCE_CONSTANTS depend
on the code used to compute them.

The resulting FC file will be saved as 'fc'.

NOTE THAT THE fc FILE IS NOT SUITABLE FOR USE IN THERMAL2
BECAUSE IT IS NOT OPTIMIZED FOR FOURIER INTERPOLATION
YOU WILL HAVE TO CONVERT IT TO THERMAL2 USING THE SCRIPT
fc2mat2R.sh FROM THE TOOLS DIRECTORY!

"""

MASS_DALTON_TO_RY = 0.5*1822.88839

HARTREE_SI        = 4.35974394E-18
ELECTRONVOLT_SI   = 1.602176487E-19
AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
RYTOEV           = AUTOEV / 2.
AUTOEV            = HARTREE_SI / ELECTRONVOLT_SI

BOHR_RADIUS_ANGS  = .52917721067
ANGSTROM_AU       = 1./BOHR_RADIUS_ANGS

# Force constants seem to be already in Ry/bohr^2, at least when
# produced from QE calculation
F_FACTOR = 1. #/RYTOEV/ANGSTROM_AU**2

print "Force constant unit conversion factor:", F_FACTOR

def atoms2species(listat):
  nat = len(listat)
  idxat = [0]*nat
  listsp = []
  invsp = []
  #
  for i in range(nat):
    try:
      idxat[i] = listsp.index(listat[i])
    except ValueError:
      listsp.append(listat[i])
      idxat[i] = len(listsp)-1
      invsp.append(i)
  return listsp, invsp, idxat

def supercell2cell(nat, atuc, natsc, atsc, cell, rcell, nq1,nq2,nq3):
  from numpy import dot, transpose, zeros, int
  listrbig = zeros(shape=(natsc,3),dtype=int)
  idxsc = [0]*natsc
  trcell = transpose(rcell)
  for i in range(natsc):
    for j in range(nat):
      rbig = atsc[i,:] - atuc[j,:]
      rbigx =  dot(trcell, rbig)
      rbigi = [round(rbigx[0]), round(rbigx[1]), round(rbigx[2])]
      delta = norm(rbigx-rbigi)
      if delta < 1.e-6:
        listrbig[i][0] = int(rbigi[0]) % nq1
        listrbig[i][1] = int(rbigi[1]) % nq2
        listrbig[i][2] = int(rbigi[2]) % nq3
        idxsc[i] = j
        break
  return idxsc, listrbig

from sys import stderr
import ase.io.espresso
fin = open("pw.in")
aseuc = ase.io.espresso.read_espresso_in(fin)
fin.close()

nat  = len(aseuc)
csym = aseuc.get_chemical_symbols()
cmass= MASS_DALTON_TO_RY * aseuc.get_masses()

listsp, invsp, idxat = atoms2species(csym)

nsp = len(invsp)


cell = aseuc.get_cell()
from  numpy.linalg import norm
alat = norm(cell[0,:])
rcell = aseuc.get_reciprocal_cell()

fout = open("fc", 'w')

fout.write( """{:6} {:6} {:6} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} \n""".format(nsp, nat, 0, alat/0.52917721067 , 0.,0.,0.,0.,0.) )
fout.write( """{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[0,0]/alat, cell[0,1]/alat, cell[0,2]/alat) )
fout.write( """{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[1,0]/alat, cell[1,1]/alat, cell[1,2]/alat) )
fout.write( """{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[2,0]/alat, cell[2,1]/alat, cell[2,2]/alat) )

for i in range(nsp):
  fout.write( """       {:6} '{:4}' {:15.10f} \n""".format( i+1, csym[invsp[i]], cmass[invsp[i]] ) )

cpos = aseuc.get_positions()
for i in range(nat):
  fout.write( """ {:6} {:6} {:15.10f} {:15.10f} {:15.10f}\n""".format( i+1, idxat[i]+1, cpos[i,0]/alat , cpos[i,1]/alat, cpos[i,2]/alat) )

# no effective charges:
fout.write("F\n")

fin = open("pw.in")
asesc =ase.io.espresso.read_espresso_in(fin)
fin.close()

scell = asesc.get_cell()
nq1 =  int(round(norm(scell[0,:])/norm(cell[0,:])))
nq2 =  int(round(norm(scell[1,:])/norm(cell[1,:])))
nq3 =  int(round(norm(scell[2,:])/norm(cell[2,:])))
ncells = nq1*nq2*nq3
print "Supercell size: ", nq1, nq2, nq3

fout.write("{:6} {:6} {:6} \n".format(nq1, nq2, nq3))

spos = asesc.get_positions()

natsc = len(asesc)

idxsc, idxcell = supercell2cell(nat, cpos, natsc, spos, cell, rcell, nq1,nq2,nq3)

from numpy import zeros
fc = zeros(shape=(3,3,nat,nat,nq1,nq2,nq3))


fin = open("FORCE_CONSTANTS")
natsc_ = int(fin.readline())
if natsc_ != natsc :
  stderr.write( "Wrong number of atoms in FC file" )
  raise Exception 

nread=0
nskip=0
for i in range(natsc):
  idat_i = idxsc[i]
  ci1,ci2,ci3 = idxcell[i]
  for j in range(natsc):
    i_,j_ = [int(x) for x in fin.readline().split()]
    if [i_-1,j_-1] != [i,j]:
      stderr.write( i,j,i_-1,j_-1 )
      stderr.write( "Unexpected atom" )
      raise Exception 

    fxx, fxy, fxz = [float(x) for x in fin.readline().split()]
    fyx, fyy, fyz = [float(x) for x in fin.readline().split()]
    fzx, fzy, fzz = [float(x) for x in fin.readline().split()]
    if ci1 == 0 and ci2 == 0 and ci3 == 0:
      nread += 1
      idat_j = idxsc[j]
      c1,c2,c3 = idxcell[j]
      fc[0][0][idat_i][idat_j][c1][c2][c3] = fxx
      fc[1][0][idat_i][idat_j][c1][c2][c3] = fyx
      fc[2][0][idat_i][idat_j][c1][c2][c3] = fzx
      fc[0][1][idat_i][idat_j][c1][c2][c3] = fxy
      fc[1][1][idat_i][idat_j][c1][c2][c3] = fyy
      fc[2][1][idat_i][idat_j][c1][c2][c3] = fzy
      fc[0][2][idat_i][idat_j][c1][c2][c3] = fxz
      fc[1][2][idat_i][idat_j][c1][c2][c3] = fyz
      fc[2][2][idat_i][idat_j][c1][c2][c3] = fzz
    else:
      nskip += 1
fin.close()

if (nskip+nread)/nread != natsc/nat:
  stderr.write( "Skipped the wrong amount {} {}".format( nskip, nread) )
  raise Exception

for alpha in range(3):
  for beta in range(3):
    for i in range(nat):
      for j in range(nat):
        fout.write( "{:4}{:4}{:4}{:4}\n".format(alpha+1, beta+1, i+1, j+1) )
        for iq3 in range(nq3):
          for iq2 in range(nq2):
            for iq1 in range(nq1):
              fout.write( "{:4}{:4}{:4}{:20.10e}\n".format(iq1+1, iq2+1, iq3+1, fc[alpha][beta][i][j][iq1][iq2][iq3]*F_FACTOR) )

fout.close()
print "Force constants written to file 'fc'"


