#!/bin/bash

export EXAMPLE=$(pwd)
export ESPRESSO=$(cd ../../../; pwd)
export ESPRESSO_BIN=$ESPRESSO/bin

pref="mpirun -np 4"

cp input.QHA0 input.QHA

# Run phonon calculation for lattice parameters between 10.00 and 10.40
for i in $(seq 10.0  0.1 10.4);do
  mkdir -p a$i;
 ( cd a$i;
   echo $i;

   sed -e s/__alat/$i/ ../pw.in > pw.in;
   $pref $ESPRESSO_BIN/pw.x -in pw.in > pw.out 2>&1

   cp ../ph.in .;
   $pref $ESPRESSO_BIN/ph.x -in ph.in > ph.out 2>&1

   rm -rf tmp;

  $ESPRESSO_BIN/d3_q2r.x < ../q2r.in > q2r.out 2>&1

  # Add this mat2R file and the electronic energy to the QHA input file
  echo \"a$i/mat2R\" $(grep ^! pw.out| awk '{print $5}') >> ../input.QHA 
  );
done

$ESPRESSO_BIN/d3_qha.x  > qha.out 2>&1

$ESPRESSO_BIN/d3_qha.x --outdir 5kbar --press_kbar 5 > qha_5kbar.out 2>&1

gnuplot -p plot-voft.gp 
gnuplot -p plot-gofv.gp

