#!/bin/bash

export EXAMPLE=$(pwd)
export ESPRESSO=$(cd ../../../; pwd)
export ESPRESSO_BIN=$ESPRESSO/bin

pref="mpirun -np 4"

cp input.QHA0 input.QHA

# Run phonon calculation for lattice parameters between 10.00 and 10.40
for i in $(seq 10.00 0.05 10.40);do
  mkdir -p a$i;
 ( cd a$i;
   echo $i;

   sed -e s/__alat/$i/ ../pw.in > pw.in;
   $pref $ESPRESSO_BIN/pw.x -in pw.in > pw.out;

   cp ../ph.in .;
   $pref $ESPRESSO_BIN/ph.x -in ph.in > ph.out;

   rm -rf tmp;

  $ESPRESSO_BIN/d3_q2r.x < ../q2r.in 

  # Add this mat2R file and the electronic energy to the input file
  echo \"a$i/mat2R\" $(grep ^! | '{awk print $5}') >> ../input.QHA 
  );
done

$ESPRESSO_BIN/d3_qha.x 



