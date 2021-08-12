#!/usr/bin/gnuplot
reset

set title "Phonon dispersion projected on atomic species"
set autoscale fix

nat=2  # number of atoms
ntyp=2 # number of atomic species
nat3=3*nat


unset xtics

turning="0.0 1.000000 1.5 2.000000 3.414214 4.280239"
labels = "Γ X W X Γ L"
nar = words(turning)
do for [a=1:nar]{
 set label a word(labels,a)  at word(turning,a), -30 center front
}

rgb(r,g,b) = 65536 * int(256*r) + 256 * int(256*g) + int(256*b)

# The file contains nat3 columns with the frequencies,
#  then nat3 columns with the projection on the first atomic type
#  then nat3 columns with the projection on the secodn atomic type
p [][0:800] for [i=1:nat3] 'LW/freq.out' \
      u 2:(column(i+5)):(rgb(\
        column(ntyp*(i-1)+nat3+6),\
        0,\
        column(ntyp*(i-1)+nat3+7))) \
      w l lw 2 lc rgb variable not, \
      NaN lc rgb "#ff0000" lw 2 t 'Mg', \
      NaN lc rgb "#7f007f" lw 2 t ' ', \
      NaN lc rgb "#0000ff" lw 2 t 'O ', \
      0 lt -1 not

      #NaN lc rgb "#555555" lw 2 t '⅓⅓⅓', \
