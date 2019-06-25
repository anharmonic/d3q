set xlabel "Temperature (K)"
set ylabel "Volume (a_0^3)" textcolor "black"
set y2label "Volumetric thermal expansion (ppm/K)" textcolor "blue"
set title "Silicon volume/temperature" 
set ytics nomirror

set y2tics auto

p 'qha/k20.dat' w l lw 2 lc "black" not, \
  'qha/k20.dat' u 1:($3*1.e+6) axis x1y2  w l lw 2 lc "blue" not

