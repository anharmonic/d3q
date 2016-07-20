set encoding iso_8859_15
set terminal postscript enhanced solid color "Helvetica" 20
set output "gnuplot.ps"
#
set key off

set xrange [0:4.280239]
set yrange [0:600]
set arrow from 1,0. to 1,600 nohead  lw 3
set arrow from 2,0. to 2,600 nohead  lw 3
set arrow from 1.5,0. to 1.5,600 nohead  lw 3
set arrow from 3.4142,0. to 3.4142,600 nohead  lw 3
set ylabel "{/Symbol w} (cm^{-1})"
unset xtics
set label "{/Symbol G}" at -0.05,-20
set label "X" at 0.95,-20
set label "W" at 1.45,-20
set label "X" at 1.95,-20
set label "{/Symbol G}" at 3.37,-20
set label "L" at 4.1897,-20

plot "freq.plot" u 1:2 w l lw 2 ; 
replot "frequencies.dat" u 1:2 w p pt 3
