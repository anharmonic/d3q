#!/usr/bin/gnuplot -p
reset

set view map
set pm3d
set cbrange [5.e-6:0.05]
set log cb
set palette defined ( 20 "#ffffff", 30 "orange", 50 "red") 
set xlabel "Path"
set ylabel "Energy (cm^{-1})"

set yrange [0:800]
set ar 1 from   1.000000, 0 to 1.000000, 800 nohead  front
set ar 2 from   1.500000, 0 to 1.500000, 800 nohead  front  dt "-"
set ar 3 from   2.000000, 0 to 2.000000, 800 nohead  front
set ar 4 from   3.414214, 0 to 3.414214, 800 nohead  front

set ar 5 from 0,682. to 4.280239,682 lt 2 nohead front

unset xtics

turning="0.0 1.000000 1.5 2.000000 3.414214 4.280239"
labels = "Γ X W X Γ L"
nar = words(turning)
do for [a=1:nar]{
 set label a word(labels,a)  at word(turning,a), -30 center front
}


sp 'LW/final_q0_qresolved_T300_s5.out' u 1:2:(0):6 w pm3d not,\
    for [i=6:11] 'LW/freq.out' u 2:i:(0) w l lc "black" dt "." not, \

