#!/usr/bin.gnuplot -p
set ar 1 from 682./2,0 to 682./2, graph 1 lt 8 nohead
p [200:500] 'LW/final_q0_T300_s5.out' w l t 'Final state energy', \
             NaN lt 8 t 'Initial energy/2'

