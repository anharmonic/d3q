set log y
set ylabel "Spectral weight"
set xlabel "Energy cm^{-1}"
p 'LW/spfh_T300_s5.out' u 2:7 w l t 'TO', 'LW/spfh_T300_s5.out' u 2:9 w l t 'LO'
