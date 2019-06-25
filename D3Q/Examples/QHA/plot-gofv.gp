reset

list=system("ls qha/k20_T*.dat|sort -nk1.6")
temperatures=system("ls qha/k20_T*.dat|sort -nk1.6|sed -e s/.*_T// -e s/.dat//")
nt=words(list)

set xlabel "Volume (a_0^3)"
set ylabel "Gibbs free energy (Ry)"
set title "Silicon G(T) at different volumes"
set ytics format "%.3f"
set key center top

set cbrange [0:1000]
set cblabel "Temperature (K)"
set cbtics

p [250:280] for [i=1:nt ] word(list,i) u 1:2 smooth csp lw 0.5 lc palette cb word(temperatures,i) not, \
  'qha/k20.dat' u 2:4 w l lc "black" lw 2 t "V_0(T)"

