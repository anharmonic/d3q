nat3=2*3
# functions to select the column of phono frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

wf=2

f1="LW/lw_NK_T10_s50.out"
f2="LW/lw_NK_T10_strt.out"
set autoscale fix
unset xtics

set ar 1 from 0.577350,graph 0   to 0.577350, graph 1 nohead lt -1 front
set ar 2 from 1.244017,graph 0   to 1.244017, graph 1 nohead lt -1 front

set key bottom right

set style fill transparent solid 0.2 noborder
p [][] \
  for [i=1:nat3] f1 u 2:(c(i)-wf*w(i)):(c(i)+wf*w(i)) w filledcurve ls 2 not , \
  for [i=1:nat3] f2 u 2:(c(i)-wf*w(i)):(c(i)+wf*w(i)) w filledcurve ls 4 not , \
  for [i=1:nat3] f1 u 2:(c(i)) w l ls i lw 2 not, \
  0 lt -1 not, NaN ls 2 w filledcu t 'gauss', NaN ls 4 w filledcu t 'tetra'


