nat3=2*3
# functions to select the column of phono frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

wf=2

file="LW/lw_NK_T10_s50.out"
set autoscale fix
unset xtics

set ar 1 from 0.5000,graph 0   to 0.5000, graph 1 nohead lt -1 front
set ar 2 from 0.9714,graph 0   to 0.9714, graph 1 nohead lt -1 front

p [][] \
  for [i=1:nat3] file u 2:(c(i)-wf*w(i)):(c(i)+wf*w(i)) w filledcurve not , \
  for [i=1:nat3] file u 2:(c(i)) w l ls  -1 lw 2    not, \
  0 lt -1

