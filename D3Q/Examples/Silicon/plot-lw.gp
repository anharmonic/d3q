nat3=2*3
# functions to select the column of phono frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

# multiply the linewidth by 100 (+50 to -50)
wf=50

file="LW/lw_NK_T300_s10.out"
set autoscale fix
unset xtics

set ar 1 from 1,graph 0   to 1, graph 1 nohead lt -1 front
set ar 2 from 1.7071,graph 0   to 1.7071, graph 1 nohead lt -1 front
set ar 3 from 2.2071,graph 0   to 2.2071, graph 1 nohead lt -1 front
set ar 4 from 3.0731,graph 0   to 3.0731, graph 1 nohead lt -1 front

p [][] \
  for [i=1:nat3] file u 2:(c(i)-wf*w(i)):(c(i)+wf*w(i)) w filledcurve not , \
  for [i=1:nat3] file u 2:(c(i)) w l ls  -1 lw 2    not, \
  0 lt -1 not

