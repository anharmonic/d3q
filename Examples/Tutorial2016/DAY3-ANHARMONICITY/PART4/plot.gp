reset

nat3=2*3
#freq
c(i)=column(i+5)
#width
w(i)=column(i+5+nat3)

smw(i)=c(i)-wf*w(i)
spw(i)=c(i)+wf*w(i)


wf=50
file="lw_T300_s20.out"
set autoscale fix

unset xtics
set key center left


plot  \
  for [i=1:nat3] file u 2:(spw(i)):(smw(i)) w filledcurve lt i t ''.i , \
  for [i=1:nat3] file u 2:(c(i)) w l ls -1 lw 1 not


