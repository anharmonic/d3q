nat3=2*3
# functions to select the column of phono frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

set size ratio 1
unset tics
unset colorbox
my_ps = 2.7
set multiplot layout 2,3

f='LW/lw_NK.12x12x1@bz_T10_strt.out'

set title "ZA"
p f u 3:4:(w(1)) w points palette pt 7 pointsize my_ps not
set title "TA"
p f u 3:4:(w(2)) w points palette pt 7 pointsize my_ps not
set title "LA"
p f u 3:4:(w(3)) w points palette pt 7 pointsize my_ps not
set title "ZO"
p f u 3:4:(w(4)) w points palette pt 7 pointsize my_ps not
set title "TO"
p f u 3:4:(w(5)) w points palette pt 7 pointsize my_ps not
set title "LO"
p f u 3:4:(w(6)) w points palette pt 7 pointsize my_ps not

unset multiplot

