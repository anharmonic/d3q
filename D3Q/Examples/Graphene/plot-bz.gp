nat3=2*3
# functions to select the column of phono frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

set size ratio 1
unset tics
unset colorbox
my_ps = 2.7
set multiplot layout 2,3

set title "ZA"
p 'LW/lw_NK.12x12x1@bz_T10_s50.out' u 3:4:(w(1)) w points palette pt 7 pointsize my_ps not
set title "TA"
p 'LW/lw_NK.12x12x1@bz_T10_s50.out' u 3:4:(w(2)) w points palette pt 7 pointsize my_ps not
set title "LA"
p 'LW/lw_NK.12x12x1@bz_T10_s50.out' u 3:4:(w(3)) w points palette pt 7 pointsize my_ps not
set title "ZO"
p 'LW/lw_NK.12x12x1@bz_T10_s50.out' u 3:4:(w(4)) w points palette pt 7 pointsize my_ps not
set title "TO"
p 'LW/lw_NK.12x12x1@bz_T10_s50.out' u 3:4:(w(5)) w points palette pt 7 pointsize my_ps not
set title "LO"
p 'LW/lw_NK.12x12x1@bz_T10_s50.out' u 3:4:(w(6)) w points palette pt 7 pointsize my_ps not

unset multiplot

