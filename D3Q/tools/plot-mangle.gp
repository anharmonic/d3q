#!/usr/bin/gnuplot
nat3=5*3
# functions to select the column of phonon frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

set style fill transparent solid 0.8 noborder

file="LW/lw2_full_T300_s2.out"

set autoscale fix
unset xtics

#set ytics 200
set mytics 5

turning=" 0.515345 2.054931 2.899002 4.285886 5.899750 6.777038"
nar = words(turning)
do for [a=1:nar]{
 set ar a from word(turning,a), graph 0 to word(turning,a), graph 1 nohead lt -1
}


colors="#006ddb #920000 #009292 #db6d00 #924900 #490092 #24ff24 #004949 #ff6db6 #b66dff #6db6ff #b6dbff #ffff6d #ffb6db #eeeeee"

factor=5
do for [c=1:nat3]{
 idxa=c+5
 idxw=c+5+nat3
 eval(system('python mangle.py '.file.' 2 '.idxa.' '.idxw.' '.c.' '.factor ))
 set object c polygon fillstyle transparent solid .5 noborder fillcolor rgb word(colors,c)
}

#set term pdfcairo font "Tex Gyre Pagella,16"
#set out "PATH.pdf"

#set term pngcairo font "Open Sans,16" size 1024,768
#set out "PATH.png"
set ylabel "Frequency +/- linewidth (x".factor.")"
p [][0:210] \
  for [i=1:nat3] file u 2:(c(i)) w l ls  -1 lc rgb word(colors,i) lw 2    not, \
  0 lt -1 not
#set out

