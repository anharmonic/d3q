#!/usr/bin/gnuplot
nat3=2*3
# functions to select the column of phonon frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

set style fill transparent solid 0.8 noborder

file="LW/lw_T300.0_s2.out"

ti=system("echo ".file."| sed -e 's/[-._]/ /g' -e 's/^.*[/]//' -e 's/ out//'")
print ti
set title ti

set autoscale fix
unset xtics

set key opaque

#set ytics 200
set mytics 5

turning=" 0.707107 1.013293 1.931852 2.797877 3.151430 "
nar = words(turning)
do for [a=1:nar]{
 set ar a from word(turning,a), graph 0 to word(turning,a), graph 1 nohead lt -1 
}

turning=" 1.013293 "
nar_old = nar
nar = words(turning)
do for [a=1:nar]{
 set ar a+nar_old from word(turning,a), graph 0 to word(turning,a), graph 1 nohead lt -1 lw 3 front 
}

turning="0.0 0.707107 1.013293 1.931852 2.797877 3.151430 3.504984"
labels="Γ X U-K Γ L W X"
nar = words(turning)
do for [a=1:nar]{
 set label a word(labels,a)  at word(turning,a), character .5 center
}


colors="#006ddb #920000 #009292 #db6d00 #924900 #490092 #24ff24 #004949 #ff6db6 #b66dff #6db6ff #b6dbff #ffff6d #ffb6db #eeeeee"

factor=10
do for [c=nat3:1:-1]{
 idxa=c+5
 idxw=c+5+nat3
 eval(system('python mangle.py '.file.' 2 '.idxa.' '.idxw.' '.c.' '.factor ))
 set object c polygon fillstyle transparent solid .5 noborder fillcolor rgb word(colors,c)
}

#set term pdfcairo font "Tex Gyre Pagella,16"
#set out "PATH.pdf"

#set term pngcairo font "Open Sans,16" size 1024,768
#set out file.".png"
set ylabel "Frequency +/- linewidth (x".factor.")"
p [][:750] \
  for [i=nat3:1:-1] file u 2:(c(i)) w l ls  -1 lc rgb word(colors,i) lw 2    not, \
  0 lt -1 not, 292,395
#set out

