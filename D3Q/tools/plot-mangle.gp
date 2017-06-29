#!/usr/bin/gnuplot
nat3=2*3
# functions to select the column of phono frequency and linewidth
c(i)=column(i+5)
w(i)=column(i+5+nat3)

wf=1

#set term pngcairo font "Tex Gyre Pagella,36" size 2560,1920
#set term png font "Tex Gyre Pagella,36" size 2560,1920

set style fill transparent solid 0.8 noborder


file="LW/lw_NK_T300_s10.out"
set autoscale fix
unset xtics

 set ar 1 from 1.401188, graph 0 to 1.401188, graph 1 nohead lt -1 lw 3 front 
 set ar 2 from 1.917140, graph 0 to 1.917140, graph 1 nohead lt -1 lw 3 front 
 set ar 3 from 3.279029, graph 0 to 3.279029, graph 1 nohead lt -1 lw 3 front 
 set ar 4 from 4.673116, graph 0 to 4.673116, graph 1 nohead lt -1 lw 3 front 

colors="#000000 #006ddb #920000 #009292 #db6d00 #924900 #490092 #24ff24 #004949 #ff6db6 #b66dff #6db6ff #b6dbff #ffff6d #ffb6db"

factor=50
do for [c=1:nat3]{
 idxa=c+5
 idxw=c+5+nat3
 eval(system('python mangle.py '.file.' 2 '.idxa.' '.idxw.' '.c.' '.factor ))
 set object c polygon fillstyle transparent solid .5 noborder fillcolor rgb word(colors,c)
}

#set term pdfcairo font "Tex Gyre Pagella,16"
#set out "PATH.pdf"
set ylabel "Frequency +/- linewidth (x".factor.")"
p [][:600] \
  for [i=1:nat3] file u 2:(c(i)) w l ls  -1 lw 2    not, \
  0 lt -1 not
#set out

