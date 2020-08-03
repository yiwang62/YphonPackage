reset
#set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

#set output "vdis.eps"

qunit=1.0
eunit=1.000000
funit=1.000000
p0 = 0.000000
pp0 = 0.266667
p1 = 0.266667
pp1 = 0.377124
p2 = 0.643790
pp2 = 0.230940
p3 = 0.874730

set xtics ( '{/Symbol G}' qunit*0.000000, 'M' qunit*0.266667, '{/Symbol G}' qunit*0.643790, 'R' qunit*0.874730)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p3*1.0001] [funit*0.000000:funit*7.000000] \
'vline.dat' using (qunit*$1):(funit*$2) notitle w l lt 4, \
 'vdis.out' index 0 using (qunit*p0+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$7) notitle w l lt -1, \
 'exp01.dat' index 0 using (qunit*p0+qunit*pp0*($1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 1 using (qunit*p1+qunit*pp1*(1-$1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 2 using (qunit*p2+qunit*pp2*($1)):(eunit*$2) notitle w p pt 6 lt 1

