reset
set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

set output "Exercise_4-1.eps"

qunit=1.0
eunit=1.000000
funit=1.000000
p0 = 0.000000
pp0 = 0.221938
p1 = 0.221938
pp1 = 0.312471
p2 = 0.534408
pp2 = 0.181127
p3 = 0.715535
pp3 = 0.212931
p4 = 0.928467

set xtics ( '{/Symbol G}' qunit*0.000000, 'X' qunit*0.221938, '{/Symbol G}' qunit*0.534408, 'L' qunit*0.715535, 'X' qunit*0.928467)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p4*1.0001] [funit*0.000000:funit*17.000000] \
'vline.dat' using (qunit*$1):(funit*$2) notitle w l lt 4, \
 'vdis.out' index 0 using (qunit*p0+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$16) notitle w l lt -1, \
 'exp05.dat' index 0 using (qunit*p0+qunit*pp0*($1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 3 using (qunit*p0+qunit*pp0*($1)):(eunit*$2) notitle w p pt 7 lt 2, \
 '' index 6 using (qunit*p0+qunit*pp0*($1*2)):(eunit*$2) notitle w p pt 8 lt 3, \
 '' index 1 using (qunit*p1+qunit*pp1*(1-$1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 4 using (qunit*p1+qunit*pp1*(1-$1)):(eunit*$2) notitle w p pt 7 lt 2, \
 '' index 2 using (qunit*p2+qunit*pp2*($1*2)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 5 using (qunit*p2+qunit*pp2*($1)):(eunit*$2) notitle w p pt 7 lt 2, \
 '' index 7 using (qunit*p2+qunit*pp2*($1)):(eunit*$2) notitle w p pt 8 lt 3

