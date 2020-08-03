reset
set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

set output "Exercise_1_8.eps"

qunit=1.0
eunit=1.000000
funit=1.000000
p0 = 0.000000
pp0 = 0.313659
p1 = 0.313659
pp1 = 0.181091
p2 = 0.494750
pp2 = 0.096406
p3 = 0.591156
pp3 = 0.313659
p4 = 0.904815
pp4 = 0.096406
p5 = 1.001221

set xtics ( '{/Symbol G}' qunit*0.000000, 'K' qunit*0.313659, '{/Symbol G}' qunit*0.494750, 'A' qunit*0.591156, 'L' qunit*0.904815, 'M' qunit*1.001221)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p5*1.0001] [funit*0.000000:funit*9.000000] \
'vline.dat' using (qunit*$1):(funit*$2) notitle w l lt 4, \
 'vdis.out' index 0 using (qunit*p0+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$10) notitle w l lt -1, \
 'exp01.dat' index 0 using (qunit*p0+qunit*pp0*($1*2)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 1 using (qunit*p1+qunit*pp1*((0.5-$1)*2)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 2 using (qunit*p2+qunit*pp2*(($1)*2)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 3 using (qunit*p3+qunit*pp3*($1*2)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 4 using (qunit*p4+qunit*pp4*((0.5-$1)*2)):(eunit*$2) notitle w p pt 6 lt 1

