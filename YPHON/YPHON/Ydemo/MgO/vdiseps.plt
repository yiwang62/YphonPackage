reset
set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

set output "Exercise_2_8.eps"

qunit=1.0
eunit=1.000000
funit=1.000000
p0 = 0.000000
pp0 = 0.237417
p1 = 0.237417
pp1 = 0.237417
p2 = 0.474834
pp2 = 0.335758
p3 = 0.810592
pp3 = 0.205609
p4 = 1.016201
pp4 = 0.205609
p5 = 1.221810

set xtics ( '{/Symbol G}' qunit*0.000000, 'X' qunit*0.237417, 'X' qunit*0.474834, '{/Symbol G}' qunit*0.810592, 'L' qunit*1.016201, 'X' qunit*1.221810)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p5*1.0001] [funit*0.000000:funit*23.000000] \
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
 'exp01.dat' index 1 using (qunit*p0+qunit*pp0*$1):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 3 using (qunit*p1+qunit*pp1*$1):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 0 using (qunit*p2+qunit*pp2*(1-$1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 2 using (qunit*p3+qunit*pp3*(2*$1)):(eunit*$2) notitle w p pt 6 lt 1

