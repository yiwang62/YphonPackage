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
pp0 = 0.108165
p1 = 0.108165
pp1 = 0.109588
p2 = 0.217753
pp2 = 0.126133
p3 = 0.343886
pp3 = 0.179245
p4 = 0.523131
pp4 = 0.177508
p5 = 0.700640

set xtics ( '{/Symbol G}' qunit*0.000000, 'R' qunit*0.108165, 'D' qunit*0.217753, '{/Symbol G}' qunit*0.343886, 'XX' qunit*0.523131, '{/Symbol G}' qunit*0.700640)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p5*1.0001] [funit*-1.000000:funit*21.000000] \
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
 '' index 0 using (qunit*p0+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$34) notitle w l lt -1, \
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
 '' index 1 using (qunit*p1+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$34) notitle w l lt -1, \
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
 '' index 2 using (qunit*p2+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$34) notitle w l lt -1, \
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
 '' index 3 using (qunit*p3+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$34) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$34) notitle w l lt -1, \
 'exp03.dat' index 0 using (qunit*p0+qunit*pp0*($1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 1 using (qunit*p2+qunit*pp2*(1-$1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 2 using (qunit*p3+qunit*pp3*($1)):(eunit*$2) notitle w p pt 6 lt 1, \
 '' index 3 using (qunit*p4+qunit*pp4*(1-$1)):(eunit*$2) notitle w p pt 6 lt 1

