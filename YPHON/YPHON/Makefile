.SUFFIXES: .inc .f .f90 .F .c .cc
CC = g++
#CC = c++

YGSL = gsl-1.13
CFLAGS = -Bstatic -O3 -m64 -I ${YGSL}/include -DGSL -Wall
#CFLAGS = -Bstatic -g -m64 -I ${YGSL}/include -DGSL -Wall
GSLLIB = ${YGSL}/lib
CLDFLAGS = ${GSLLIB}/libgsl.a ${GSLLIB}/libgslcblas.a -lm

#EXE = Ycell Yphon vaspfijxml Ymode Ymix yphon YphonN
EXE = Ycell Yphon vaspfijxml Ymode Ymix

all: $(EXE)

YphonQ: YphonQ.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

YphonP: YphonP.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

YphonN: YphonN.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

Yphon: Yphon.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

yphon: yphon.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

Ycell: Ycell.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

Ymode: Ymode.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

vaspfijxml: vaspfijxml.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@

Ymix: Ymix.cpp
	$(CC) $(CFLAGS) $? $(CLDFLAGS) -o $@
