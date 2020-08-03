#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <fstream>
#include <iostream>

/*
supporting vasp_fij script
A C++ code to convert the Hessian matrix from the vasprun.xml file into the force constant matrix. vaspfijxml is only called by vasp_fij when it sees the vasprun.xml file in the current folder.
vasp_fij script calls this code.
*/

int foo; //convert a field into integer

/* VASP.5 elemental symbols */

typedef struct {
    char ElementName[20];
} VASPELEMENT;

const char* vaspelement [] = {
                           "Ru","Re","Ra","Rb","Rn","Rh",
                           "Be","Ba","Bi","Bk","Br","H","P",
                            "Os","Hg","Ge","Gd","Ga","Pr","Pt",
                            "Pu","C","Pb","Pa","Pd","Cd","Po",
                            "Pm","Ho","Hf","K","He","Mg","Mo",
                            "Mn","O","S","W","Zn","Eu","Zr",
                            "Er","Ni","Na","Nb","Nd","Ne","Np",
                            "Fr","Fe","B","F","Sr","N","Kr",
                            "Si","Sn","Sm","V","Sc","Sb","Se",
                            "Co","Cm","Cl","Ca","Cf","Ce","Xe",
                            "Lu","Cs","Cr","Cu","La","Li","Tl",
                            "Tm","Th","Ti","Te","Tb","Tc","Ta",
                            "Yb","Dy","I","U","Y","Ac","Ag",
                            "Ir","Am","Al","As","Ar","Au","At",
                            "In"
};

/* VASP.5 atomic mass */

double vaspmass[] = {
                             101.07,186.2,226.0, 85.47, 222.0, 102.90,
                             9.0122, 137.34, 208.98, 247.0,79.91, 1.0079, 30.974,
                             190.20, 200.59, 72.59, 157.25, 69.72, 140.91, 195.09,
                             244.0, 12.01, 207.19, 231.0, 106.40, 112.40, 210.0,
                             145.0, 164.93, 178.49, 39.09, 4.0026, 24.305, 95.94,
                             54.938, 15.999, 32.064, 183.85, 65.38, 151.96, 91.22,
                             167.26, 58.71, 22.9898, 92.91, 144.24, 20.18, 237.05,
                             223.0, 55.85, 10.81, 18.998, 87.62, 14.007, 39.948,
                             28.086, 118.69, 150.35, 50.942, 44.956, 121.75, 78.96,
                             58.93, 247.0, 35.453, 40.08, 251.0, 140.12, 131.30,
                             174.97, 132.91, 52.00, 63.55, 138.91, 6.941, 204.37,
                             168.93, 232.04, 47.90, 127.60, 158.92, 98.91, 180.95,
                             173.04, 162.50, 126.90, 238.03, 88.91, 227.0, 107.87,
                             192.22, 243.0, 26.982, 74.922, 39.948, 196.97, 210.0,
                             114.82
};

/* mode to return the atomic mass by the atomic symbol "e" */

double eMass (char *e) {
    double Mass = -999.0;
    int nElement = sizeof(vaspelement);

    for (int i=0; i<nElement; i++) {
	if (!strcmp(e, vaspelement[i])) {
	    //printf ("Found mass of element %s = %lf\n", vaspelement[i], vaspmass[i]);
	    Mass = vaspmass[i];
	    break;
   	}
    }

    if (Mass <0.e0) {
	printf ("CANNOT FIND mass of element %s, make sure the elemental Symbol correct\n", e);
	exit(1);
    }

    return  Mass;
}

/*
   argv[1] is the number of atoms in the supercell
   argv[2] is file contains the atomic symbols in the supercell
   argv[3] is file contains the Hessians from the vasprun.xml file (body of dynmat)
   All these keys have been defined in the vasp_fij script
*/

int main(int argc, char* argv[])
{
    int natom = atoi(argv[1]);
    FILE *cfp = fopen(argv[2], "r");
    char **eatom = new char* [natom];
    for (int i=0; i<natom; i++) {
	eatom[i] = new char[80];
	foo = fscanf (cfp, "%s", eatom[i]);
    }
    fclose (cfp);

    cfp = fopen(argv[3], "r");
    for (int i=0; i<natom; i++) {
	//printf ( "%d\n", i);
	char* iE = eatom[i];
	double iM = eMass(iE);
	for (int ii=0; ii<3; ii++) {
	    for (int j=0; j<natom; j++) {
		char* jE = eatom[j];
		double fac = sqrt(iM*eMass(jE));
		double hij;
		for (int jj=0; jj<3; jj++) {
		     foo = fscanf (cfp, "%lf", &hij);
		     printf (" %.8lg", hij*fac);
		}
	    }
	    char tmp[80];
	    foo = fscanf (cfp, "%s", tmp);
	    printf ("\n");
	}
    }
    fclose (cfp);
}
