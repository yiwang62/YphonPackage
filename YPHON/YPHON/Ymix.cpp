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
A C++ code to mix the force constant
*/

/*
   argv[1] is the number of atoms in the supercell
   argv[2] is file contains the atomic symbols in the supercell
   argv[3] is file contains the Hessians from the vasprun.xml file (body of dynmat)
   All these keys have been defined in the vasp_fij script
*/

int foo; //for removing the noising g++ Warnings

double platvec[3][3], slatvec[3][3];

double THR = 0.05e0;

class ATOM {

    private:

    public:
	double x, y, z;
	char s[80];
	int jx;
	double xj, yj, zj;

	static double dmax;
	static int jmax;
	static double xf, yf, zf;

	ATOM() {}
	void set(char *line) {
	    sscanf(line,  "%lf%lf%lf%s", &x, &y, &z, s);
	}

	void jAt(ATOM *spos, int N) {
	    int i = -1;
	    double big = 1.e30;
	    for(i=0; i<N; i++) {
		double xx = x - spos[i].x;
		if (xx <-0.5) xx += 1.0; if (xx >0.5) xx -= 1.0;
		double yy = y - spos[i].y;
		if (yy <-0.5) yy += 1.0; if (yy >0.5) yy -= 1.0;
		double zz = z - spos[i].z;
		if (zz <-0.5) zz += 1.0; if (zz >0.5) zz -= 1.0;
		double dx = xx*slatvec[0][0] + yy*slatvec[1][0] + zz*slatvec[2][0];
		double dy = xx*slatvec[0][1] + yy*slatvec[1][1] + zz*slatvec[2][1];
		double dz = xx*slatvec[0][2] + yy*slatvec[1][2] + zz*slatvec[2][2];
		double dd = sqrt(dx*dx+dy*dy+dz*dz);
		//fprintf(stderr, "dd= %lf\n", dd);
		if (dd < big && !strcmp(s, spos[i].s)) {
		    jx = i;
		    xj = spos[i].x;
		    yj = spos[i].y;
		    zj = spos[i].z;
		    big = dd;
		}
	    }
	    if (big==1.e30) 
		fprintf(stderr, "********Coor not matched for %s at %lf %lf %lf\n", s, x, y, z);
		//fprintf(stderr, "big= %lf\n", big);

	    if (big > dmax) {
		xf = x;
		yf = y;
		zf = z;
		dmax = big;
		jmax = jx;
//fprintf(stderr, "********jmax= %d big=%lf\n", jmax, big);
	    }
	}

	void printpos(double ff1, double ff2, int mode) {
	    double xd = xj - x;
	    if (xd <-0.5) xj += 1.0; if (xd >0.5) xj -= 1.0;
	    double yd = yj - y;
	    if (yd <-0.5) yj += 1.0; if (yd >0.5) yj -= 1.0;
	    double zd = zj - z;
	    if (zd <-0.5) zj += 1.0; if (zd >0.5) zj -= 1.0;
	    if (mode==1) printf (" %12.8lf %12.8lf %12.8lf   %s\n", x, y, z, s);
	    else if (mode==2) printf (" %12.8lf %12.8lf %12.8lf   %s\n", xj, yj, zj, s);
	    else printf (" %12.8lf %12.8lf %12.8lf   %s\n", ff1*x+ff2*xj, ff1*y+ff2*yj, ff1*z+ff2*zj, s);
	}

	static void printmax(ATOM *spos) {
	    fprintf(stderr, "# max atom displacement = %lf  at '%s' No. %d\n      f=%12.8lf %12.8lf %12.8lf  s=%12.8lf %12.8lf %12.8lf\n",
		dmax, spos[jmax].s, jmax, xf, yf, zf, spos[jmax].x,  spos[jmax].y,  spos[jmax].z);
	    if (dmax > 0.5) {
	        fprintf(stderr, "*********be carefule, max atom displacement = %lf is too large,  \n", dmax);
	        fprintf(stderr, "*********be carefule, max atom displacement = %lf is too large,  \n", dmax);
	        fprintf(stderr, "*********be carefule, max atom displacement = %lf is too large,  \n", dmax);
	    }
	}
};

double ATOM::dmax=0.e0;
double ATOM::xf=0.e0;
double ATOM::yf=0.e0;
double ATOM::zf=0.e0;
int ATOM::jmax=0;

int main(int argc, char* argv[])
{
    double ff1=-1.0, ff2=-1.0;
    FILE *f1 = NULL;
    FILE *f2 = NULL;
    double E=0.0, T=300.0;
    double BkB = 8.6173324e-5;
    int alat = 0, blat = 0, mlat = 0;
    double S = 0.e0, NN=-1.e0, T0 = 0.0;
    
        for (int i=1; i<argc; i++) {
            if (!strcmp(argv[i], "-E")) {
                E = atof(argv[++i]); //energy barrier
            } else if (!strcmp(argv[i], "-NN")) {
                NN = atof(argv[++i]); //multiplicity of the second phase
            } else if (!strcmp(argv[i], "-T")) {
                T = atof(argv[++i]); //Temperature
            } else if (!strcmp(argv[i], "-S")) {
                S = atof(argv[++i]); //the number of different Irreps at Gamma point
            } else if (!strcmp(argv[i], "-T0")) {
                T0 = atof(argv[++i]); //the number of different Irreps at Gamma point
            } else if (!strcmp(argv[i], "-f")) {
                ff1 = atof(argv[++i]); //the fraction of the second phase
            } else if (!strcmp(argv[i], "-thr")) {
                THR = atof(argv[++i]); //threshold
            } else if (!strcmp(argv[i], "-alat")) {
                alat = 1; //using the first lattice
            } else if (!strcmp(argv[i], "-blat")) {
                blat = 1; //using the secomd lattice
            } else if (!strcmp(argv[i], "-mlat")) {
                mlat = 1; //using the secomd lattice
	    } else if (f1==NULL) {
		f1 = fopen(argv[i], "r");
		if (f1==NULL) {
		    fprintf(stderr, "\n*******ERROR, CANNOT OPEN FILE '%s'\n", argv[i]);
		    exit(1);
		}
		printf ("# used first file '%s'\n", argv[i]);
	    } else if (f2==NULL) {
		f2 = fopen(argv[i], "r");
		if (f2==NULL) {
		    fprintf(stderr, "\n*******ERROR, CANNOT OPEN FILE '%s'\n", argv[i]);
		    exit(1);
		}
		printf ("# used secnd file '%s'\n", argv[i]);
	    }
        }

    if (ff1 < 0.0 ) {
        if (T <1.e-5) ff1 = 1.0;
        else if ( NN > 0.e0) {
	    double tmp = exp(E/BkB/T)/NN;
	    ff1 = tmp/(1.0+tmp);
    	   printf ("# option: -E %lf -NN %lf -T %lf\n", E, NN, T);
    	   fprintf (stderr, "# option: -E %lf -NN %lf -T %lf\n", E, NN, T);
	} else {
	    //double tmp = exp(-E/kB/T)*exp(S*(1.0-T0/T));
	    double tmp = exp(-T0/T)*exp(S);
	    ff1 = tmp/(1.0+tmp);
	}
    }
    ff2 = 1.0 - ff1;

    fprintf (stderr, "# used mixing parameter ff1= %lf,  ff2= %lf\n", ff1,ff2);
    printf ("# used mixing parameter ff1= %lf,  ff2= %lf\n", ff1,ff2);

    char * line = NULL;
    size_t len = 0;

        //skip comment line leading by '#'
        for (int i=0; ;i++) {
            long tpos = ftell(f1);
            getline(&line, &len, f1);
            if (line[0] != '#') {
                fseek(f1, tpos, SEEK_SET);
                break;
            }
        }

        //skip comment line leading by '#'
        for (int i=0; ;i++) {
            long tpos = ftell(f2);
            getline(&line, &len, f2);
            if (line[0] != '#') {
                fseek(f2, tpos, SEEK_SET);
                break;
            }
        }

    double *pvec;
    for (int i=0; i<6; i++) {
        foo = getline(&line, &len, f1);
	double x1, y1, z1;
	sscanf (line, "%lf%lf%lf", &x1, &y1, &z1);
        foo = getline(&line, &len, f2);
	double x2, y2, z2;
	sscanf (line, "%lf%lf%lf", &x2, &y2, &z2);
	if (i<3) pvec = platvec[i];
	else pvec = slatvec[i-3];
	if (alat) {
	    pvec[0] = x1;
	    pvec[1] = y1;
	    pvec[2] = z1;
	    //printf(" %16.10lf %16.10lf %16.10lf\n", x1,y1,z1);
	} else if (blat) {
	    pvec[0] = x2;
	    pvec[1] = y2;
	    pvec[2] = z2;
	    //printf(" %16.10lf %16.10lf %16.10lf\n", x2,y2,z2);
	} else if (mlat) {
	    pvec[0] = ff1*x1+ff2*x2;
	    pvec[1] = ff1*y1+ff2*y2;
	    pvec[2] = ff1*z1+ff2*z2;
	    //printf(" %16.10lf %16.10lf %16.10lf\n", ff1*x1+ff2*x2,ff1*y1+ff2*y2,ff1*z1+ff2*z2);
	} else {
	    pvec[0] = x1;
	    pvec[1] = y1;
	    pvec[2] = z1;
	    //printf(" %16.10lf %16.10lf %16.10lf\n", x1,y1,z1);
	}
	printf(" %16.10lf %16.10lf %16.10lf\n", pvec[0], pvec[1], pvec[2]);

	//if (alat) printf(" %16.10lf %16.10lf %16.10lf\n", x1,y1,z1);
	//else if (blat) printf(" %16.10lf %16.10lf %16.10lf\n", x2,y2,z2);
	//else if (mlat) printf(" %16.10lf %16.10lf %16.10lf\n", ff1*x1+ff2*x2,ff1*y1+ff2*y2,ff1*z1+ff2*z2);
	//else printf(" %16.10lf %16.10lf %16.10lf\n", x1,y1,z1);
    }

    char * aline = NULL;
    char * bline = NULL;
    foo = getline(&aline, &len, f1);
    foo = getline(&bline, &len, f2);
    int nA, kA;
    sscanf (aline, "%d%d", &nA, &kA);
    int nB, kB;
    sscanf (bline, "%d%d", &nB, &kB);
    int natom, ncell;
    if (nA != nB) {
	fprintf(stderr, "\n *********FETAL ERROR, not matched numbers of atoms between %d and %d\n", nA, nB);
	exit(1);
    } else natom = nA;
    if (alat) {
      ncell = kA;
    } else if (blat) {
      ncell = kB;
    } else if (mlat) {
      if (kA > kB) ncell = kB;
      else ncell = kA;
    } else {
      ncell = kA;
    }

    printf("%d %d\n", natom, ncell);

    foo = getline(&aline, &len, f1);
    foo = getline(&bline, &len, f2);
    printf("%s", aline);

    ATOM *spos = new ATOM[natom];
    ATOM *fpos = new ATOM[natom];

    for (int i=0; i<natom; i++) {
        foo = getline(&line, &len, f1);
	fpos[i].set(line);
        foo = getline(&line, &len, f2);
	spos[i].set(line);
    }

    for (int i=0; i<natom; i++) fpos[i].jAt(spos, natom);

    for (int i=0; i<natom; i++) {
	fpos[i].jAt(spos, natom);
	if (alat) fpos[i].printpos(ff1, ff2, 1);
	else if (blat) fpos[i].printpos(ff1, ff2, 2);
	else fpos[i].printpos(ff1, ff2, -1);
    }
    fpos[0].printmax(spos);

    double *fijf = new double[3*natom*3*natom];
    double *fijs = new double[3*natom*3*natom];

    for (int i=0; i<3*natom; i++) {
        for (int j=0; j<3*natom; j++) {
	    double fc1, fc2;
	    fscanf(f1, "%lf", &fc1);
	    fijf[i*3*natom+j] = fc1;
	    fscanf(f2, "%lf", &fc2);
	    fijs[i*3*natom+j] = fc2;
	}
    }

    for (int i=0; i<3*natom; i++) {
	//int ii = jAt(fpos[i/3], spos, natom);
	int ii = fpos[i/3].jx;
 	    //fprintf(stderr, " i=%d, ii=%d\n", i/3, ii);
        for (int j=0; j<3*natom; j++) {
	    //int jj = jAt(fpos[j/3], spos, natom);
	    int jj = fpos[j/3].jx;
	    double fc1 = fijf[i*3*natom+j];
	    double fc2 = fijs[(3*ii+i%3)*3*natom+3*jj + j%3];
	    //double fc2 = fijs[i*3*natom+j];
	    //fprintf (stderr, "ii=%d, jj=%d\n", i*3*natom+j, (3*ii+i%3)*3*natom+3*jj + j%3);
 	    printf(" %12.6lf", ff1*fc1+ff2*fc2);
	}
 	printf("\n");
    }
}
