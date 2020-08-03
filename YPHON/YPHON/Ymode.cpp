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
  A C++ code to support the script pos2s
  The purpose of this code is to reorganize the output of the ISOTROPY smodes and 
  the ISOTROPY findsym code into a the symmetry adopted basis of the Irrep of phonon modes at 
  Gamma point, and pick out the group operation matrix.
  All the options are defined by pos2s script and users are not needed to prepare inputs
  for this Ymode code
  The standard input of Ymode is the output file "smodesfile.out" the ISOTROPY smodes code
  rotfile is the body part of the output file "findsymout.out" of the ISOTROPY findsym code between
  "_space_group_symop_operation_xyz" and "loop_"
  vecfile is the three lines of the output of the ISOTROPY findsym code after
  "Lattice vectors in cartesian coordinates:"

  By default, the outputs from ISOTROPY codes will be deleted after Ymode exited. User can provide
  the -debug option to pos2s to not delate the  output of the ISOTROPY codes
*/

int foo; //for removing the noising g++ Warnings

double THR = 1.e-4;
int debug = 0;

FILE *iopipe;

/* replacing '\n' with '\000' for a line */
void zeroline(char *line) {
    for (unsigned int i=0; i<strlen(line); i++) {
	if (line[i] != '\n') continue;
	line[i] = '\000';
	break;
    }
}

/* some misc character handling routine */
bool isLower (char ch)
{ // check if a character is lower case (a-z)
    if (ch >= 'a' && ch <= 'z')
        return true;
    else return false;
}
 
bool isUpper (char ch)
{ // check if a character is upper case (A-Z)
    if (ch >= 'A' && ch <= 'Z')
        return true;
    else return false;
}
 
bool isLetter (char ch)
{ // check if a character a character (A-Z || a-z)
    if ((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z'))
        return true;
    else return false;
}
 
bool isDigit (char ch)
{ // check if a character is a number (0-9)
    if (ch >= '0' && ch <= '9')
        return true;
    else return false;
}
 
bool isSymbol (char ch)
{
    if ((ch >= '!' && ch <= '/') || (ch >= ':' && ch <= '@')
      ||(ch >= '[' && ch <= '`') || (ch >= '{' && ch <= '~'))
        return true;
    else return false;
}

#define MAX_STR_LEN 1024

void mproduct(double a[3][3], double b[3][3], double c[3][3])
{
int i, j, k;
double tmp;

        for (i=0; i <3; i++) {
                for (j=0; j <3; j++) {
                        tmp = 0.e0;
                        for (k=0; k<3; k++) tmp += a[i][k]*b[k][j];
                        c[i][j] = tmp;
                }
        }
}

/* routines for misc matrix operations */
void mproduct(double a[3], double b[3][3], double c[3])
{
int i, k;
double tmp;

        for (i=0; i <3; i++) {
                tmp = 0.e0;
                for (k=0; k<3; k++) tmp += a[k]*b[k][i];
                c[i] = tmp;
        }
}

void mproduct(double a[3], double *bb, double c[3])
{
double b[3][3];
int i, k; 
	for (i=0; i <3; i++) 
	    for (k=0; k<3; k++) b[i][k] = bb[i*3+k];
	mproduct (a, b, c);
}   

void mproduct(double a[3][3], double b[3], double c[3])
{
int i, k;
double tmp;

        for (i=0; i <3; i++) {
                tmp = 0.e0;
                for (k=0; k<3; k++) tmp += a[i][k]*b[k];
                c[i] = tmp;
        }
}

void mproductR(double *bb, double a[3], double c[3])
{
double b[3][3];
int i, k;
        for (i=0; i <3; i++)
            for (k=0; k<3; k++) b[i][k] = bb[i*3+k];
        mproduct (b, a, c);
}

/* routines for misc vector operations */
void Anormal(double *a)
{
	double N = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	if (N <= THR) N = 1.e0;
        for (int i=0; i <3; i++) a[i] /= N;
}

double normal(double a)
{
	a = fmod(a, 1.e0);
	if (fabs(a) < THR) a = 0.e0;
	else if (a < 0.e0) a += 1.e0;
	if (fabs(a-1.e0) < THR) a = 0.e0;
	return a;
}

double normal(double *a)
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

double normal2(double *a)
{
	return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
}

double volume(double a[3][3])
{
        return a[0][0]*a[1][1]*a[2][2]
            + a[0][1]*a[1][2]*a[2][0]
            + a[0][2]*a[1][0]*a[2][1]
            - a[0][2]*a[1][1]*a[2][0]
            - a[0][0]*a[1][2]*a[2][1]
            - a[0][1]*a[1][0]*a[2][2];
}

double volume(double *a, double *b, double *c)
{
	return a[0]*b[1]*c[2]
            + a[1]*b[2]*c[0]
            + a[2]*b[0]*c[1]
            - a[2]*b[1]*c[0]
            - a[0]*b[2]*c[1]
            - a[1]*b[0]*c[2];
}

void cross_product(double a[3], double b[3], double c[3])
{
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
}

double dotproduct(double a[3], double b[3])
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

void inv_m(double m[3][3], double mout[3][3]) {
      double d = volume(m[0], m[1], m[2]);
      mout[0][0]=(+m[1][1]*m[2][2]-m[1][2]*m[2][1])/d;
      mout[0][1]=(-m[0][1]*m[2][2]+m[0][2]*m[2][1])/d;
      mout[0][2]=(+m[0][1]*m[1][2]-m[0][2]*m[1][1])/d;
      mout[1][0]=(-m[1][0]*m[2][2]+m[1][2]*m[2][0])/d;
      mout[1][1]=(+m[0][0]*m[2][2]-m[0][2]*m[2][0])/d;
      mout[1][2]=(-m[0][0]*m[1][2]+m[0][2]*m[1][0])/d;
      mout[2][0]=(+m[1][0]*m[2][1]-m[1][1]*m[2][0])/d;
      mout[2][1]=(-m[0][0]*m[2][1]+m[0][1]*m[2][0])/d;
      mout[2][2]=(+m[0][0]*m[1][1]-m[0][1]*m[1][0])/d;
}

void inv_m(double *mm, double *mmout) {
      double m[3][3], mout[3][3];
      for (int i=0; i<3; i++) {
	for (int j=0; j<3;j++) {
	  m[i][j] = mm[i*3+j];
	}
      }
      inv_m(m, mout);
      for (int i=0; i<3; i++) {
	for (int j=0; j<3;j++) {
	  mmout[i*3+j] = mout[i][j];
	}
      }
}

/* class for splitting a line of characters into strings delimited by space, comma, or tab */

class SplitItem {

    private:
    int len, len_p;
    char *list;
    int N;
    char **ListString;

    public:
    SplitItem(char *dlist) { //dlist contains the line of characters
        list = dlist;
        ListString = (char **) malloc( (size_t) (sizeof(char *)) );
        N = 0;
        len = (int) strlen(list);
        len_p = 0;
        while (len_p < len) {
        int x0 = p0();
        int x1 = p1();
        if (x1 == x0) break;
        ListString[N] = (char *) malloc( (size_t) ((x1-x0+1)*sizeof(char)) );
        strncpy(ListString[N], list+x0, x1 - x0);
        ListString[N++][x1 - x0] = '\000';
        ListString = (char **) 
	    realloc(ListString, (size_t) ((N+1)*sizeof(char *)) );
        }
    }

    /* finding the starting point of a field */
    int p0() {
    for (; len_p<len; len_p++) {
        if (list[len_p]!=',' && list[len_p]!=' ' && list[len_p]!='\t') break;
    }
    return len_p;
    }

    /* finding the end point of a field */
    int p1() {
    for (; len_p<len; len_p++) {
        if (list[len_p]==',' ||
                list[len_p]==' ' ||
                list[len_p]=='\t' ||
                list[len_p]=='\n') break;
    }
    return len_p;
    }

    int GetN() {return N;} //return number of fields
    char *str(int i) {return ListString[i];}  //return a field
    int operator[] (int i) {return atoi(ListString[i]);} //convert a field into integer
};

/* class to handle the coefficients of irreducible representations for phonons at the Gamma point. */

class PCOEF {
    private:
	char *line;
	size_t len;
	SplitItem *word;

    public:
	double *c; //save the displacements
	int N;

    PCOEF(int _N) { // N is number of phonon branches - 3 by N (number of atoms in PUC)
	N = _N;
	line = NULL;
	c = new double[N];
	for (int i=0; i<N; i++) c[i] = 0.e0;

        while (1) {
            foo = getline(&line, &len, iopipe);
	    zeroline(line);
            if (!strncmp(line,"--------",8)) { //finding the line of "------------" in smodesfile.out
                break;
            }
        }

        while (1) {
	    long fpos = ftell(iopipe);
            foo = getline(&line, &len, iopipe);
	    zeroline(line);
            if (!strncmp(line,"--------",8)) {
		fseek(iopipe, fpos, SEEK_SET);
                break;
            } else if (!strncmp(line,"********",8)) {
		fseek(iopipe, fpos, SEEK_SET);
                break;
	    } else {
		word = new SplitItem(line);
		double *p = c+3*(atoi(word->str(0))-1);
	        *p++ = atof(word->str(2));
	        *p++ = atof(word->str(3));
	        *p++ = atof(word->str(4));
            }
        }
    }
};

/* 
  class for the vibrational mode analysis, in terms of irreducible representations for phonons at the Gamma 
  point. The outputs is the symmetry.mode file,  containing the information about the group theory 
  representations of different vibrational modes and the The Rotation.sym file, containing the rotation 
  operations of the high symmetry structure from group theory.    
*/

class GammaPhonon {

    private:
	char *line;
	size_t len;
	SplitItem *word;

	double inv_avec[3][3];

    public:
	int Natom;
	char *irrep;
	int IR;
	int Raman;
	int tranMode;
	int Multi;
	int Nmode;
	PCOEF **coef;
	char **sym;
	double *pos;
	double avec[3][3];
	int nTHz;

    GammaPhonon(int _Natom) { //Natom is the number of atoms in primitive unit cell
	line = NULL;
	Natom = _Natom;
	nTHz = 3*Natom;
	IR = 0;
	Raman = 0;
	tranMode = 0;

	while (1) {
            foo = getline(&line, &len, iopipe);
	    zeroline(line);
	    word = new SplitItem(line);
	    if (!strcmp(word->str(0),"Irrep")) { //finding the Irrep like T1g, A2u etc which is the
			//third field of a line in smodesfile.out file initiated by "Irrep"
		irrep = new char[strlen(word->str(2))+1];
		strcpy(irrep,word->str(2));
		break;
	    }
	}


        while (1) {
            foo = getline(&line, &len, iopipe);
	    zeroline(line);
	    word = new SplitItem(line);
	    if (!strcmp(word->str(0), "Degeneracy:")) { //get Degeneracy from smodesfile.out file
		Multi = atoi(word->str(1));
	    } else if (!strncmp(line,"Total number of modes:", 18)) {
		if (strstr(irrep, "*")) Nmode = atoi(word->str(4));
		else Nmode = atoi(word->str(4))/Multi;
	    } else if (!strncmp(line,"Raman active", 12)) {
		Raman = 1;
	    } else if (!strncmp(line,"IR active", 9)) {
		IR = 1;
	    } else if (strstr(line,"translational")) {
		tranMode = 1;
	    } else if (!strncmp(line,"Vectors defining superlattice", 29)) {
		break;
	    }
	}
	//fprintf(stderr, "%4s %3d IR= %d Raman= %d translational= %d\n", irrep, Nmode, IR, Raman, tranMode);

	/* get lattice vector from smodesfile.out file */
	for (int i=0; i<3; i++) {
            foo = getline(&line, &len, iopipe);
	    zeroline(line);
	    word = new SplitItem(line);
	    for (int j=0; j<3; j++) avec[i][j] = atof(word->str(j));
	}

	sym = new char*[Natom];
        foo = getline(&line, &len, iopipe);
	pos = new double[3*Natom];
	/* get atomic positions in the lattice */
	double *p = pos;
	for (int i=0; i<Natom; i++) {
            foo = getline(&line, &len, iopipe);
//printf("%s", line);
	    zeroline(line);
	    word = new SplitItem(line);
	    sym[i] = new char[strlen(word->str(1))+1]; //get atomic symbol
	    strcpy(sym[i], word->str(1)); //get atomic positions
	    *p++ = atof(word->str(2));
	    *p++ = atof(word->str(3));
	    *p++ = atof(word->str(4));
	}

	coef = new PCOEF*[Nmode]; //Nmode is the number of phonon modes at Gamma point

	for (int i=0; i<Nmode; i++) {
	    coef[i] = new PCOEF(nTHz); //read in the displacement pattern of the Irrep
	}

	Normalize(coef[0]->c);

	for (int i=1; i<Nmode; i++) {
	    for (int j=0; j<i; j++) projout(coef[j]->c, coef[i]->c);
	    if (!Normalize(coef[i]->c)) {
		fprintf(stderr, "********FETAL ERROR, can not normalize %s for mode= %i\n", irrep, i);
	    }
	}
    }

    /* serve to othogonalize the coefficients
	gp - a point to a struc of a specivic phonon Irrep
    */
    void overlapping(GammaPhonon *gp) {
	fprintf (stderr, "\n Overlapping between %s and %s\n", irrep, gp->irrep);
	for (int i=0; i<Nmode; i++) {
	    double *pi = coef[i]->c;
	    for (int j=0; j<gp->Nmode; j++) {
		PCOEF *PJ = gp->coef[j];
		double *pj = PJ->c;
	    	double tmp = 0.e0;
	    	for (int x=0; x<nTHz; x++) tmp += pi[x]*pj[x];
		fprintf (stderr, " %6.4lf", tmp);
	    }
	    fprintf (stderr, "\n");
	}
    }

    /* projected out vector pout from  vector p
    */
    void projout(double *pout, double *p) {
	double tmp = 0.e0;
	for (int i=0; i<nTHz; i++) tmp += pout[i]*p[i];
	for (int i=0; i<nTHz; i++) p[i] -= tmp*pout[i];
    }

    /* Normalize the coefficient p
       the original coefs from smodesfile.out file may not normalized or othorgonal each other */
    int Normalize(double *p) {
	double tmp = 0.e0;
	for (int i=0; i<nTHz; i++) tmp += p[i]*p[i];
	if (fabs(tmp) > THR) {
	    tmp = 1.e0/sqrt(tmp);
	    for (int i=0; i<nTHz; i++) p[i] *= tmp;
	    return 1;
	} else return 0;
    }

    /* print out IR, Raman, or translational infor of a mode 
	cfp - output stream to
    */
    void print(FILE *cfp) {
	fprintf(cfp, "%s %d IR= %d Raman= %d translational= %d\n", irrep, Nmode, IR, Raman, tranMode);
	for (int i=0; i<Nmode; i++) {
	    double *p = coef[i]->c;
	    for (int j=0; j<Natom; j++) {
	    	fprintf(cfp, " %20.16lf %20.16lf %20.16lf\n", p[0], p[1], p[2]);
		p += 3;
	    }
	    fprintf(cfp, "\n");
	}
    }

    /* print out IR, Raman, or translational infor of a mode 
	cfp - output stream to
	blankln - an extra line of text
    */
    void printmode(FILE *cfp, const char *blankln) {
        fprintf(cfp, "Irrep %4s Multi= %d ", irrep, Nmode);
        if (IR) fprintf(cfp, " IR active ");
	else if (Raman) fprintf(cfp, " Raman active ");
	else fprintf(cfp, " Silent ");
        if (tranMode) fprintf(cfp, " with one tranlational mode ");
        fprintf(cfp, "%s", blankln);
    }

    /* print out the coeficients for use by Yphon
	cfp - output stream to
	Rot - ratation matrix
	iMap - atomic correspondence between VASP.5 POSCAR and the smodesfile.out file
	Inverse - not USEd
    */
    void print(FILE *cfp, char **posym, double Rot[3][3], int *iMap, double *Inverse) {
        fprintf(cfp, "%s %d IR= %d Raman= %d translational= %d\n", irrep, Nmode, IR, Raman, tranMode);
//	fprintf (stderr,"\n The rotation matrix is:\n");
//	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[0][0], Rot[0][1], Rot[0][2]);
//	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[1][0], Rot[1][1], Rot[1][2]);
//	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[2][0], Rot[2][1], Rot[2][2]);

	printmode(stdout, "\n\n");
        for (int i=0; i<Nmode; i++) {
            double *p = coef[i]->c;
            for (int j=0; j<Natom; j++) {
		double v[3], tmp[3];
		tmp[0] = p[3*iMap[j]];
		tmp[1] = p[3*iMap[j]+1];
		tmp[2] = p[3*iMap[j]+2];
		mproduct(tmp, Rot, v);
		//mproduct(Rot, tmp, v);
                fprintf(cfp, " %9.6lf %9.6lf %9.6lf\n", v[0], v[1], v[2]);
                if (normal(v) > THR) {
		    printf(" %9.6lf %9.6lf %9.6lf %4s %4d\n", v[0], v[1], v[2], posym[j], j);
		}
            }
            fprintf(cfp, "\n");
            printf("\n");
        }
    }

    /* print out POSCAR in VASP.5 format
	cfp - output stream to
    */
    void printPOS(FILE *cfp) {
	fprintf(cfp, "YW POSCAR from Ymode for irrep %s\n 1.00\n", irrep);
	for (int i=0; i<3; i++) {
	    for (int j=0; j<3; j++) fprintf(cfp, " %20.12lf", avec[i][j]);
	    fprintf(cfp, "\n");
	}

	int *n = new int[Natom];
	int N = 1;
	n[0] = 1;

	fprintf(cfp, " %4s", sym[0]);
	for (int i=1; i<Natom; i++) {
	    if (strcmp(sym[i-1], sym[i])) {
		fprintf(cfp, " %4s", sym[i]);
		n[N] = 1;
		N++;
	    } else {
		n[N-1]++;
	    }
	}
	fprintf(cfp, "\n");

	for (int i=0; i<N; i++) fprintf(cfp, " %4d", n[i]);
	fprintf(cfp, "\nD\n");
	inv_m(avec, inv_avec);
	double v[3];

	for (int i=0; i<Natom; i++) {
	    mproduct(pos+3*i, inv_avec, v);
	    for (int j=0; j<3; j++) {
		pos[3*i+j] = v[j];
		fprintf(cfp, " %20.12lf", v[j]);
	    }
	    fprintf(cfp, "\n");
	}
    }
};

/*
   handling the VASP.5 POSCAR file and finding out the rotation symmetries
*/
class RefPOSCAR {

    private:
	double *posM;
	int ii[6];
	double rr[6];
	double gp_avec[3][3];
	double *posgp;
	char **symgp;

    public:
	double Shift[3];
	double avec[3][3], Rot[3][3];
	int iRot[3];
	int Natom;
	char **sym;
	char *line;
	size_t len;
	SplitItem *word;
	double *pos;
	double *posR;
	int *iMap;
	double Inverse[3];

    /* file - the VASP.5 POSCAR file name
       atom - a list of atomic symbols separated by comma representing atomic species follow VASP.5 */
    RefPOSCAR(char *file, char *atom) {
	if (!file) return;
	line = NULL;
	FILE *fpos = fopen(file, "r");
	foo = getline(&line, &len, fpos);
	foo = getline(&line, &len, fpos);
	zeroline(line);
	word = new SplitItem(line);
	/* get lattice vectors */
	double a = atof(word->str(0));
	for (int i=0; i<3; i++) {
            foo = getline(&line, &len, fpos);
	    zeroline(line);
	    word = new SplitItem(line);
	    for (int j=0; j<3; j++) avec[i][j] = a*atof(word->str(j));
//printf("%lf %lf %lfs\n", avec[i][0], avec[i][1], avec[i][2]);
	}

	char *str = NULL;
        foo = getline(&str, &len, fpos);
	zeroline(str);
	SplitItem S(str);
	char *c = S.str(0);

	int NoA = 0;
	if (!isDigit(c[0])) {
            foo = getline(&line, &len, fpos);
	    zeroline(line);
	    atom = str;
	} else {
	    line = str;
  	    if (!atom) {
		atom = line;
		NoA = 1;
	    }
	}

	SplitItem A(atom);
	SplitItem Num(line);

	Natom = 0;
	for (int i=0; i<Num.GetN(); i++) {
	    Natom += Num[i];
	}

	sym = new char* [Natom];
	int n = 0;
	for (int i=0; i<Num.GetN(); i++) {
	    for (int j=0; j<Num[i]; j++) {
		sym[n] = new char[strlen(A.str(i))+4];
		if (NoA) {
		    sprintf(sym[n], "A%d",i+1);
		} else {
		    strcpy(sym[n], A.str(i));
		}
		n++;
	    }
	}
        foo = getline(&line, &len, fpos);

	pos = new double [3*Natom];
	posR = new double [3*Natom];
	posM = new double [3*Natom];
	iMap = new int[Natom];

	double *p = pos;
	for (int i=0; i<Natom; i++) {
            foo = getline(&line, &len, fpos);
	    word = new SplitItem(line);
	    *p++ = atof(word->str(0));
	    *p++ = atof(word->str(1));
	    *p++ = atof(word->str(2));
	}
    }

    /* finding the rotational matrix 
       ISOTROPY smodes made some rational transformation, I can find out that rotation matrix
	gp - a point to a struc of a specivic phonon Irrep
    */
    int Restore(GammaPhonon *gp) {
	double mvec[3][3];
	if (Natom != gp->Natom) { //Natom is the number of atoms in the primitive unit cell 
	    fprintf (stderr,"\n********FETAL ERROR, The input POSCAR is not a primitive unit cell\n");
	    fprintf (stderr,"try use the one from symmetry.pos by pos2s Symmetry.pos -debug\n");
	    fprintf (stderr,"and then run 'mv symmetry.pos p.pos; Ymode <smodesfile.out -POSCAR p.pos \n");
	    exit(1);
	}


	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++) gp_avec[i][j] = gp->avec[i][j];
//	for (int i=0; i<3; i++)
//	    for (int j=0; j<3; j++) fprintf(stderr,"%lf, ", gp_avec[i][j]);
//	fprintf(stderr,"\n");
	posgp = gp->pos;
	symgp = gp->sym;

	if (!FindRM(mvec)) {
	    fprintf (stderr,"\n********FETAL ERROR, CANNOT restore lattice vectors\n");
	    return 0;
	}

	double tmp[3][3];

	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++) mvec[i][j] = Inverse[i]*mvec[i][j];
	inv_m(mvec, tmp);
	mproduct(tmp, avec, Rot);
	//mproduct(avec, tmp, Rot);

	fprintf (stderr,"\n The rotation matrix is:\n");
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[0][0], Rot[0][1], Rot[0][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[1][0], Rot[1][1], Rot[1][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[2][0], Rot[2][1], Rot[2][2]);

	fprintf (stderr,"\n The rotation index are:");
	for (int i=0; i<6; i++) fprintf (stderr, " %4.1lf", rr[i]);
        fprintf (stderr," and Axix relation %d %d %d\n\n", iRot[0], iRot[1], iRot[2]);

//	fprintf (stderr,"\n\n Atomic positios after rotation\n");
//	for (int i=0; i<Natom; i++)
//	    fprintf (stderr,"%10.6lf %10.6lf %10.6lf %s\n", posR[3*i], posR[3*i+1], posR[3*i+2], sym[i]);
	return 1;
    }

    /* finding the rotational matrix index 
 	pos - VASP.5 atomic postion
        sym - VASP.5 atomic symbol
	avec - lattice vector
    */
    int Restore(double *_pos, char **_sym, double _avec[3][3]) {
	double mvec[3][3];

	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++) gp_avec[i][j] = _avec[i][j];
//	for (int i=0; i<3; i++)
//	    for (int j=0; j<3; j++) fprintf(stderr,"%lf, ", gp_avec[i][j]);
//	fprintf(stderr,"\n");
	posgp = _pos;
	symgp = _sym;

	if (!FindRM(mvec)) {
	    fprintf (stderr,"\n********FETAL ERROR, CANNOT restore lattice vectors\n");
	    return 0;
	}

	double tmp[3][3];

	for (int i=0; i<3; i++)
	    for (int j=0; j<3; j++) mvec[i][j] = Inverse[i]*mvec[i][j];
	inv_m(mvec, tmp);
	mproduct(avec, tmp, Rot);
	fprintf (stderr,"\n The POSCAR recovering rotation matrix is:\n");
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[0][0], Rot[0][1], Rot[0][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[1][0], Rot[1][1], Rot[1][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", Rot[2][0], Rot[2][1], Rot[2][2]);

/*
	fprintf (stderr,"\n The POSCAR matrix is:\n");
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", avec[0][0], avec[0][1], avec[0][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", avec[1][0], avec[1][1], avec[1][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", avec[2][0], avec[2][1], avec[2][2]);

	fprintf (stderr,"\n The mvec matrix is:\n");
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", mvec[0][0], mvec[0][1], mvec[0][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", mvec[1][0], mvec[1][1], mvec[1][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", mvec[2][0], mvec[2][1], mvec[2][2]);

	fprintf (stderr,"\n The gp_avec matrix is:\n");
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", gp_avec[0][0], gp_avec[0][1], gp_avec[0][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", gp_avec[1][0], gp_avec[1][1], gp_avec[1][2]);
	fprintf (stderr,"%10.6lf %10.6lf %10.6lf\n", gp_avec[2][0], gp_avec[2][1], gp_avec[2][2]);
*/


	fprintf (stderr,"\n The rotation index are:");
	for (int i=0; i<6; i++) fprintf (stderr, " %4.1lf", rr[i]);
        fprintf (stderr," and Axix relation %d %d %d\n\n", iRot[0], iRot[1], iRot[2]);

//	fprintf (stderr,"\n\n Atomic positios after rotation\n");
//	for (int i=0; i<Natom; i++)
//	    fprintf (stderr,"%10.6lf %10.6lf %10.6lf %s\n", posR[3*i], posR[3*i+1], posR[3*i+2], sym[i]);
	return 1;
    }
   
    /* looping control 
	nn - loop index
    */
    int loop(int nn) {
	int ib = nn,i;
	
	for (i=0; i<6; i++) {
	    int x = ib%5;
	    ii[i] = x;
	    ib /= 5;
	    if (ib==0) break;
	}
	if (i==6) return 0;

	for (i=0; i<6; i++) {
	     if (ii[i] == 0) rr[i] = 0.e0;
	     else if (ii[i] == 1) rr[i] = 1.e0;
	     else if (ii[i] == 2) rr[i] = -1.e0;
	     else if (ii[i] == 3) rr[i] = 2.e0;
	     else if (ii[i] == 4) rr[i] = -2.e0;
	}

	return 1;
    }

    /* finding the rotational matrix
       ISOTROPY smodes made some rational transformation, I can find out that rotation matrix
       mvec - the rresulted transformation matrix
       */
    int FindRM(double mvec[3][3]) {
	double b[3][3], bb[3][3];
	for (int i=0; i<6; i++) {
	    ii[i] = 0;
	    rr[i] = 0.e0;
	}

	int nn = 0;
	while (loop(nn++)) {
		/* reshape the lattice ivectors of the smodes to match them with thos of POSCAR */
		for (int j=0; j<3; j++) {
		    b[0][j] = gp_avec[0][j] + rr[0]*gp_avec[1][j] + rr[1]*gp_avec[2][j];
		    b[1][j] = gp_avec[1][j] + rr[2]*gp_avec[2][j] + rr[3]*gp_avec[0][j];
		    b[2][j] = gp_avec[2][j] + rr[4]*gp_avec[0][j] + rr[5]*gp_avec[1][j];
		}

		if (volume(b) < THR) continue;
		if (FindRM(avec,  b,  mvec)) return 1;

		for (int j=0; j<3; j++) {
		    bb[0][j] = b[0][j];
		    bb[1][j] = -b[1][j];
		    bb[2][j] = -b[2][j];
		}
		if (FindRM(avec,  bb,  mvec)) return 1;

		for (int j=0; j<3; j++) {
		    bb[0][j] = -b[0][j];
		    bb[1][j] = b[1][j];
		    bb[2][j] = -b[2][j];
		}
		if (FindRM(avec,  bb,  mvec)) return 1;

		for (int j=0; j<3; j++) {
		    bb[0][j] = -b[0][j];
		    bb[1][j] = -b[1][j];
		    bb[2][j] = b[2][j];
		}
		if (FindRM(avec,  bb,  mvec)) return 1;
	}
	return 0;
    }

    /* mapping/match atoms between smodesfile.out and Symmetry.pos */
    int mapping() {
	double shift[3];
	if (debug) fprintf (stderr, " mapping() called\n");
	/* compare atomic postion, if OK, print postions out in screen */
	int key = mapping(shift);
        if (!key) return 0;
	mproduct(shift, avec, Shift);

	fprintf (stderr,"\n Atomic mapping with shift %10.6lf %10.6lf %10.6lf Inverse= %10.6lf %10.6lf %10.6lf\n", shift[0], shift[1], shift[2], Inverse[0], Inverse[1], Inverse[2]);
	for (int i=0; i<Natom; i++) {
	    fprintf (stderr," %4d  %10.6lf %10.6lf %10.6lf %4s", i, posgp[3*i], posgp[3*i+1], posgp[3*i+2], symgp[i]);
	    for (int j=0; j<Natom; j++) {
		if (iMap[j]!=i) continue;
		fprintf (stderr," originated from POSCAR atom %4d\n", j);
		break;
	    }
	}
	return 1;
    }

    int mapping(double *s) { //return the translational shift s
	int key = 0;
	for (int i=0; i<Natom; i++) {
	    if (strcmp(sym[0], symgp[i])) continue;
	    for (int j=0; j<3; j++) s[j] = posR[i*3+j] - pos[j];
	    if ((key = mapping(s, symgp))) break;
	}
	return key;
    }

    /* check atomic position matching by symmetry operation 
	s - translation shift
	symgp - Irrep for debug print
    */
    int mapping(double *s, char **symgp) {
	if (debug) fprintf (stderr, " mapping(s, sym) called\n");
        if (debug) {
	    fprintf (stderr,"\n Axix relation %d %d %d\n\n", iRot[0], iRot[1], iRot[2]);
	    for (int i=0; i<6; i++) fprintf (stderr, " %5.1lf", rr[i]);
	    fprintf (stderr, "\n");
	}
	if (debug) {
	    for (int i=0; i<Natom; i++) fprintf (stderr, " %lf %lf %lf %s  %lf %lf %lf %s\n", 
		pos[3*i], pos[3*i+1], pos[3*i+2], sym[i],
		posR[3*i], posR[3*i+1], posR[3*i+2], symgp[i]);
	}
	for (int i=0; i<Natom; i++) {
	    int j;
	    for (j=0; j<Natom; j++) {
	if (debug) fprintf (stderr, "i= %d, j=, %d\n", i, j);
		if (strcmp(sym[i], symgp[j])) continue;
		int k;
		for (k=0; k<3; k++) {
		    double tmp = fabs(pos[i*3+k] + s[k] - posR[j*3+k]);
		    if (debug) fprintf (stderr, "pos i=, %d, %d, %lf ", i, j, tmp);
		    tmp = fmod(tmp, 1.e0);
		    if (debug) fprintf (stderr, "%lf\n", tmp);
		    if (tmp>THR && tmp<1.e0 - THR) break;
		}
		if (k==3) break;
	    }
	    if (j==Natom) return 0;
	    iMap[i] = j;
	    if (debug) fprintf (stderr, " iMap[%d] = %d\n", i,  j);
	}
	return 1;
    }

    /* served to FindRM to find the rotation translation matrix
      to make a match for atoms between smodesfile.out and Symmetry.pos
       ISOTROPY smodes made some rotational transformation, I can find out that rotation matrix.
	For different Irrep, smodes may have made different rotation transformation
	b - smodes lattice vectors in the smodesfile.out file
    */
    int mapping(double b[3][3]) {
		double inv_b[3][3];
                inv_m(b, inv_b);
                for (int i=0; i<Natom; i++) {
                    double v[3];
                    mproduct(posgp+i*3, gp_avec, v);
                    mproduct(v, inv_b, posM+i*3);
                }
                //for (int j=0; j<3*Natom; j++) posM[j] = normal(posM[j]);
                for (int j=0; j<3*Natom; j++) posM[j] = fmod(posM[j],1.e0);

		/* compare atomic postion, if OK, print postions out in screen */
                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = posM[i*3+iRot[0]];
                    posR[i*3+1] = posM[i*3+iRot[1]];
                    posR[i*3+2] = posM[i*3+iRot[2]];
                }
		Inverse[0] = 1.e0;
		Inverse[1] = 1.e0;
		Inverse[2] = 1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = 1.e0 - posM[i*3+iRot[0]];
                    posR[i*3+1] = posM[i*3+iRot[1]];
                    posR[i*3+2] = posM[i*3+iRot[2]];
                }
                Inverse[0] = -1.e0;
                Inverse[1] = 1.e0;
                Inverse[2] = 1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = posM[i*3+iRot[0]];
                    posR[i*3+1] = 1.e0 - posM[i*3+iRot[1]];
                    posR[i*3+2] = posM[i*3+iRot[2]];
                }
                Inverse[0] = 1.e0;
                Inverse[1] = -1.e0;
                Inverse[2] = 1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = posM[i*3+iRot[0]];
                    posR[i*3+1] = posM[i*3+iRot[1]];
                    posR[i*3+2] = 1.e0 - posM[i*3+iRot[2]];
                }
                Inverse[0] = 1.e0;
                Inverse[1] = 1.e0;
                Inverse[2] = -1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = posM[i*3+iRot[0]];
                    posR[i*3+1] = 1.e0 - posM[i*3+iRot[1]];
                    posR[i*3+2] = 1.e0 - posM[i*3+iRot[2]];
                }
                Inverse[0] = 1.e0;
                Inverse[1] = -1.e0;
                Inverse[2] = -1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = 1.e0 - posM[i*3+iRot[0]];
                    posR[i*3+1] = posM[i*3+iRot[1]];
                    posR[i*3+2] = 1.e0 - posM[i*3+iRot[2]];
                }
                Inverse[0] = -1.e0;
                Inverse[1] = 1.e0;
                Inverse[2] = -1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = 1.e0 - posM[i*3+iRot[0]];
                    posR[i*3+1] = 1.e0 - posM[i*3+iRot[1]];
                    posR[i*3+2] = posM[i*3+iRot[2]];
                }
                Inverse[0] = -1.e0;
                Inverse[1] = -1.e0;
                Inverse[2] = 1.e0;
                if (mapping()) return 1;

                for (int i=0; i<Natom; i++) {
                    posR[i*3]   = 1.e0 - posM[i*3+iRot[0]];
                    posR[i*3+1] = 1.e0 - posM[i*3+iRot[1]];
                    posR[i*3+2] = 1.e0 - posM[i*3+iRot[2]];
                }
                Inverse[0] = -1.e0;
                Inverse[1] = -1.e0;
                Inverse[2] = -1.e0;
                if (mapping()) return 1;

		return 0;
    }

    /* finding the rotational matrix
       ISOTROPY smodes made some rotational transformation, I can find out that rotation matrix
       this overlay is to switch the axis of the lattice vectors 
       to match the VASP POSCAR and that rotated ones by smodes
       a - VASP POSCAR lattice vectors
       b - smomdes lattice vectors
       m - the rresulted transformation matrix
       */

    int FindRM(double a[3][3], double b[3][3],  double m[3][3]) {
	if (eq(a[0], a[1], a[2], b[0], b[1], b[2], m)) {
	    iRot[0] = 0;
	    iRot[1] = 1;
	    iRot[2] = 2;
	    if (mapping(b)) return 1;
	}
	if (eq(a[0], a[1], a[2], b[0], b[2], b[1], m)) {
	    iRot[0] = 0;
	    iRot[1] = 2;
	    iRot[2] = 1;
	    if (mapping(b)) return 1;
	}
	if (eq(a[0], a[1], a[2], b[1], b[2], b[0], m)) {
	    iRot[0] = 1;
	    iRot[1] = 2;
	    iRot[2] = 0;
	    if (mapping(b)) return 1;
	}
	if (eq(a[0], a[1], a[2], b[1], b[0], b[2], m)) {
	    iRot[0] = 1;
	    iRot[1] = 0;
	    iRot[2] = 2;
	    if (mapping(b)) return 1;
	}
	if (eq(a[0], a[1], a[2], b[2], b[0], b[1], m)) {
	    iRot[0] = 2;
	    iRot[1] = 0;
	    iRot[2] = 1;
	    if (mapping(b)) return 1;
	}
	if (eq(a[0], a[1], a[2], b[2], b[1], b[0], m)) {
	    iRot[0] = 2;
	    iRot[1] = 1;
	    iRot[2] = 0;
	    if (mapping(b)) return 1;
	}

	return 0;
    }

    /* check if a and be are equal */
    int eq(double a, double b) {
	if (fabs(a-b) < THR)return 1;
	else {
	    if (debug) fprintf (stderr,"a= %10.6lf b=%10.6lf a-b=%10.6lf\n", a, b, fabs(a-b));
	    return 0;
	}
    }

    /* check if vector a and b are in the same length */
    int eq(double *a, double *b) {
//	if (debug) fprintf (stderr,"a %10.6lf %10.6lf %10.6lf\n", a[0], a[1], a[2]);
//	if (debug) fprintf (stderr,"b %10.6lf %10.6lf %10.6lf\n", b[0], b[1], b[2]);
	if (fabs(normal(a)-normal(b)) < THR) {
	    //if (debug) fprintf (stderr,"a= %10.6lf b=%10.6lf\n", fabs(normal(a)-normal(b)), THR);
	    return 1;
	} else return 0;
    }

    /* check if vector a1 and b1 are in the same length */
    int eq(double *a1, double *a2, double *a3, double *b1, double *b2, double *b3, double m[3][3]) {
	if (!eq(a1,b1)) return 0;
	if (!eq(a2,b2)) return 0;
	if (!eq(a3,b3)) return 0;

	if (debug) {
	    fprintf (stderr,"a %10.6lf %10.6lf %10.6lf b %10.6lf %10.6lf %10.6lf\n", a1[0], a1[1], a1[2], b1[0], b1[1], b1[2]);
	    fprintf (stderr,"a %10.6lf %10.6lf %10.6lf b %10.6lf %10.6lf %10.6lf\n", a2[0], a2[1], a2[2], b2[0], b2[1], b2[2]);
	    fprintf (stderr,"a %10.6lf %10.6lf %10.6lf b %10.6lf %10.6lf %10.6lf\n", a3[0], a3[1], a3[2], b3[0], b3[1], b3[2]);
	}

	if (!eq(dotproduct(a1,a2), dotproduct(b1,b2))) return 0;
	if (!eq(dotproduct(a2,a3), dotproduct(b2,b3))) return 0;
	if (!eq(dotproduct(a3,a1), dotproduct(b3,b1))) return 0;
	if (debug) fprintf (stderr,"checking angle done\n");

	for (int j=0; j<3; j++) m[0][j] = b1[j];
	for (int j=0; j<3; j++) m[1][j] = b2[j];
	for (int j=0; j<3; j++) m[2][j] = b3[j];
	return 1;
    }
};

/* class for producing the rotation matrix of the space group */

class RMATRIX {

    private:
	char line[80];
	double (*m)[3][3];
	double (*t)[3];
	int Nrot;
	int natom, nposw;

    public:
	double vectorabc[3][3];
	double prim[3][3];
	double *pos, *posw;
	char **sym, **symw;

	/* print out the rotation and translational operations 
	   Rot - rotation matrix
	   Shift - translation vector
	*/
	void print(double Rot[3][3], double Shift[3]) {
	    double tmp[3][3];
	    mproduct (Rot, vectorabc, tmp);
	    FILE *cfp = fopen("Rotation.sym", "w");
	    for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) fprintf(cfp, " %16.10lf", tmp[i][j]); fprintf(cfp, "\n"); }
	    for (int j=0; j<3; j++) fprintf(cfp, " %16.10lf", Shift[j]); fprintf(cfp, "\n");

	    //int NrotP = 0;
	    //for (int n=0; n<Nrot; n++) if (normal(t[n])<THR) NrotP++;

	    fprintf(cfp, "\n%d\n", Nrot);

	    for (int n=0; n<Nrot; n++) {
		//if (normal(t[n])>THR) continue;
	        for (int i=0; i<3; i++) {
		    for (int j=0; j<3; j++) fprintf(cfp, " %16.10lf", m[n][j][i]); fprintf(cfp, "\n"); }
		for (int j=0; j<3; j++) fprintf(cfp, " %16.10lf", t[n][j]); fprintf(cfp, "\n");
	        fprintf(cfp, "\n");
	    }
	    fclose(cfp);
	    
	}

	/* rotfile is the body part of the output file "findsymout.out" of the ISOTROPY findsym code between
	  "_space_group_symop_operation_xyz" and "loop_"
	  vecfile is the three lines of the output of the ISOTROPY findsym code after
	  "Lattice vectors in cartesian coordinates:" */

	RMATRIX (char *rotfile, char *vecfile) {
	    double slat[3][3];

	    FILE *cfp = fopen(vecfile, "r");
	    if (!cfp) {
		fprintf(stderr, "\n********FETAL ERROR, CANNOT FIND file: %s\n\n", vecfile);
		exit(1);
	    }
	    for (int i=0; i<3; i++) 
		for (int j=0; j<3; j++) foo = fscanf(cfp, "%lf", &prim[i][j]);
	    for (int i=0; i<3; i++) 
		for (int j=0; j<3; j++) foo = fscanf(cfp, "%lf", &slat[i][j]);
	    foo = fscanf(cfp, "%d", &nposw);

	    symw = new char*[nposw];
	    posw = new double[3*nposw];

	    double *p = posw;
	    for (int i=0; i<nposw; i++) {
		char tmp[80], tmp1[80];
		//foo = fscanf(cfp, "%s %s %lf %lf %lf", tmp, tmp1, p++, p++, p++);
		foo = fscanf(cfp, "%s %s %lf %lf %lf", tmp, tmp1, p, p+1, p+2); p += 3;
		symw[i] = new char[strlen(tmp) + 1];
		strcpy (symw[i], tmp);
	        //fprintf(stderr, "input W postion: %lf %lf %lf %s\n", posw[i*3], posw[i*3+1], posw[i*3+2], symw[i]);
	    }

	    fclose(cfp);
	    mproduct(slat, prim, vectorabc);

	    cfp = fopen(rotfile, "r");
	    if (!cfp) {
		fprintf(stderr, "\n********FETAL ERROR, CANNOT FIND file: %s\n\n", rotfile);
		exit(1);
	    }
	    Nrot = 0;
	    while (1) {
		foo = fscanf(cfp, "%s", line);
		if (feof(cfp)) break;
		trim(line);
		if (strlen(line)) Nrot++;
//printf("Nrot=%d, line=<%s>\n", Nrot, line);
	    }
	    fclose(cfp);

	    m = new double[Nrot][3][3];
	    t = new double[Nrot][3];

	    cfp = fopen(rotfile, "r");
	    int index = 0;
	    while (1) {
		foo = fscanf(cfp, "%s", line);
		if (feof(cfp)) break;
		SplitItem Item(line);
//printf ("%s %s %s\n", Item.str(0), Item.str(1), Item.str(2));
		m[index][0][0] = getxyz(Item.str(0), 'x');
		m[index][0][1] = getxyz(Item.str(0), 'y');
		m[index][0][2] = getxyz(Item.str(0), 'z');
		t[index][0] = getxyz(Item.str(0));
//printf ("%s %s %s, m[index][0][0]=%lf %lf %lf %lf\n", Item.str(0), Item.str(1), Item.str(2), m[index][0][0], m[index][0][1], m[index][0][2], t[index][0]);
		m[index][1][0] = getxyz(Item.str(1), 'x');
		m[index][1][1] = getxyz(Item.str(1), 'y');
		m[index][1][2] = getxyz(Item.str(1), 'z');
		t[index][1] = getxyz(Item.str(1));
//printf ("%s %s %s, m[index][0][0]=%lf %lf %lf %lf\n", Item.str(0), Item.str(1), Item.str(2), m[index][1][0], m[index][1][1], m[index][1][2], t[index][1]);
		m[index][2][0] = getxyz(Item.str(2), 'x');
		m[index][2][1] = getxyz(Item.str(2), 'y');
		m[index][2][2] = getxyz(Item.str(2), 'z');
		t[index][2] = getxyz(Item.str(2));
//printf ("%s %s %s, m[index][0][0]=%lf %lf %lf %lf\n", Item.str(0), Item.str(1), Item.str(2), m[index][2][0], m[index][2][1], m[index][2][2], t[index][2]);
//printf ("m[index][0][0]=%lf %lf %lf %lf\n", m[index][0][0], m[index][0][1], m[index][0][2], t[index][0]);
//printf ("m[index][0][0]=%lf %lf %lf %lf\n", m[index][1][0], m[index][1][1], m[index][1][2], t[index][1]);
//printf ("m[index][0][0]=%lf %lf %lf %lf\n", m[index][2][0], m[index][2][1], m[index][2][2], t[index][2]);
//exit(1);
		index++;
	    }
	    fclose(cfp);
	}

	/* making POSCAR file
	   Natom - number of atoms in the PUC
	*/

	void makepos(int Natom) {
	    natom = 0;
	    pos = new double[3*Natom];
	    sym = new char*[Natom];
	    for (int w=0; w<nposw; w++) {
	        for (int r=0; r<Nrot; r++) {
		    double v[3];
		    mproduct(m[r], posw+3*w, v);
		    for (int j=0; j<3; j++) {
			v[j] += t[r][j];
			v[j] = normal(v[j]);
		    }
		    addpos(v, symw[w]);
		    if (natom > Natom) {
			double *p = pos;
			fprintf(stderr, "\n\n********FETAL ERROR, inconsistant wyckoff postion, in =%d, out = %d!\n\n", Natom, natom);
			for (int n=0; n<natom; n++) {
			    fprintf(stderr, " %lf %lf %lf %s\n", p[0], p[1], p[2], sym[n]);
			    p += 3;
			}
			exit(1);
		    }
		}
	    }
	    if (natom != Natom) {
		fprintf(stderr, "\n\n********FETAL ERROR, inconsistant wyckoff postion, in =%d, out = %d!\n\n", Natom, natom);
		exit(1);
	    }
	}

	/* adding new position to the POSCAR 
	  v - a relative atomic position
	  s - atomic symbol
	*/

	void addpos(double *v, char *s) {
	    double tmp[3][3], pp[3], vv[3];
	    mproduct(v, vectorabc, pp);
	    inv_m(prim, tmp);
	    //inv_m(vectorabc, tmp);
	    mproduct(pp, tmp, vv);
	    for (int j=0; j<3; j++) vv[j] = normal(vv[j]);
	    for (int n=0; n<natom; n++) {
	    	double *p = pos+3*n;
		if ((fabs(vv[0]-p[0]) + fabs(vv[1]-p[1]) + fabs(vv[2]-p[2])) < THR) return;
	    }
	    sym[natom] = new char[strlen(s) + 1];
	    strcpy (sym[natom], s);
	    double *p = pos+3*natom;
	    for (int j=0; j<3; j++) p[j] = vv[j];
	    //fprintf(stderr, "Atomic position: %lf %lf %lf %s\n", *p++, *p++, *p++, sym[natom]);
	    natom++;
	}

	/* get the atomic positions from an expression
	  str - an expression like x/y
	  token - 'x', 'y'. or 'z'
        */

	double getxyz(char *str, char token) {
	    int l = strlen (str);
	    int i = 0;
	    for (i=0; i<l; i++) {
		if (str[i]==token) {
//printf("token=%c, str=%s\n", token, str);
		    double val = getxyz(str, i);
		    for (int j=0;j<i;j++) str[j] = ' ';
		    return val;
		}
	    }
	    return 0.e0;
	}

	/* delete the tail space of str */
	void trim(char *str) {
	    for (int i=strlen(str); i >= 0; i--) {
		if (str[i]==' ') str[i] = '\000';
		else break;
	    }
	}

	/* return a double from a "str" with length of "end" */

	double getxyz(char *str, int end) {
	    char s[80];
//printf("ssssstr=%s\n", str);
	    str[end] = ' ';
//printf("ffffstr=%s\n", str);
	    if (end==0) return 1.e0;
	    strncpy(s, str, end);
	    s[end] = '\000';
//printf("s=%s, str=%s\n", s, str);
	    trim(s);
	    if (strlen(s)==0) return 1.e0;
	    else if (s[strlen(s)-1] == '-') return -1.e0;
	    else if (s[strlen(s)-1] == '+') return 1.e0;
	    for (unsigned int i=0; i<strlen(s); i++) {
		if (s[i] == '/') {
		    //s[i] == ' ';
		    s[i] = ' ';
	    	    SplitItem frac(s);
//printf("end= %s, ff=%s\n", frac.str(0), frac.str(1));
	    	    return atof(frac.str(0))/atof(frac.str(1));
		}
	    }
	    return atof(s);
	}

	/* return a double a string expression of "str" */
        double getxyz(char *s) {
            trim(s);
            if (strlen(s)==0) return 0.e0;
            for (unsigned int i=0; i<strlen(s); i++) {
                if (s[i] == '/') {
                    s[i] = ' ';
                    SplitItem frac(s);
//printf("end= %s, ff=%s\n", frac.str(0), frac.str(1));
                    return atof(frac.str(0))/atof(frac.str(1));
                }
            }
            return atof(s);
        }

};

int main(int argc, char* argv[])
{
char *line=NULL;
size_t len;
int Nrep=0, Natom=0;
RefPOSCAR *refposcar;
char* poscar=0;
char* rotfile=0;
char* vecfile=0;
char* atom=0;

        for (int i=0; i<argc; i++) {
            if (!strcmp(argv[i], "-mode")) {
                Nrep = atoi(argv[++i]); //the number of different Irreps at Gamma point
            } else if (!strcmp(argv[i], "-POSCAR")) {
		poscar = argv[++i]; //poscar is the Symmetry.pos file outputed by Yphon in the VASP.5 POSCAR format
            } else if (!strcmp(argv[i], "-Rfile")) {
		rotfile = argv[++i];
            } else if (!strcmp(argv[i], "-Vfile")) {
		vecfile = argv[++i];
            } else if (!strcmp(argv[i], "-atom")) {
		atom = argv[++i]; //ilist of atomic symbols the different atom species follow VASP.5
            } else if (!strcmp(argv[i], "-Natom")) {
                Natom = atoi(argv[++i]); //Number atoms
            } else if (!strcmp(argv[i], "-debug")) {
                debug = 1;
            } else if (!strcmp(argv[i], "-THR")) {
                THR = atof(argv[++i]);
            }
        }

        refposcar = new RefPOSCAR(poscar, atom);

	iopipe = fopen("symmetry.iopipe", "w+");
	while (1) {
	    foo = getline(&line, &len, stdin);
	    if (feof(stdin)) break;
	    fprintf(iopipe, "%s", line);
	}
	fflush(iopipe);
	fseek(iopipe, 0, SEEK_SET);

        /* find the line "atom  type   position" in the "smodesfile.out" file */
	int n = 0;
	while (!Natom) {
	    foo = getline(&line, &len, iopipe);
	    zeroline(line);
	    if (!strcmp(line, "atom  type   position")) break;
	}

	/* find the line started by "Symmetry modes:" */
	while (!Natom) {
	    foo = getline(&line, &len, iopipe);
	    zeroline(line);
	    n++;
	    if (strncmp(line, "Symmetry modes:", 15)) continue;
	    Natom = --n;
	    break;
	}
	fseek(iopipe, 0, SEEK_SET);

	n = 0;
        while (!Nrep) {
            foo = getline(&line, &len, iopipe);
	    if (feof(iopipe)) break;
            if (!strncmp(line, "Irrep GM", 8)) n++;
        }
	if (!Nrep) Nrep = n; //so that we get the right Nrep
	fseek(iopipe, 0, SEEK_SET);

	/* then we can initilize the gp by the number of different Irrep by symmetry */
	GammaPhonon **gp = new GammaPhonon* [Nrep];
	for (int i=0; i<Nrep; i++) {
	    gp[i] = new GammaPhonon(Natom); //gp potinter to each type of Irrep
	    if (feof(iopipe)) break;
	    //for (int j=0; j<=i; j++) gp[i]->overlapping(gp[j]);
	}
	fclose(iopipe);


	FILE *fpos = fopen("symmetry.pos", "w"); //symmetry.pos is for error check
	FILE *fmode = fopen("symmetry.mode", "w"); //"symmetry.mode" is the the file containing the information about the group theory representations of different vibrational modes. 

        printf("\n\n");

	int nn = 0;
	for (int i=0; i<Nrep; i++) {
	    gp[i]->printPOS(fpos);
	    nn += gp[i]->Nmode;
	    int key = 0;
	    if (poscar) key = refposcar->Restore(gp[i]); //check if the rotation matrix is OK for ith Irrep
	    if (!key)gp[i]->print(fmode);
    	    else gp[i]->print(fmode, refposcar->sym, refposcar->Rot, refposcar->iMap, refposcar->Inverse);
	}
	fclose(fmode);
	fclose(fpos);

	fprintf (stderr, "\n Found %d Modes are :", nn);
	for (int i=0; i<Nrep; i++) {
	    if (i==0) fprintf(stderr, " %d%s", gp[i]->Nmode, gp[i]->irrep);
	    else fprintf(stderr, " + %d%s", gp[i]->Nmode, gp[i]->irrep);
	}
	fprintf(stderr, ", with translational modes :");
	for (int i=0; i<Nrep; i++) {
	    if (gp[i]->tranMode) fprintf(stderr, " %s", gp[i]->irrep);
	}
	fprintf(stderr, "\n\n");

	if (!debug) {
	    remove("symmetry.iopipe");
	    remove("symmetry.pos");
	}


	/* rotfile is the body part of the output file "findsymout.out" of the ISOTROPY findsym code between
	  "_space_group_symop_operation_xyz" and "loop_"
	  vecfile is the three lines of the output of the ISOTROPY findsym code after
	  "Lattice vectors in cartesian coordinates:" */
	   
	if (rotfile && vecfile) {
	    RMATRIX *rmatrix = new RMATRIX(rotfile, vecfile);
            if (poscar) {
		rmatrix->makepos(refposcar->Natom);
		refposcar->Restore(rmatrix->pos, rmatrix->sym, rmatrix->prim);
		//refposcar->Restore(rmatrix->pos, rmatrix->sym, rmatrix->vectorabc);
	    }
	    rmatrix->print(refposcar->Rot, refposcar->Shift);
	}

	for (int i=0; i<Nrep; i++) gp[i]->printmode(stderr, "\n");
	fprintf (stderr, "\n found %d atoms and %d irreps\n\n", Natom, Nrep);
}

