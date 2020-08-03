#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <complex>
#include <iostream>
#define PI 3.141592653589793

/*
 A C++ code to build supercell in the VASP.5 POSCAR format.

Usage: Ycell [options] <yourprimitiveposcarfile >yoursupercellposcarfile

Where yourprimitiveposcarfile is a file from the standard input and yoursupercellposcarfile is a file for the standard output. Both files are in the VASP.5 POSCAR format.  

The options are:

-bc2
Make a supercell by a kind of doubling the primitive unit cell.

-bc3 
Make a supercell by a kind of tripling the primitive unit cell.

-bc4 
Make a supercell with a size of four times of the primitive unit cell.

-mat matrix 3 by 3 
Make a supercell by transforming the primitive unit cell with a 3 by 3 matrix (9 parameters), for example '-mat 2 -2 0 2 2 0 0 2' where the numbers represent directions in the order of "xx xy xz yx yy yz zx zy zz"

-ss n
Make a n n n  supercell of the primitive unit cell 
*/


using namespace std;

int foo; //for removing the noising g++ Warnings

/* default thresholds */

double ETHR = 1.e-3;
double THR = 1.e-3;
double THR2 = 1.e-10;

double scc[9] = { 1.,0.,0.,0., 1.,0.,0.,0.,1.}; 
double bcc[9] = { 0.,1.,1.,1., 0.,1.,1.,1.,0.}; //double cells
double rho[9] = {1.,-1.,0.,0.,1.,-1.,1.,1.,1.}; //triple cells
double fcc[9] = {-1.,1.,1.,1.,-1.,1.,1.,1.,-1}; //quadraple cells
double tran[3] = { 0.,0.,0.};

/* check if "c" is an number */
int isNum( char c )
{
    if ( c < '0' || c > '9' ) return -1; 
    return c - '0';
} 

double fraction(double xx) 
{
	double fac = 1.e-5;
	for (int i=1; i<21; i++) {
		double x = xx*i;
		double y = rint(x);
		if (abs(y-x)<=fac) return(xx-y/i);
		//if (abs(y-x)<=fac) return(xx-y/i);
	}
	return(0.0);
}

/* misc vector operations */
double *vplus(double *a, double *b)
{
        double *c = new double[3];
	for (int i=0; i<3; i++) c[i] = a[i]+b[i];
	return c;
}

void vRange(double *a, double *b)
{
	for (int i=0; i<3; i++) {
		a[i] = min(a[i],b[i]);
		a[i+3] = max(a[i+3],b[i]);
	}
}

double normal(double a)
{
        a = fmod(a, 1.e0);
        if (fabs(a) < THR*THR) a = 0.e0;
        else if (a < 0.e0) a += 1.e0;
        if (fabs(a-1.e0) < THR*THR) a = 0.e0;
        return a;
}

double normal(double *a)
{
        return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

/* read in a matrix */
void r3(double a[3][3]) {
int i, j;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			foo = scanf("%lf", &a[i][j]);
		}
	}
}

/* read in atomic positions and symbols */
void v3a(double a[3], char *as) {
int i;
	for (i=0; i<3; i++) {
		foo = scanf("%lf", &a[i]);
	}
	foo = scanf("%s", as);
}

/* print out atomic positions and symbols */
void v3w(double a[3], char *as) {
int i;
	for (i=0; i<3; i++) {
		printf("%13.8lf", a[i]);
	}
	printf("  %s\n", as);
}

/* write out a matrix */
void w3(double a[3][3]) {
int i, j;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			printf(" %20.10lf", a[i][j]);
		}
		printf("\n");
	}
}

/* misc matrix operations */
double det(double matrix[3][3]) {
      return
        matrix[0][0] * (matrix[1][1]*matrix[2][2] -matrix[1][2]*matrix[2][1])+
        matrix[0][1] * (matrix[1][2]*matrix[2][0] -matrix[1][0]*matrix[2][2])+
        matrix[0][2] * (matrix[1][0]*matrix[2][1] -matrix[1][1]*matrix[2][0]);
}

int det(int matrix[3][3]) {
      return
        matrix[0][0] * (matrix[1][1]*matrix[2][2] -matrix[1][2]*matrix[2][1])+
        matrix[0][1] * (matrix[1][2]*matrix[2][0] -matrix[1][0]*matrix[2][2])+
        matrix[0][2] * (matrix[1][0]*matrix[2][1] -matrix[1][1]*matrix[2][0]);
}

void inverse_matrix(double a[3][3], double inv[3][3]) {
double d=det(a);
      inv[0][0]=(+a[1][1]*a[2][2]-a[1][2]*a[2][1])/d;
      inv[0][1]=(-a[0][1]*a[2][2]+a[0][2]*a[2][1])/d;
      inv[0][2]=(+a[0][1]*a[1][2]-a[0][2]*a[1][1])/d;
      inv[1][0]=(-a[1][0]*a[2][2]+a[1][2]*a[2][0])/d;
      inv[1][1]=(+a[0][0]*a[2][2]-a[0][2]*a[2][0])/d;
      inv[1][2]=(-a[0][0]*a[1][2]+a[0][2]*a[1][0])/d;
      inv[2][0]=(+a[1][0]*a[2][1]-a[1][1]*a[2][0])/d;
      inv[2][1]=(-a[0][0]*a[2][1]+a[0][1]*a[2][0])/d;
      inv[2][2]=(+a[0][0]*a[1][1]-a[0][1]*a[1][0])/d;
}

void mproduct(double a, double b[3][3], double c[3][3])
{
int i, j;
        for (i=0; i <3; i++) {
                for (j=0; j <3; j++) {
                        c[i][j] = a*b[i][j];
                }
        }
}

void mproduct(double a[3], double b[3][3], double c[3])
{
int i, j;
        for (i=0; i <3; i++) {
		c[i] = 0.0;
                for (j=0; j <3; j++) {
                        c[i] += a[j]*b[j][i];
                }
        }
}


void mproduct(double *a, double b[3][3], double c[3][3])
{
int i, j, k;
double tmp;
        for (i=0; i <3; i++) {
                for (j=0; j <3; j++) {
                        tmp = 0.e0;
                        for (k=0; k<3; k++) {
                                tmp += a[i*3+k]*b[k][j];
                        }
                        c[i][j] = tmp;
                }
        }
}

void mproduct(double a[3][3], double *b, double c[3][3])
{
int i, j, k;
double tmp;
        for (i=0; i <3; i++) {
                for (j=0; j <3; j++) {
                        tmp = 0.e0;
                        for (k=0; k<3; k++) {
                                tmp += a[i][k]*b[k*3+j];
                        }
                        c[i][j] = tmp;
                }
        }
}

void mproduct(double *a, double *b, double c[3][3])
{
int i, j, k;
double tmp;
        for (i=0; i <3; i++) {
                for (j=0; j <3; j++) {
                        tmp = 0.e0;
                        for (k=0; k<3; k++) {
                                tmp += a[i*3+k]*b[k*3+j];
                        }
                        c[i][j] = tmp;
                }
        }
}

void mproduct(double a[3][3], double b[3][3], double c[3][3])
{
int i, j, k;
double tmp;
        for (i=0; i <3; i++) {
                for (j=0; j <3; j++) {
                        tmp = 0.e0;
                        for (k=0; k<3; k++) {
                                tmp += a[i][k]*b[k][j];
                        }
                        c[i][j] = tmp;
                }
        }
}

void copy(double a[3][3], double b[3][3])
{
int i, j;
        for (i=0; i <3; i++)
                for (j=0; j <3; j++)
                        b[i][j] = a[i][j];
}

void copy(double *a, double scale, double b[3][3])
{
	for (int i=0;i<3;i++) {
		for (int j=0;j<3;j++) {
	      		b[i][j] = a[i*3+j]*scale;
		}
	}		
}

void vproduct(double a[3], double b[3][3], double c[3])
{
int i, k;
double tmp;

        for (i=0; i <3; i++) {
                tmp = 0.e0;
                for (k=0; k<3; k++) {
                         tmp += a[k]*b[k][i];
                }
                c[i] = tmp;
        }
}

double volume(double *a, double *b, double *c)
{
double tmp;

        tmp = a[0]*b[1]*c[2]
            + a[1]*b[2]*c[0]
            + a[2]*b[0]*c[1]
            - a[2]*b[1]*c[0]
            - a[0]*b[2]*c[1]
            - a[1]*b[0]*c[2];
        return tmp;
}

int vecm(double *a, double*b, double*c, double *d, double l0)
{
int i, jj=0;
double tmp, e[3];
        for(i=0;i<3;i++) e[i] = a[i] - d[i];
        tmp = volume(e, b, c);
        if (tmp <= 0.e0) return jj;
        tmp = normal(e);
        if (tmp < l0) {
                for(i=0;i<3;i++) a[i] = e[i];
                jj = 1;
        }
        return jj;
}

int vecp(double *a, double*b, double*c, double *d, double l0)
{
int i,jj=0;
double tmp, e[3];
        for(i=0;i<3;i++) e[i] = a[i] + d[i];
        tmp = volume(e, b, c);
        if (tmp <= 0.e0) return jj;
        tmp = normal(e);
        if (tmp < l0) {
                for(i=0;i<3;i++) a[i] = e[i];
                jj = 1;
        }
        return jj;
}

int reduc(double *a, double*b, double*c)
{
double l0;
int jj = 0;
        l0 = normal(a);
        jj += vecm(a, b, c, b, l0);
        l0 = normal(a);
        jj += vecp(a, b, c, b, l0);
        l0 = normal(a);
        jj += vecm(a, b, c, c, l0);
        l0 = normal(a);
        jj += vecp(a, b, c, c, l0);
        return jj;
}

/* class defining the atomic posttions and symbols */

class ATM {

    private:
    
    public:
        double x, y, z;
        char *sym;
        char *comment;
        
        /* set up by the explicit atomic posttions and symbols */
        ATM(double _x, double _y, double _z, char *_sym) {
            x = _x;
            y = _y;
            z = _z;
            int s = strlen(_sym);
            sym = new char[s+1];
            strcpy(sym, _sym);
        }   
        
        /* set up by the explicit atomic posttions and symbols with comment */
        ATM(double _x, double _y, double _z, char *_sym, char *_comment) {
            x = _x;
            y = _y;
            z = _z;
            int s = strlen(_sym);
            sym = new char[s+1];
            strcpy(sym, _sym);
            s = strlen(_comment);
            comment = new char[s+1];
            strcpy(comment, _comment);
        }   

        /* check if var is very small */
        int isq(double var) {
            if (fabs(var)<THR) return 1;
            else return 0;
        }   
        
        /* check if "xx,yy,zz" and "x,y,z" are the same position */
        int operator==(ATM s) {
            double xx = x - s.x;
            double yy = y - s.y;
            double zz = z - s.z;
//            printf ("eq %d\n", isq(xx)+isq(yy)+isq(zz));
            if (isq(xx)+isq(yy)+isq(zz) == 3 ) return 1;
            else return 0;
        }   
};

/* class for splitting a line of characters into strings delimited by space, comma, or tab */

class SplitItem {

    private:
    int len, len_p;
    char *list;
    int N;
    char **ListString;
    char sep;

    public:
    SplitItem(char *dlist, char _sep=',') { //dlist contains the line of characters, default delimiter is ',' */
	sep = _sep;
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
        if (list[len_p]!=',' && 
		list[len_p]!=' ' && 
		list[len_p]!=sep && 
		list[len_p]!='\t') break;
    }
    return len_p;
    }

    /* finding the end point of a field */
    int p1() {
    for (; len_p<len; len_p++) {
        if (list[len_p]==',' ||
                list[len_p]==' ' ||
                list[len_p]=='\t' ||
                list[len_p]==sep ||
                list[len_p]=='\n') break;
    }
    return len_p;
    }
    
    int GetN() {return N;} //return number of fields
    char *str(int i) {return ListString[i];} //return a filed
    int operator[] (int i) {return atoi(ListString[i]);} //convert a field into integer
};

void sccube (double SN, double s[3][3], int *ss) {
double *scell;
	scell = scc;
        int scale = 1;
	//clog <<"SN="<<SN<<"\n";
        if (SN<3.0) {
		scell = bcc;
        } else if (SN<4.0) {
		scell = rho;
        } else if (SN<8.0) {
		scell = fcc;
        } else if (SN<16.0) {
                scale = 2;
		//clog <<"hid SN="<<SN<<"\n";
        } else if (SN<24.0) {
		scell = bcc;
                scale = 2;
        } else if (SN<27.0) {
		scell = rho;
                scale = 2;
        } else if (SN<32.0) {
                scale = 3;
        } else if (SN<54.0) {
		scell = fcc;
                scale = 2;
        } else if (SN<64.0) {
		scell = bcc;
                scale = 3;
        } else if (SN<81.0) {
                scale = 4;
        } else if (SN<108.0) {
		scell = rho;
                scale = 3;
        } else if (SN<125.0) {
		scell = fcc;
                scale = 3;
        } else if (SN<128.0) {
                scale = 5;
        } else if (SN<196.0) {
		scell = bcc;
                scale = 4;
        } else if (SN<216.0) {
		scell = rho;
                scale = 4;
        } else {
                scale = 6;
        }

	copy(scell,scale,s);
	*ss = scale;
	//clog <<"SN="<<SN<<"  scale"<<scale<<"\n";
}

void bccube (double SN, double s[3][3], int *ss) {
        *ss = 1;
	if (SN<2.0) {
		copy(scc,1.0,s);
		return;
	}
	SN /= 2.0;
	double tmp0[3][3];
	sccube(SN,tmp0,ss);
	mproduct(tmp0,bcc,s);
}

void fccube (double SN, double s[3][3], int *ss) {
        *ss = 1;

	if (SN<2.0) {
		copy(scc,1.0,s);
	} else if (SN<4.0) {
		copy(bcc,1.0,s);
	} else {
		SN /= 4.0;
		double tmp0[3][3];
		sccube(SN,tmp0,ss);
		mproduct(tmp0,fcc,s);
	}
//w3(tmp0);
//w3(tmp1);
//w3(s);
}

void rhcube (double SN, double s[3][3], int *ss) {
double *scell;
        *ss = 1;
	double tmp0[3][3];

	if (SN<2.0) {
		scell = scc;
	} else if (SN<3.0) {
		scell = bcc;
	} else if (SN<4.0) {
		scell = rho;
	} else {
		scell = fcc;
		SN /= 4.0;
		sccube(SN,tmp0,ss);
	}
	mproduct(tmp0,scell,s);
}


int main(int argc, char *argv[]) {
FILE * fp = stdin;
char * line = NULL;
size_t len = 0;

double scale;
double alat[3][3], a[3][3], s[3][3], vlat[3][3], inv_a[3][3], inv_vlat[3][3];
double v[3], vout[3], pos[3], tmp[3];
int i, j, k, ss=1, natomS, shift = 1, oshift=1;
char as[128], sym[128], natline[256], v52[256], comment[1024];
char *MAGMOM=0;
int *magmom=0;
double SN,SN0=0.0;
int DN=6;
int esym = 1;
ATM **atomS;
double *scell = scc;
bool cartisian = false;

        atomS = (ATM **) malloc ( (size_t) (sizeof(ATM *)) );
	natomS = 0;

        i = 1;
	
        while (i < argc) {
                if (!strcmp(argv[i], "-ss")) ss = atoi(argv[++i]);
                else if (!strcmp(argv[i], "-bcc")) scell = bcc;
                else if (!strcmp(argv[i], "-bc2")) scell = bcc;
                else if (!strcmp(argv[i], "-rho")) scell = rho;
                else if (!strcmp(argv[i], "-bc3")) scell = rho;
                else if (!strcmp(argv[i], "-fcc")) scell = fcc;
                else if (!strcmp(argv[i], "-bc4")) scell = fcc;
                else if (!strcmp(argv[i], "-noshift")) shift = 0;
                else if (!strcmp(argv[i], "-oshift")) oshift = 0;
                else if (!strcmp(argv[i], "-esym")) esym = 0;
                else if (!strcmp(argv[i], "-cartisian")) cartisian = true;
                else if (!strcmp(argv[i], "-thr")) THR = atof(argv[++i]);
                else if (!strcmp(argv[i], "-SN")) SN0 = atof(argv[++i]);
                else if (!strcmp(argv[i], "-DN")) {
                        DN = atoi(argv[++i]);
			DN = max(1,DN);
                } else if (!strcmp(argv[i], "-ethr")) ETHR = atof(argv[++i]);
                else if (!strcmp(argv[i], "-thr2")) THR2 = atof(argv[++i]);
                else if (!strcmp(argv[i], "-MAGMOM")) MAGMOM = argv[++i];
                else if (!strcmp(argv[i], "-tran")) {
			for (int k=0; k<3; k++) tran[k] = atof(argv[++i]);
		}
                else if (!strcmp(argv[i], "-mat")) {
			scell = new double[9];
			for (int k=0; k<9; k++) scell[k] = atof(argv[++i]);
		}
                i++;
        }

	/* read in POSCAR */
	foo = getline(&line, &len, fp);
	foo = scanf ("%lf",&scale);
	r3(alat);
	mproduct(scale, alat, a);
	inverse_matrix(a, inv_a);

	foo = getline(&line, &len, fp);
	foo = getline(&line, &len, fp);
	SplitItem *el = new SplitItem(line);
	char *str = el->str(0);	
	int key = isdigit(str[0]);
	if (key == 0) {
	    foo = getline(&line, &len, fp);
	}	
	SplitItem *nat = new SplitItem(line);
	int natom = 0;
	for (int i=0; i<nat->GetN(); i++) {
		natom += atoi(nat->str(i));
	}

	/*
	int k = 1;
       	while (k != 0) {
         	k = reduc(a[0], a[1], a[2]);
               	k += reduc(a[1], a[2], a[0]);
               	k += reduc(a[2], a[0], a[1]);
       	}	
	*/

        double SN2 = 0.0;
	if (SN0>0.0) {
            double a0 = normal(a[0]);
            double b0 = normal(a[1]);
            double c0 = normal(a[2]);
            double ang0 = acos((a[1][0]*a[2][0]+a[1][1]*a[2][1]+a[1][2]*a[2][2])/b0/c0)/PI*180.;
            double bng0 = acos((a[2][0]*a[0][0]+a[2][1]*a[0][1]+a[2][2]*a[0][2])/c0/a0)/PI*180.;
            double cng0 = acos((a[0][0]*a[1][0]+a[0][1]*a[1][1]+a[0][2]*a[1][2])/a0/b0)/PI*180.;
	    SN = rint(SN0/(double)natom);
	    SN = max(1.0,SN);
            if (abs(a0-b0)<1.e-4 && abs(a0-c0)<1.e-4) {
		if (abs(ang0-bng0)<1.e-2 && abs(ang0-cng0)<1.e-2) {
			if (abs(ang0-90.) < 1.e-2) { //simple cubic
				sccube (SN, s, &ss);
			} else if (abs(ang0-60.) < 1.e-2) { //fcc
				fccube (SN, s, &ss);
			} else if (abs(ang0-109.4712206) < 1.e-2) {//bcc
				bccube (SN, s, &ss);
			} else {//rhombohedral
				rhcube (SN, s, &ss);
			}
		} else { //might be tetragonal
			SN2 = 1.0;
		}
	    } else {
			SN2 = 1.0;
	    }
/*
 	    if (recalc) {
		double tmp1[3][3];
		natom = (int)(det(s)*natom+0.1);
		//clog <"natom="<<natom<<"\n";
		mproduct (s, a, tmp1);
		copy(tmp1,a);
		inverse_matrix(a, inv_a);
		//clog <<"a0,b0,c0,ang0...=" <<a0 <<" " <<b0 <<" "<<c0 <<" "<<ang0 <<" "<<bng0 <<" "<<cng0 <<" "<<"\n";
        	a0 = normal(a[0]);
        	b0 = normal(a[1]);
        	c0 = normal(a[2]);
        	ang0 = acos((a[1][0]*a[2][0]+a[1][1]*a[2][1]+a[1][2]*a[2][2])/b0/c0)/PI*180.;
        	bng0 = acos((a[2][0]*a[0][0]+a[2][1]*a[0][1]+a[2][2]*a[0][2])/c0/a0)/PI*180.;
        	cng0 = acos((a[0][0]*a[1][0]+a[0][1]*a[1][1]+a[0][2]*a[1][2])/a0/b0)/PI*180.;
	        SN = rint(SN0/(double)natom);
	        SN = max(1.0,SN);
	    } 
*/

	    if (SN2 > 0.0) {
			ss = 1;
			for (int i1=4; i1>0; i1--)
				for (int j1=4; j1>0; j1--)
					for (int k1=4; k1>0; k1--)
						if (SN>=i1*j1*k1) {
							double tmp = i1*j1*k1;
							SN2 = max(tmp,SN2);
			}
	    }
	    //fprintf(stderr,"****SN2= %.3f, SN=%.3f natom=%d\n", SN2, SN, natom);
        }

	if (SN2 == 0.0) {
	    if (SN0==0.0) {
		/* regular mode */
		for (i=0;i<3;i++) {
			for (j=0;j<3;j++) 
			      s[i][j] = scell[i*3+j]*(double)ss;
		}	
	    }
	    mproduct(s, a, vlat);
	}

	if (SN2 > 0.0) {
		/* specified supercel szie mode */
		//inverse_matrix(a, inv_a);
        	double a0 = normal(a[0]);
        	double b0 = normal(a[1]);
        	double c0 = normal(a[2]);
        	double abc000 = max(a0,b0);
        	abc000 = max(abc000,c0);

		SN2 = rint(SN2);
		SN2 = max(1.0,SN2);
		//printf("SN@= %.3f, SN=%.3f\n", SN2, SN);
		DN = min((int)SN2,DN);
                int iSN2 = (int(SN2+0.5));
		//clog <<"iSN2="<<iSN2<<"\n";
	    //fprintf(stderr,"****SN2= %.3f, SN=%.3f\n", SN2, SN);
	    //w3(a);
		double lmax = 1.e36;
		int tmp0[3][3];
		double tmp1[3][3];
		double tmp2[3][3];
		for (int i00=0; i00<=DN; i00++) {
			tmp0[0][0] = i00;
		for (int i01=0; i01<=DN; i01++) {
			tmp0[0][1] = i01;
		for (int i02=0; i02<=DN; i02++) {
			tmp0[0][2] = i02;
			double x1 = i00*a[0][0]+i01*a[1][0]+i02*a[2][0];
			double y1 = i00*a[0][1]+i01*a[1][1]+i02*a[2][1];
			double z1 = i00*a[0][2]+i01*a[1][2]+i02*a[2][2];
			if (sqrt(x1*x1+y1*y1+z1*z1) > 0.1+DN*abc000) continue;
		for (int i10=-DN; i10<=DN; i10++) {
			tmp0[1][0] = i10;
		for (int i11=-DN; i11<=DN; i11++) {
			tmp0[1][1] = i11;
		for (int i12=-DN; i12<=DN; i12++) {
			tmp0[1][2] = i12;
			x1 = i10*a[0][0]+i11*a[1][0]+i12*a[2][0];
			y1 = i10*a[0][1]+i11*a[1][1]+i12*a[2][1];
			z1 = i10*a[0][2]+i11*a[1][2]+i12*a[2][2];
			if (sqrt(x1*x1+y1*y1+z1*z1) > 0.1+DN*abc000) continue;
		for (int i20=-DN; i20<=DN; i20++) {
			tmp0[2][0] = i20;
		for (int i21=-DN; i21<=DN; i21++) {
			tmp0[2][1] = i21;
		for (int i22=-DN; i22<=DN; i22++) {
			tmp0[2][2] = i22;
			if (det(tmp0)!=iSN2) continue;
			x1 = i20*a[0][0]+i21*a[1][0]+i22*a[2][0];
			y1 = i20*a[0][1]+i21*a[1][1]+i22*a[2][1];
			z1 = i20*a[0][2]+i21*a[1][2]+i22*a[2][2];
			if (sqrt(x1*x1+y1*y1+z1*z1) > 0.1+DN*abc000) continue;
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++) tmp1[i][j] = (double)tmp0[i][j];
			mproduct(tmp1, a, tmp2);
			int k = 1;
        		while (k != 0) {
                //printf("k = %d\n", k);
              			k = reduc(tmp2[0], tmp2[1], tmp2[2]);
                		k += reduc(tmp2[1], tmp2[2], tmp2[0]);
                		k += reduc(tmp2[2], tmp2[0], tmp2[1]);
        		}	
			double l0 = normal(tmp2[0]);
			double l1 = normal(tmp2[1]);
			l0 = max(l0,l1);
			l1 = normal(tmp2[2]);
			l0 = max(l0,l1);
			if (l0 < lmax) {
				lmax = l0;
				for (int i=0; i<3; i++)
					for (int j=0; j<3; j++) vlat[i][j] = tmp2[i][j]*(double)ss;
				mproduct(vlat, inv_a, s);
				//mproduct(inv_a, vlat, s);
			}
		}
		}
		}
		}
		}
		}
		}
		}
		}
//w3(s);
        }

	FILE *ftmp = fopen("pMatrix","w");
	for (i=0;i<3;i++)
		for (j=0;j<3;j++) 
		      fprintf (ftmp, " %.0f", s[i][j]);
        fprintf (ftmp, "\n");
        fclose(ftmp);

	inverse_matrix(vlat, inv_vlat);
	printf("Supercell by Yi Wang\n1.00\n");

	if (oshift) {
        	double olat[3][3];
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				olat[i][j] = vlat[i][j];
				if (abs(olat[i][j])<=0.0005) olat[i][j]=0.0;
			}
		w3(olat);
	} else {
		w3(vlat);
	}

	char **symbol = new char*[natom];
	int m = 0;
	if (key==0) {
	    for (int i=0; i<nat->GetN(); i++) {
		int l = atoi(nat->str(i));
		for (int k=0; k<l; k++) symbol[m++] = el->str(i);
	    } 
	} else {
		for (int i=0; i<nat->GetN(); i++) {
			int l = atoi(nat->str(i));
			//itoa(i,as,10);
			sprintf(as,"%d", i);
			for (int k=0; k<l; k++)  {
				symbol[m] = new char[strlen(as)+2];
				strcpy(symbol[m], "A");
				strcat(symbol[m++],as);
			}
	    } 
	}
	
	/* when Mag moments are provided */
	if (MAGMOM!=0) {
	    magmom = new int[natom];
	    SplitItem *mag = new SplitItem(MAGMOM);
	    int k = 0;
	    for (int i=0; i<mag->GetN(); i++) {
		SplitItem *amag = new SplitItem(mag->str(i),'*');
		if (amag->GetN()==1) magmom[k++] = atoi(amag->str(0));
		else {
		    for (int j=0; j<atoi(amag->str(0)); j++) 
			magmom[k++] = atoi(amag->str(1));
		}
	    }
	}

	/* find out supercell boundary in the unit of the input POSCAR */
	double Range[6];
	for (i=0; i<6; i++) Range[i] = 0;
	for (i=0; i<3; i++) vRange(Range, s[i]);
	vRange(Range, vplus(s[0],s[1]));
	vRange(Range, vplus(s[0],s[2]));
	vRange(Range, vplus(s[1],s[2]));
	vRange(Range, vplus(vplus(s[0],s[1]),s[2]));
	int cRange[6];
	for (i=0; i<6; i++) cRange[i] = (int)(rint(Range[i]));

	foo = getline(&line, &len, fp);
	sscanf (line, "%s", as);
	char ch = as[0];
	bool Direct = (ch=='D' || ch=='d');
	m = 0;
	/* find a supercell cube containing the supercell to be built */
	for (int natp=0; natp<natom; natp++) {
		foo = getline(&line, &len, fp);
		if (feof(fp)) break;
		SplitItem *sop = new SplitItem(line);
		if (sop->GetN()<3) continue;
		for (int ii=0; ii<3; ii++) tmp[ii] = atof(sop->str(ii));

		if ( Direct) for (int ii=0; ii<3; ii++) v[ii] = tmp[ii];
		else mproduct (tmp, inv_a, v);
		
		if (sop->GetN()>=4) {
			strcpy(as,sop->str(3));
			strcpy(comment,sop->str(3));
			for (int y=4; y<sop->GetN(); y++) {
				strcat(comment, " ");
				strcat(comment, sop->str(y));
			}
		} else {
			strcpy(as,symbol[m]);
			strcpy(comment,symbol[m++]);
		}

		/* fill the supercell */
		double voutx[3];
		for (i=cRange[0];i<=cRange[3];i++) {
			vout[0] = v[0] + (double)i;
			for (j=cRange[1]; j<=cRange[4]; j++) {
				vout[1] = v[1] + (double)j;
				for (k=cRange[2]; k<=cRange[5]; k++) {
					vout[2] = v[2] + (double)k;
					for (int x=0; x<3; x++) voutx[x] = vout[x] + tran[x];
					vproduct(voutx,a,tmp);
					vproduct(tmp, inv_vlat, pos);
					for (int jj=0; jj<3; jj++) pos[jj] = normal(pos[jj]);
                        		//ATM *tmp = new ATM(pos[0], pos[1], pos[2], as);
                        		ATM *tmp = new ATM(pos[0], pos[1], pos[2], as, comment);
					int m;
					for (m=0; m<natomS; m++) if ((*tmp)==(*atomS[m])) break;
					if (m<natomS) continue;
                        		atomS = (ATM **) realloc (atomS, (size_t) ((natomS+1)*sizeof(ATM *)) );
                        		atomS[natomS++] = tmp;
				}
			}
		}
	}

	/* shift origin */
	if (oshift!=0) {
	    int o = 0;
	    double big = 1e30;
	    for (int i=0; i<natomS; i++) {
		double x = atomS[i]->x;
		double y = atomS[i]->y;
		double z = atomS[i]->z;
		if (x*x + y*y + z*z < big) {
		   o = i;
		   big = x*x + y*y + z*z;
		}
		x = 1.e0-atomS[i]->x;
		y = 1.e0-atomS[i]->y;
		z = 1.e0-atomS[i]->z;
		if (x*x + y*y + z*z < big) {
		   o = i;
		   big = x*x + y*y + z*z;
		}
	    }
	    double ox = atomS[o]->x;
	    double oy = atomS[o]->y;
	    double oz = atomS[o]->z;
	    ox = fraction(ox);
	    oy = fraction(oy);
	    oz = fraction(oz);
	    for (int i=0; i<natomS; i++) {
		atomS[i]->x -= ox;
		atomS[i]->y -= oy;
		atomS[i]->z -= oz;
		if (atomS[i]->x < 0.e0) atomS[i]->x += 1.e0;
		if (atomS[i]->y < 0.e0) atomS[i]->y += 1.e0;
		if (atomS[i]->z < 0.e0) atomS[i]->z += 1.e0;
	    }
	}

	/* small shift of the atoms in the supercell */
	if (shift!=0) {
	    for (int i=0; i<natomS; i++) {
		if (atomS[i]->x > 1.e0-THR) atomS[i]->x -= 1.e0;
		if (atomS[i]->y > 1.e0-THR) atomS[i]->y -= 1.e0;
		if (atomS[i]->z > 1.e0-THR) atomS[i]->z -= 1.e0;
	    }
	    for (int i=0; i<3; i++) tmp[i] = 1.e72;
	    for (int i=0; i<natomS; i++) {
		tmp[0] = min(tmp[0], atomS[i]->x);
		tmp[1] = min(tmp[1], atomS[i]->y);
		tmp[2] = min(tmp[2], atomS[i]->z);
	    }
	    if (tmp[0]<THR) for (int i=0; i<natomS; i++) atomS[i]->x -= tmp[0];
	    if (tmp[1]<THR) for (int i=0; i<natomS; i++) atomS[i]->y -= tmp[1];
	    if (tmp[2]<THR) for (int i=0; i<natomS; i++) atomS[i]->z -= tmp[2];

            for (int i=0; i<3; i++) tmp[i] = -1.e72;
            for (int i=0; i<natomS; i++) {
                tmp[0] = min(tmp[0], atomS[i]->x);
                tmp[1] = min(tmp[1], atomS[i]->y);
                tmp[2] = min(tmp[2], atomS[i]->z);
            }
            if (abs(tmp[0]-1.e0) <THR) for (int i=0; i<natomS; i++) {
		atomS[i]->x -= tmp[0];
                if (atomS[i]->x <0.e0) atomS[i]->x += 1.e0;
	    }
            if (abs(tmp[1]-1.e0) <THR) for (int i=0; i<natomS; i++) {
		atomS[i]->y -= tmp[1];
                if (atomS[i]->y <0.e0) atomS[i]->y += 1.e0;
	    }
            if (abs(tmp[2]-1.e0) <THR) for (int i=0; i<natomS; i++) {
		atomS[i]->z -= tmp[2];
                if (atomS[i]->z <0.e0) atomS[i]->z += 1.e0;
	    }

	    for (int i=0; i<natomS; i++) {
		if (abs(atomS[i]->x-0.25e0) < THR2) atomS[i]->x =0.25e0;
		if (abs(atomS[i]->x-0.75e0) < THR2) atomS[i]->x =0.75e0;
		if (abs(atomS[i]->y-0.25e0) < THR2) atomS[i]->y =0.25e0;
		if (abs(atomS[i]->y-0.75e0) < THR2) atomS[i]->y =0.75e0;
		if (abs(atomS[i]->z-0.25e0) < THR2) atomS[i]->z =0.25e0;
		if (abs(atomS[i]->z-0.75e0) < THR2) atomS[i]->z =0.75e0;
	    }
	}

	/* print out supercell POSCAR */
	ftmp = fopen("t.m.p.POSCAR","w");
	strcpy(natline,"");
	strcpy(v52,"");
	for (int k=0; k<natomS;k++) {
	    strcpy(sym, atomS[k]->sym);
	    if (!strcmp(atomS[k]->sym, "")) continue;
	    strcat(v52, " ");
	    strcat(v52, sym);
	    int nat = 0;
	    for (int i=0; i<natomS;i++) {
		if (!strcmp(sym, atomS[i]->sym)) {
		    tmp[0] = atomS[i]->x;
		    tmp[1] = atomS[i]->y;
		    tmp[2] = atomS[i]->z;
		    if (cartisian) mproduct (tmp, vlat, v);
		    else for (int ii=0; ii<3; ii++) v[ii] = tmp[ii];
		    if (esym==1)
		      fprintf(ftmp, "%15.10lf %15.10lf %15.10lf %s\n", 
			v[0], v[1], v[2], atomS[i]->comment);
		    else
		      fprintf(ftmp, "%15.10lf %15.10lf %15.10lf\n", 
			v[0], v[1], v[2]);
			//atomS[i]->x, atomS[i]->y, atomS[i]->z, atomS[i]->comment);
			//atomS[i]->x, atomS[i]->y, atomS[i]->z, atomS[i]->sym);
	    	    strcpy(atomS[i]->sym, "");
		    nat++;
		}
	    }
	    sprintf(as, " %d", nat);
	    strcat(natline, as);
	}
	fclose(ftmp);
	printf("%s\n", v52);
	if (cartisian) printf("%s\nC\n", natline);
	else printf("%s\nD\n", natline);
	fflush(stdout);

	if (MAGMOM!=0) {
	    int k=natomS/natom;
	    fprintf(stderr, "MAGMOM=");
	    for (int i=0; i<natom; i++)
		fprintf(stderr, " %d*%d", k, magmom[i]);
	    fprintf(stderr, "\n");
	}

	foo = system("cat t.m.p.POSCAR");
	foo = system("\\rm t.m.p.POSCAR");
}

