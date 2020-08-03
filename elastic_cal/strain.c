#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void mproduct(a, b, c)
double a[3][3], b[3][3], c[3][3];
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

void normalize(a)
double a[3][3];
{
int i, j;
double tmp;

	tmp = a[0][0]*a[1][1]*a[2][2]
            + a[0][1]*a[1][2]*a[2][0]
            + a[0][2]*a[1][0]*a[2][1]
            - a[0][2]*a[1][1]*a[2][0]
            - a[0][0]*a[1][2]*a[2][1]
            - a[0][1]*a[1][0]*a[2][2];

        tmp = pow (tmp, 1.e0/3.e0);

	for (i=0; i <3; i++) {
		for (j=0; j <3; j++) {
			a[i][j] /= tmp;
		}
	}
}

main(int argc, char **argv)
{
FILE * fp_in = stdin;
char * line = NULL;
size_t len = 0;
int i, j , k;
double ee = 0.01e0;
double avec[3][3], evec[3][3], new_avec[3][3];

	if (argc == 3 ) ee = atof(argv[2]);
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			evec[i][j] = 0.e0;
		}
		evec[i][i] = 1.e0;
	}
	if (!strcmp(argv[1],"e1")) evec[0][0] += ee;
	else if (!strcmp(argv[1],"e2")) evec[1][1] += ee;
	else if (!strcmp(argv[1],"e3")) evec[2][2] += ee;
	else if (!strcmp(argv[1],"e4"))  {
		evec[1][2] += ee/2.0;
		evec[2][1] += ee/2.0;
	} else if (!strcmp(argv[1],"e5"))  {
		evec[0][2] += ee/2.0;
		evec[2][0] += ee/2.0;
	} else if (!strcmp(argv[1],"e6"))  {
		evec[0][1] += ee/2.0;
		evec[1][0] += ee/2.0;
	} else if (!strcmp(argv[1],"e36"))  {
		evec[0][1] += ee/2.0;
		evec[1][0] += ee/2.0;
		evec[2][2] += ee;
	} else if (!strcmp(argv[1],"e45"))  {
		evec[0][2] += ee/2.0;
		evec[2][0] += ee/2.0;
		evec[1][2] += ee/2.0;
		evec[2][1] += ee/2.0;
	} else if (!strcmp(argv[1],"e56"))  {
		evec[0][1] += ee/2.0;
		evec[1][0] += ee/2.0;
		evec[0][2] += ee/2.0;
		evec[2][0] += ee/2.0;
	} else if (!strcmp(argv[1],"e456"))  {
		evec[0][1] += ee/2.0;
		evec[1][0] += ee/2.0;
		evec[0][2] += ee/2.0;
		evec[2][0] += ee/2.0;
		evec[1][2] += ee/2.0;
		evec[2][1] += ee/2.0;
	} 
        getline(&line, &len, fp_in);
	printf("%s", line);
        getline(&line, &len, fp_in);
	printf("%s", line);
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			fscanf(fp_in, "%lf", &avec[i][j]);
		}
	}
	mproduct(avec, evec, new_avec);
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			printf("%18.12f  ",new_avec[i][j]);
		}
		printf("\n");
	}
        getline(&line, &len, fp_in);	
	while (!feof(fp_in)) {
        	getline(&line, &len, fp_in);	
		if (feof(fp_in)) break;
		printf("%s", line);
	}
}
