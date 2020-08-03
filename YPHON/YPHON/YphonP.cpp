#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <fstream>
#include <iostream>

/* for mirosoft C++ */
#ifdef WIN32
#define R_OK 1
#include <windows.h>
#endif

#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

/*
Main  C++ code to calculate the phonon density of states and the phonon dispersion based on the force constants calculated by VASP.5. Yphon is the central code of the YPHON package doing the phonon calculation.

The usage of YPHON follows the command style of the Linux operating system. YPHON are composed of a number of Linux commands. A YPHON command is followed by a series of keywords and parameters. Different keywords and parameters are separated by a space character in the command. All keywords in YPHON commands are case sensitive.

Usage: Yphon [options] <superfij.out

Where superfij.out is the name of the file created by vasp_fij containing crystal structure as well as force constant information.  Be default, Yphon calculate the PDOS which is outputted into the file of vdos.out.

The options are:

Note: You may see many other options in the Yphon source code not explained below. The unexplained options are either for the debug purpose or in the test stage.

-nq nqx nqy nqz
Instruct Yphon to calculate the PDOS using a nqx   nqy   nqz mesh in the wave vector space. The default values of nqx nqy nqz are those can provide ~3,000,000 phonon frequencies. Do not be afraid! The high efficiency of Yphon makes the calculation of millions of frequencies done in minutes. 

-DebCut f
Instruct Yphon to fit the PDOS using the Debye expression in the form of C*f**2 for the phonon frequency range lower than the given frequency f. The typical value of f should be around 0.1 ~ 1.0 (in the unit of THz), depending your system and the used q mesh. One should always plot the calculated PDOS and inspect the region at around the frequency f for the smoothness of the curve to decide the proper value of f. This option is useful for the case if one wants to use the PDOS to calculate the heat capacity or Debye temperature at temperature lower than ~10 K, in particular, for superconductor etc.

-Born dielecfij.out 
Consider the vibration-induced dipole-dipole interaction (called LO-TO splitting in the literature) in the calculation, where dielecfij.out is the name of the file created by vasp_BE contianing all information about Born effective charge and high frequency dielectric tensors.

-pdis yourdisfile
This option instructs Yphon to calculate the phonon dispersion instead of the PDOS. yourdisfile is a file defining the directions for the dispersion calculation, see the subsection of "File for dispersion calculation" for instruction on how to prepare the file yourdisfile.

-pvdos 
Calculate also the generalized PDOS (GPDOS, the neutron scattering cross section weighted PDOS), followed Zbiri et al. [21], i.e.,  where  i and pDOSi represent respectively the atomic scattering cross section [22] and the partial phonon density-of-states projected into the individual atoms. The results are saved in the file pvdos.out.

-bvec 
For polar materials only, use lattice vectors the primitive unit cell from the dielecfij.out file to define the wave vector space. 

Note: sometimes, the primitive lattice vectors in the dielecfij.out file may be different from those defined in superfij.out.

 -noNA -nof 
Use together with -bvec. Use the primitive lattice vectors from the dielecfij.out file to define wave vector space, but not consider the effects of the vibration-induced dipole-dipole interaction.

Note: In some cases it makes sense to calculate the phonons of a high symmetry structure using the force constant matrix calculated under a low symmetry structure. For this purpose, one can use a different lattice structure through the dielecfij.out file from that from the superfij.out file. The option -noNA is for the case of conductor where no dielectric information available. The option -nof tells Yphon to ignore some error check steps such as the atomic type mismatch in the case of SQS calculation explained below.

-sqs 
Used together with -noNA -nof. Calculate the phonon dispersions of a random alloy with respect to the wave vector space of the ideal lattice, by averaging over the force constants calculated using a special quasirandom structure (SQS). The detailed formulism can be found in our previous publication.

-mall 
Make an average of the dynamical matrix over all the primitive unit cells within the supercell

When there are no atomic position distortions within the supercell with respect to the primitive unit cell, this option does not change the results, except for taking more computer time. This option is only useful when the supercell is allowed to relax or to calculate the phonon dispersions of high symmetry structure using the force constants calculated from a low symmetry structure. The detailed formulism can be found in our previous publications.

-thr2 parameter 
Define the threshold on how to determine the atomic position relation between the high symmetry structure and the low symmetry structure. Care should be taken in this kind of calculation. One should gradually increase the value of the parameter from 0.01 to 0.15.

-Mass yourmassfile
This option tells Yphon to redefine the atomic mass, being required for the SQS phonon dispersion calculation. The context in the yourmassfile file contains lines like

Cu 96.8975
Au 96.8975

-Gfile  symmetry.mode
This option tells Yphon to make the vibrational mode analysis, in terms of irreducible representations for phonons at the   point. The symmetry.mode file must be made by using the enclosed script pos2s, containing the information about the group theory representations of different vibrational modes.   

-Rfile  Rotation.sym
In some cases it makes sense to calculate the phonons of a high symmetry structure using the force constant matrix calculated under a low symmetry structure. This option tells Yphon to "restore" the symmetry of the high symmetry crystal


Note: Caution should be taken for this kind of calculation. In particular, large supercell is needed to delimit the effects of symmetry loss of the supercell. 

-plot 
Instructs Yphon to display the plot in the terminal using gnuplot for one to check the calculated results.

-expt exp01.dat 
Instructs Yphon to plot the experimental data contained in the file "exp01.dat" together with the calculations. 

*/

/*
To calculate the phonon dispersion, one is needed to use the command line option of "-pdis yourdisfile" with Yphon. yourdisfile is a file defining the direction for the dispersion calculation, and typically its format is (from our unpublished calculation of GaAs)

0 0 0 0 0 .5 Gamma X 0 $1 2 0 $1 3
0 .5 .5 0 0 0 X Gamma 1 (1.-$1) 2 1 (1.-$1) 3 1 (1.-$1) 4
0 0 0 .25 .25 .25 Gamma L 2 (2*$1) 2 2 (2*$1) 3
.25 .25 .25 0 0 .5 L X 3 (2*$1) 2 3 (2*$1) 3 3 (2*$1) 4
0 0 .5  0 .25 .5 X W 4 (2*$1) 2 4 (2*$1) 3 4 (2*$1) 4
0 .25 .5 .25 .25 .25 W L 5 (2*$1) 2 5 (2*$1) 3 5 (2*$1) 4

Data column 1-3: the reciprocal reduced coordinate of the starting q point along the dispersion path 
Data column 4-6: the reciprocal reduced coordinate of the end q point along the dispersion path
Data column 7: the label of the starting q point along the dispersion path
Data column 8: the label of the end q point along the dispersion path
Data column 9-: multi sets of data each set containing three columns wherein the first column is the index of the data group (separated by two blank lines following the convention of the gnuplot) in the experimental data file, the second column tells gnuplot which column and how to transform the column into the q point following the convention of the gnuplot, and the third column tells gnuplot  which column of experimental data will be used as the frequency data. This can help you save a lot of time if you can learn it. The calculated phonon dispersions are contained in the file of vdis.out.

Note: care should be taken about the suffix of the yourdisfile file. The suffix .fcc, .bcc, .hcp, and .tet2 are reserved for the fcc, bcc, hcp, and tetragonal  crystals only. For these crystal, Yphon internally converts the q vector of the primitive unit cell into that of the conventional unit cell using the following C++ statements

double bcc[9] = {-1.,1.,1.,1.,-1.,1.,1.,1.,-1};
double fcc[9] = { 0.,1.,1.,1., 0.,1.,1.,1.,0.};
double hcp[9] = {1., 0.,0.,1.,-1.,0.,0.,0.,1.};
double tet2[9] = {-1.,1.,1.,1.,-1.,1.,1.,1.,-1};

	For the cases of bcc and fcc crystals, most neutron scattering data are reported with respect to the cubic conventional cell instead of the primitive unit cell. Yphon did the conversion internally, assumed that you defined the shape of the primitive unit cell of the fcc crystal in the POSCAR file as

0 .5 .5
.5 0 .5
.5 .5 0 

and you defined the shape of the primitive unit cell of the bcc crystal in the POSCAR file as

-.5 .5 .5
.5 -.5 .5
.5 .5 -.5

	For the case of hcp crystal, Yphon assumed that you defined the shape of the primitive unit cell of the crystal in the POSCAR file as (note that c in the third line below is the relative lattice parameter in c direction)

0.8660254037844 -.5 0.
0.0000000000000  1. 0.
0.0000000000000  0. c

	For the case of tetragonal crystal, Yphon assumed that you defined the shape of the primitive unit cell of the crystal in the POSCAR file as

1.0	0. 0.
0.	1.0 0.
0.5 0.5 c

	If you do not use the suffix .fcc, .bcc, .hcp, and .tet2 for the yourdisfile file, Yphon will define the direction of wave vector using the reciprocal lattice vector of the primitive unit cell. The reciprocal lattice vector is printed out in the screen as the last three lines like

   -0.000000000000     -0.362182366065     -0.000000000000  
   -0.313659560187     -0.181091751042     -0.000000000000  
   -0.000000000000     -0.000000000000     -0.192811598003  

	You can refer these (only the direction is important) to define the direction of your phonon dispersion calculation.

	For the dispersion calculation, we strongly recommend one refer the web site http://www.cryst.ehu.es/ for the definition of the KVEC, i.e., the k-vector types and Brillouin zones of the space groups.
*/

using namespace std;

int foo; //for removing the noising g++ Warnings

const char *ostype; //type of oerating system
const char *shell; //linux shell type

/* to make a plot using gnuplot 
	plt -  the name of a  file containing the gnuplot script */

void ynuplot (const char *plt) {
    if (!strcmp(ostype,"windows")) return;
    char *cmd = new char[256];
    try {
        strcpy(cmd, shell);
        strcat(cmd, " -c \"gnuplot -persist ");
        strcat(cmd, plt);
        strcat(cmd, "\"");
        foo = system(cmd);
        printf ("run gnuplot on %s ON under %s\n", ostype, cmd);
    }
    catch (char * str ) {
      cout << "Exception raised in ynuplot: " << str << '\n';
    }
}

//#define _hypot hypot

/* to check if a file named "pdisfile" existed (for compatibility with microsoft C++) 
  the input "key" is not used */

int access(char *pdisfile, int key)
{
	FILE *cfp = fopen(pdisfile, "r");

	if (cfp) {
		fclose(cfp);
		return 0;
	} else {
		return 1;
	}
}

/* not used when using GSL, for using with ISML or INTEL linear algebra solver */

#ifdef IMSL
#include <imsl.h>
#endif

#ifdef INTEL
#include <mkl_lapack.h>
#endif


typedef complex<double> dcmplx;

/* default thresholds */

double THR = 1.e-4;
double THR2 = 1.e-4;
double THR3 = 1.e-4;
double THRE = 5.e-3;
int NRDIS = 10;
int DMODE = 0;
int CTLDM=1;

double inv_dielec[9];

/* some physical constants */
#define PI 3.141592653589793238462643e0
#define PI2 (2.e0*3.141592653589793238462643e0)
#define meVtoTHz 0.2417991077061 /* 1.60217733e-19/6.626068e-34/1.e15 */
#define LIGHTSPEED 299792458.e0
#define EVTOJ 1.60217733e-19
#define AMTOKG 1.6605402e-27
#define toTHz 15.63330230023 /* sqrt(EVTOJ/1e-20/AMTOKG)*1.e-12/PI2 */
#define UtoGPa 1.60217733e2 /* EVTOJ*1.e30*1.e-9 */
#define KeVac 8.9875517873681764e9

/* converstion matrix for bcc, fcc phonon diperstions in cubic directions */

double bcc[9] = {-1.,1.,1.,1.,-1.,1.,1.,1.,-1};
double fcc[9] = { 0.,1.,1.,1., 0.,1.,1.,1.,0.};
double scc[9] = { 1.,0.,0.,0., 1.,0.,0.,0.,1.};
double hcp[9] = {1., 0.,0.,1.,-1.,0.,0.,0.,1.};
//double rho[9] = {1.,-1.,0.,0.,1.,-1.,1.,1.,1.};
double tet2[9] = {-1.,1.,1.,1.,-1.,1.,1.,1.,-1};

/* a memory storing the direction of electric field for IR LO-TO splitting */

double dirE[3] = {0.,0.,0.};

/* struct for atomic mass infor */

typedef struct {
    int Number; //Atomic number
    char Element[4]; //elemental symbol
    char ElementName[20]; //name of the element
    double mass; //mass 0f elememt
} AMASS;

//.wxdragon/prefs> tail +4 dragon_atoms_pref | awk '{printf "    \{\"%s\",%f\},\n", $1,$7}'
//awk '{print "\{"$1","" \""$2"\", \""$3"\", "$4"\},"}' periodictable

AMASS amass_avg[] = {
{1, "H", "Hydrogen", 1.00794},
{2, "He", "Helium", 4.002602},
{3, "Li", "Lithium", 6.941},
{4, "Be", "Beryllium", 9.012182},
{5, "B", "Boron", 10.811},
{6, "C", "Carbon", 12.0107},
{7, "N", "Nitrogen", 14.0067},
{8, "O", "Oxygen", 15.9994},
{9, "F", "Fluorine", 18.9984032},
{10, "Ne", "Neon", 20.1797},
{11, "Na", "Sodium", 22.98976928},
{12, "Mg", "Magnesium", 24.3050},
{13, "Al", "Aluminium", 26.9815386},
{14, "Si", "Silicon", 28.0855},
{15, "P", "Phosphorus", 30.973762},
{16, "S", "Sulfur", 32.065},
{17, "Cl", "Chlorine", 35.453},
{18, "Ar", "Argon", 39.948},
{19, "K", "Potassium", 39.0983},
{20, "Ca", "Calcium", 40.078},
{21, "Sc", "Scandium", 44.955912},
{22, "Ti", "Titanium", 47.867},
{23, "V", "Vanadium", 50.9415},
{24, "Cr", "Chromium", 51.9961},
{25, "Mn", "Manganese", 54.938045},
{26, "Fe", "Iron", 55.845},
{27, "Co", "Cobalt", 58.933195},
{28, "Ni", "Nickel", 58.6934},
{29, "Cu", "Copper", 63.546},
{30, "Zn", "Zinc", 65.38},
{31, "Ga", "Gallium", 69.723},
{32, "Ge", "Germanium", 72.64},
{33, "As", "Arsenic", 74.92160},
{34, "Se", "Selenium", 78.96},
{35, "Br", "Bromine", 79.904},
{36, "Kr", "Krypton", 83.798},
{37, "Rb", "Rubidium", 85.4678},
{38, "Sr", "Strontium", 87.62},
{39, "Y", "Yttrium", 88.90585},
{40, "Zr", "Zirconium", 91.224},
{41, "Nb", "Niobium", 92.90638},
{42, "Mo", "Molybdenum", 95.96},
{43, "Tc", "Technetium", 98},
{44, "Ru", "Ruthenium", 101.07},
{45, "Rh", "Rhodium", 102.90550},
{46, "Pd", "Palladium", 106.42},
{47, "Ag", "Silver", 107.8682},
{48, "Cd", "Cadmium", 112.411},
{49, "In", "Indium", 114.818},
{50, "Sn", "Tin", 118.710},
{51, "Sb", "Antimony", 121.760},
{52, "Te", "Tellurium", 127.60},
{53, "I", "Iodine", 126.90447},
{54, "Xe", "Xenon", 131.293},
{55, "Cs", "Caesium", 132.9054519},
{56, "Ba", "Barium", 137.327},
{57, "La", "Lanthanum", 138.90547},
{58, "Ce", "Cerium", 140.116},
{59, "Pr", "Praseodymium", 140.90765},
{60, "Nd", "Neodymium", 144.242},
{61, "Pm", "Promethium", 145},
{62, "Sm", "Samarium", 150.36},
{63, "Eu", "Europium", 151.964},
{64, "Gd", "Gadolinium", 157.25},
{65, "Tb", "Terbium", 158.92535},
{66, "Dy", "Dysprosium", 162.500},
{67, "Ho", "Holmium", 164.93032},
{68, "Er", "Erbium", 167.259},
{69, "Tm", "Thulium", 168.93421},
{70, "Yb", "Ytterbium", 173.054},
{71, "Lu", "Lutetium", 174.9668},
{72, "Hf", "Hafnium", 178.49},
{73, "Ta", "Tantalum", 180.94788},
{74, "W", "Tungsten", 183.84},
{75, "Re", "Rhenium", 186.207},
{76, "Os", "Osmium", 190.23},
{77, "Ir", "Iridium", 192.217},
{78, "Pt", "Platinum", 195.084},
{79, "Au", "Gold", 196.966569},
{80, "Hg", "Mercury", 200.59},
{81, "Tl", "Thallium", 204.3833},
{82, "Pb", "Lead", 207.2},
{83, "Bi", "Bismuth", 208.98040},
{84, "Po", "Polonium", 210},
{85, "At", "Astatine", 210},
{86, "Rn", "Radon", 222},
{87, "Fr", "Francium", 223},
{88, "Ra", "Radium", 223},
{89, "Ac", "Actinium", 227},
{90, "Th", "Thorium", 232.03806},
{91, "Pa", "Protactinium", 231.03588},
{92, "U", "Uranium", 238.02891},
{93, "Np", "Neptunium", 237},
{94, "Pu", "Plutonium", 244},
{95, "Am", "Americium", 243},
{96, "Cm", "Curium", 247},
{97, "Bk", "Berkelium", 247},
{98, "Cf", "Californium", 251},
{99, "Es", "Einsteinium", 252},
{100, "Fm", "Fermium", 257},
{101, "Md", "Mendelevium", 258},
{102, "No", "Nobelium", 259},
{103, "Lr", "Lawrencium", 262},
{104, "Rf", "Rutherfordium", 261},
{105, "Db", "Dubnium", 262},
{106, "Sg", "Seaborgium", 266},
{107, "Bh", "Bohrium", 264},
{108, "Hs", "Hassium", 267},
{109, "Mt", "Meitnerium", 268},
{110, "Ds", "Darmstadtium", 271},
{111, "Rg", "Roentgenium", 272},
{112, "Cn", "Copernicium", 285},
{113, "Uut", "Ununtrium", 284},
{114, "Uuq", "Ununquadium", 289},
{115, "Uup", "Ununpentium", 288},
{116, "Uuh", "Ununhexium", 292},
{117, "Uus", "Ununseptium", 295},
{118, "Uuo", "Ununoctium", 294}
}; 
/*
    {"X",0.000000},
    {"H",1.008000},
    {"He",4.003000},
    {"Li",6.939000},
    {"Be",9.012000},
    {"B",10.810000},
    {"C",12.001000},
    {"N",14.007000},
    {"O",15.999400},
    {"F",19.000000},
    {"Ne",20.183000},
    {"Na",22.989800},
    {"Mg",24.312000},
    {"Al",26.980000},
    {"Si",28.090000},
    {"P",30.974000},
    {"S",32.064000},
    {"Cl",35.045300},
    {"Ar",39.948000},
    {"K",39.102000},
    {"Ca",40.080000},
    {"Sc",44.960000},
    {"Ti",47.900000},
    {"V",50.940000},
    {"Cr",52.000000},
    {"Mn",54.940000},
    {"Fe",55.850000},
    {"Co",58.930000},
    {"Ni",58.710000},
    {"Cu",63.540000},
    {"Zn",65.370000},
    {"Ga",69.720000},
    {"Ge",72.590000},
    {"As",74.920000},
    {"Se",78.960000},
    {"Br",79.909000},
    {"Kr",83.800000},
    {"Rb",85.470000},
    {"Sr",87.620000},
    {"Y",88.910000},
    {"Zr",91.220000},
    {"Nb",92.910000},
    {"Mo",95.940000},
    {"Tc",99.000000},
    {"Ru",101.100000},
    {"Rh",102.910000},
    {"Pd",106.400000},
    {"Ag",107.870000},
    {"Cd",112.400000},
    {"In",114.820000},
    {"Sn",118.690000},
    {"Sb",121.750000},
    {"Te",127.600000},
    {"I",126.900000},
    {"Xe",131.300000},
    {"Cs",132.910000},
    {"Ba",137.360000},
    {"La",138.910000},
    {"Ce",140.120000},
    {"Pr",140.910000},
    {"Nd",144.240000},
    {"Pm",147.000000},
    {"Sm",150.350000},
    {"Eu",151.960000},
    {"Gd",157.250000},
    {"Tb",158.920000},
    {"Dy",162.500000},
    {"Ho",164.930000},
    {"Er",167.260000},
    {"Tm",168.930000},
    {"Yb",173.040000},
    {"Lu",174.970000},
    {"Hf",178.490000},
    {"Ta",180.950000},
    {"W",183.850000},
    {"Re",186.230000},
    {"Os",190.200000},
    {"Ir",192.200000},
    {"Pt",195.090000},
    {"Au",196.970000},
    {"Hg",200.590000},
    {"Tl",204.370000},
    {"Pb",207.190000},
    {"Bi",208.980000},
    {"Po",210.000000},
    {"At",210.000000},
    {"Rn",222.000000},
    {"Fr",223.000000},
    {"Ra",226.000000},
    {"Ac",227.000000},
    {"Th",232.040000},
    {"Pa",231.000000},
    {"U",238.030000},
    {"Np",237.000000},
    {"Pu",242.000000},
    {"Am",243.000000},
    {"Cm",247.000000},
    {"Bk",249.000000},
    {"Cf",251.000000},
    {"Es",254.000000},
    {"Fm",253.000000},
    {"Md",256.000000},
    {"No",253.000000},
    {"Lr",257.000000},
    {"Db",259.000000},
    {"Jl",200.000000},
    {"Rf",200.000000},
    {"Bh",200.000000},
    {"Hn",200.000000},
    {"Mt",200.000000}
*/

AMASS *amass;
int nElement;

/* module setting up the atomic mass 
  doing so, amass canbe realloced for the purpose of adding more atoms */
void setamass()
{
    nElement = sizeof(amass_avg)/sizeof(AMASS);
    amass = (AMASS *) malloc( (size_t) (nElement*sizeof(AMASS )) );
    for (int i=0; i<nElement; i++) {
	    strcpy(amass[i].Element, amass_avg[i].Element);
//printf ("Ele %d = %s\n", i, amass[i].Element);
	    amass[i].Number = amass_avg[i].Number;
	    amass[i].mass = amass_avg[i].mass;
	    strcpy(amass[i].ElementName, amass_avg[i].ElementName);
    }
}

/* reset/add atomic mass for element e */
void setamass(char e[4], double mass)
{
    //int nElement = sizeof(amass);
    //printf ("nElement=%d, Check Element %s with mass %.6lf\n", nElement, e, mass);
    for (int i=0; i<nElement; i++) {
    //printf ("Check Element %s with %d with %s\n", e, i, amass[i].Element);
	if (!strcmp(amass[i].Element, e)) {
	    amass[i].mass = mass;
	    printf ("mass of Element %s is changed to %.6lf\n", amass[i].Element, mass);
	    return;
	}
    }
    amass = (AMASS *) realloc(amass, (size_t) ((nElement+1)*sizeof(AMASS )) );
    strncpy(amass[nElement].Element, e, 4);
    amass[nElement].mass = mass;
    printf ("add new Element %s with mass %.6lf\n", amass[nElement].Element, mass);
    nElement++;
}

/* reset/add atomic masses for elements read in from file fmass 
   this for particular use of SQS calculation */

void addamass(char *fmass) {
    char e[256];
    double mass;
    FILE *fp = fopen(fmass, "r");
    while (1) {
	foo = fscanf(fp, "%s %lf", e, &mass);
	if (feof(fp)) break;
	setamass(e, mass);
    }
    fclose(fp);
}

/* return atomic mass by elemental symbol e */

double getmass(char e[4]){
	//for (int i=0; i<sizeof(amass)/sizeof(AMASS); i++) 
	for (int i=0; i<nElement; i++) 
		if (!strcmp(amass[i].Element, e)) return amass[i].mass;
	printf("\n********FETAL ERROR! CANNOT FIND MASS of element \"%s\"\n\n", e);
	exit(1);
}

/* return atomic number by elemental symbol e */

int getElementNumber(char e[4]){
        //for (int i=0; i<sizeof(amass)/sizeof(AMASS); i++)
        for (int i=0; i<nElement; i++) 
                if (!strcmp(amass[i].Element, e)) {
			//printf ("EleNum %d %s %d\n", i, amass[i].Element, amass[i].Number);
			return amass[i].Number;
		}
        printf("\n********FETAL ERROR! CANNOT FIND MASS of element \"%s\"\n\n", e);
        exit(1);
}

/* struct for atomic Neutron scatering cross sections
  from V. F. Sears, Neutron News 3, 26 (1992).
  http://www.ncnr.nist.gov/resources/n-lengths/list.html 
  Yphon only uses the column Elment and Scattxs */

typedef struct {
    char Element[20]; //Yphon uses this field
    char conc[20]; //not used
    char Cohb[20]; //not used
    char incb[20]; //not used
    char Cohxs[20]; //not used
    char Incxs[20]; //not used
    char Scattxs[20]; //Yphon uses this field
    char Absxs[20]; //not used
} XCROSS;

XCROSS xcross[] = {
  {"H", "---", "-3.7390", "---", "1.7568", "80.26", "82.02", "0.3326"}, 
  {"1H", "99.985", "-3.7406", "25.274", "1.7583", "80.27", "82.03", "0.3326"}, 
  {"D", "0.015", "6.671", "4.04", "5.592", "2.05", "7.64", "0.000519"}, 
  {"3H", "(12.32a)", "4.792", "-1.04", "2.89", "0.14", "3.03", "0"}, 
  {"He", "---", "3.26(3)", "---", "1.34", "0", "1.34", "0.00747"}, 
  {"3He", "0.00014", "5.74-1.483i", "-2.5+2.568i", "4.42", "1.6", "6", "5333.(7.)"}, 
  {"4He", "99.99986", "3.26", "0", "1.34", "0", "1.34", "0"}, 
  {"Li", "---", "-1.90", "---", "0.454", "0.92", "1.37", "70.5"}, 
  {"6Li", "7.5", "2.00-0.261i", "-1.89+0.26i", "0.51", "0.46", "0.97", "940.(4.)"}, 
  {"7Li", "92.5", "-2.22", "-2.49", "0.619", "0.78", "1.4", "0.0454"}, 
  {"Be", "100", "7.79", "0.12", "7.63", "0.0018", "7.63", "0.0076"}, 
  {"B", "---", "5.30-0.213i", "---", "3.54", "1.7", "5.24", "767.(8.)"}, 
  {"10B", "20", "-0.1-1.066i", "-4.7+1.231i", "0.144", "3", "3.1", "3835.(9.)"}, 
  {"11B", "80", "6.65", "-1.3", "5.56", "0.21", "5.77", "0.0055"}, 
  {"C", "---", "6.6460", "---", "5.551", "0.001", "5.551", "0.0035"}, 
  {"12C", "98.9", "6.6511", "0", "5.559", "0", "5.559", "0.00353"}, 
  {"13C", "1.1", "6.19", "-0.52", "4.81", "0.034", "4.84", "0.00137"}, 
  {"N", "---", "9.36", "---", "11.01", "0.5", "11.51", "1.9"}, 
  {"14N", "99.63", "9.37", "2.0", "11.03", "0.5", "11.53", "1.91"}, 
  {"15N", "0.37", "6.44", "-0.02", "5.21", "0.00005", "5.21", "0.000024"}, 
  {"O", "---", "5.803", "---", "4.232", "0.0008", "4.232", "0.00019"}, 
  {"16O", "99.762", "5.803", "0", "4.232", "0", "4.232", "0.0001"}, 
  {"17O", "0.038", "5.78", "0.18", "4.2", "0.004", "4.2", "0.236"}, 
  {"18O", "0.2", "5.84", "0", "4.29", "0", "4.29", "0.00016"}, 
  {"F", "100", "5.654", "-0.082", "4.017", "0.0008", "4.018", "0.0096"}, 
  {"Ne", "---", "4.566", "---", "2.62", "0.008", "2.628", "0.039"}, 
  {"20Ne", "90.51", "4.631", "0", "2.695", "0", "2.695", "0.036"}, 
  {"21Ne", "0.27", "6.66", "(+/-)0.6", "5.6", "0.05", "5.7", "0.67"}, 
  {"22Ne", "9.22", "3.87", "0", "1.88", "0", "1.88", "0.046"}, 
  {"Na", "100", "3.63", "3.59", "1.66", "1.62", "3.28", "0.53"}, 
  {"Mg", "---", "5.375", "---", "3.631", "0.08", "3.71", "0.063"}, 
  {"24Mg", "78.99", "5.66", "0", "4.03", "0", "4.03", "0.05"}, 
  {"25Mg", "10", "3.62", "1.48", "1.65", "0.28", "1.93", "0.19"}, 
  {"26Mg", "11.01", "4.89", "0", "3", "0", "3", "0.0382"}, 
  {"Al", "100", "3.449", "0.256", "1.495", "0.0082", "1.503", "0.231"}, 
  {"Si", "---", "4.1491", "---", "2.163", "0.004", "2.167", "0.171"}, 
  {"28Si", "92.23", "4.107", "0", "2.12", "0", "2.12", "0.177"}, 
  {"29Si", "4.67", "4.70", "0.09", "2.78", "0.001", "2.78", "0.101"}, 
  {"30Si", "3.1", "4.58", "0", "2.64", "0", "2.64", "0.107"}, 
  {"P", "100", "5.13", "0.2", "3.307", "0.005", "3.312", "0.172"}, 
  {"S", "---", "2.847", "---", "1.0186", "0.007", "1.026", "0.53"}, 
  {"32S", "95.02", "2.804", "0", "0.988", "0", "0.988", "0.54"}, 
  {"33S", "0.75", "4.74", "1.5", "2.8", "0.3", "3.1", "0.54"}, 
  {"34S", "4.21", "3.48", "0", "1.52", "0", "1.52", "0.227"}, 
  {"36S", "0.02", "3.(1.)", "0", "1.1", "0", "1.1", "0.15"}, 
  {"Cl", "---", "9.5770", "---", "11.5257", "5.3", "16.8", "33.5"}, 
  {"35Cl", "75.77", "11.65", "6.1", "17.06", "4.7", "21.8", "44.1"}, 
  {"37Cl", "24.23", "3.08", "0.1", "1.19", "0.001", "1.19", "0.433"}, 
  {"Ar", "---", "1.909", "---", "0.458", "0.225", "0.683", "0.675"}, 
  {"36Ar", "0.337", "24.90", "0", "77.9", "0", "77.9", "5.2"}, 
  {"38Ar", "0.063", "3.5", "0", "1.5(3.1)", "0", "1.5(3.1)", "0.8"}, 
  {"40Ar", "99.6", "1.830", "0", "0.421", "0", "0.421", "0.66"}, 
  {"K", "---", "3.67", "---", "1.69", "0.27", "1.96", "2.1"}, 
  {"39K", "93.258", "3.74", "1.4", "1.76", "0.25", "2.01", "2.1"}, 
  {"40K", "0.012", "3.(1.)", "---", "1.1", "0.5", "1.6", "35.(8.)"}, 
  {"41K", "6.73", "2.69", "1.5", "0.91", "0.3", "1.2", "1.46"}, 
  {"Ca", "---", "4.70", "---", "2.78", "0.05", "2.83", "0.43"}, 
  {"40Ca", "96.941", "4.80", "0", "2.9", "0", "2.9", "0.41"}, 
  {"42Ca", "0.647", "3.36", "0", "1.42", "0", "1.42", "0.68"}, 
  {"43Ca", "0.135", "-1.56", "---", "0.31", "0.5", "0.8", "6.2"}, 
  {"44Ca", "2.086", "1.42", "0", "0.25", "0", "0.25", "0.88"}, 
  {"46Ca", "0.004", "3.6", "0", "1.6", "0", "1.6", "0.74"}, 
  {"48Ca", "0.187", "0.39", "0", "0.019", "0", "0.019", "1.09"}, 
  {"Sc", "100", "12.29", "-6.0", "19", "4.5", "23.5", "27.5"}, 
  {"Ti", "---", "-3.438", "---", "1.485", "2.87", "4.35", "6.09"}, 
  {"46Ti", "8.2", "4.93", "0", "3.05", "0", "3.05", "0.59"}, 
  {"47Ti", "7.4", "3.63", "-3.5", "1.66", "1.5", "3.2", "1.7"}, 
  {"48Ti", "73.8", "-6.08", "0", "4.65", "0", "4.65", "7.84"}, 
  {"49Ti", "5.4", "1.04", "5.1", "0.14", "3.3", "3.4", "2.2"}, 
  {"50Ti", "5.2", "6.18", "0", "4.8", "0", "4.8", "0.179"}, 
  {"V", "---", "-0.3824", "---", "0.0184", "5.08", "5.1", "5.08"}, 
  {"50V", "0.25", "7.6", "---", "7.3(1.1)", "0.5", "7.8(1.0)", "60.(40.)"}, 
  {"51V", "99.75", "-0.402", "6.35", "0.0203", "5.07", "5.09", "4.9"}, 
  {"Cr", "---", "3.635", "---", "1.66", "1.83", "3.49", "3.05"}, 
  {"50Cr", "4.35", "-4.50", "0", "2.54", "0", "2.54", "15.8"}, 
  {"52Cr", "83.79", "4.920", "0", "3.042", "0", "3.042", "0.76"}, 
  {"53Cr", "9.5", "-4.20", "6.87", "2.22", "5.93", "8.15", "18.1(1.5)"}, 
  {"54Cr", "2.36", "4.55", "0", "2.6", "0", "2.6", "0.36"}, 
  {"Mn", "100", "-3.73", "1.79", "1.75", "0.4", "2.15", "13.3"}, 
  {"Fe", "---", "9.45", "---", "11.22", "0.4", "11.62", "2.56"}, 
  {"54Fe", "5.8", "4.2", "0", "2.2", "0", "2.2", "2.25"}, 
  {"56Fe", "91.7", "9.94", "0", "12.42", "0", "12.42", "2.59"}, 
  {"57Fe", "2.2", "2.3", "---", "0.66", "0.3", "1", "2.48"}, 
  {"58Fe", "0.3", "15.(7.)", "0", "28", "0", "28.(26.)", "1.28"}, 
  {"Co", "100", "2.49", "-6.2", "0.779", "4.8", "5.6", "37.18"}, 
  {"Ni", "---", "10.3", "---", "13.3", "5.2", "18.5", "4.49"}, 
  {"58Ni", "68.27", "14.4", "0", "26.1", "0", "26.1", "4.6"}, 
  {"60Ni", "26.1", "2.8", "0", "0.99", "0", "0.99", "2.9"}, 
  {"61Ni", "1.13", "7.60", "(+/-)3.9", "7.26", "1.9", "9.2", "2.5"}, 
  {"62Ni", "3.59", "-8.7", "0", "9.5", "0", "9.5", "14.5"}, 
  {"64Ni", "0.91", "-0.37", "0", "0.017", "0", "0.017", "1.52"}, 
  {"Cu", "---", "7.718", "---", "7.485", "0.55", "8.03", "3.78"}, 
  {"63Cu", "69.17", "6.43", "0.22", "5.2", "0.006", "5.2", "4.5"}, 
  {"65Cu", "30.83", "10.61", "1.79", "14.1", "0.4", "14.5", "2.17"}, 
  {"Zn", "---", "5.680", "---", "4.054", "0.077", "4.131", "1.11"}, 
  {"64Zn", "48.6", "5.22", "0", "3.42", "0", "3.42", "0.93"}, 
  {"66Zn", "27.9", "5.97", "0", "4.48", "0", "4.48", "0.62"}, 
  {"67Zn", "4.1", "7.56", "-1.50", "7.18", "0.28", "7.46", "6.8"}, 
  {"68Zn", "18.8", "6.03", "0", "4.57", "0", "4.57", "1.1"}, 
  {"70Zn", "0.6", "6.(1.)", "0", "4.5", "0", "4.5(1.5)", "0.092"}, 
  {"Ga", "---", "7.288", "---", "6.675", "0.16", "6.83", "2.75"}, 
  {"69Ga", "60.1", "7.88", "-0.85", "7.8", "0.091", "7.89", "2.18"}, 
  {"71Ga", "39.9", "6.40", "-0.82", "5.15", "0.084", "5.23", "3.61"}, 
  {"Ge", "---", "8.185", "---", "8.42", "0.18", "8.6", "2.2"}, 
  {"70Ge", "20.5", "10.0", "0", "12.6", "0", "12.6", "3"}, 
  {"72Ge", "27.4", "8.51", "0", "9.1", "0", "9.1", "0.8"}, 
  {"73Ge", "7.8", "5.02", "3.4", "3.17", "1.5", "4.7", "15.1"}, 
  {"74Ge", "36.5", "7.58", "0", "7.2", "0", "7.2", "0.4"}, 
  {"76Ge", "7.8", "8.2", "0", "8.(3.)", "0", "8.(3.)", "0.16"}, 
  {"As", "100", "6.58", "-0.69", "5.44", "0.06", "5.5", "4.5"}, 
  {"Se", "---", "7.970", "---", "7.98", "0.32", "8.3", "11.7"}, 
  {"74Se", "0.9", "0.8", "0", "0.1", "0", "0.1", "51.8(1.2)"}, 
  {"76Se", "9", "12.2", "0", "18.7", "0", "18.7", "85.(7.)"}, 
  {"77Se", "7.6", "8.25", "(+/-)0.6(1.6)", "8.6", "0.05", "8.65", "42.(4.)"}, 
  {"78Se", "23.5", "8.24", "0", "8.5", "0", "8.5", "0.43"}, 
  {"80Se", "49.6", "7.48", "0", "7.03", "0", "7.03", "0.61"}, 
  {"82Se", "9.4", "6.34", "0", "5.05", "0", "5.05", "0.044"}, 
  {"Br", "---", "6.795", "---", "5.8", "0.1", "5.9", "6.9"}, 
  {"79Br", "50.69", "6.80", "-1.1", "5.81", "0.15", "5.96", "11"}, 
  {"81Br", "49.31", "6.79", "0.6", "5.79", "0.05", "5.84", "2.7"}, 
  {"Kr", "---", "7.81", "---", "7.67", "0.01", "7.68", "25.(1.)"}, 
  {"78Kr", "0.35", "---", "0", "---", "0", "---", "6.4"}, 
  {"80Kr", "2.25", "---", "0", "---", "0", "---", "11.8"}, 
  {"82Kr", "11.6", "---", "0", "---", "0", "---", "29.(20.)"}, 
  {"83Kr", "11.5", "---", "---", "---", "---", "---", "185.(30.)"}, 
  {"84Kr", "57", "---", "0", "---", "0", "6.6", "0.113"}, 
  {"86Kr", "17.3", "8.1", "0", "8.2", "0", "8.2", "0.003"}, 
  {"Rb", "---", "7.09", "---", "6.32", "0.5", "6.8", "0.38"}, 
  {"85Rb", "72.17", "7.03", "---", "6.2", "0.5", "6.7", "0.48"}, 
  {"87Rb", "27.83", "7.23", "---", "6.6", "0.5", "7.1", "0.12"}, 
  {"Sr", "---", "7.02", "---", "6.19", "0.06", "6.25", "1.28"}, 
  {"84Sr", "0.56", "7.(1.)", "0", "6.(2.)", "0", "6.(2.)", "0.87"}, 
  {"86Sr", "9.86", "5.67", "0", "4.04", "0", "4.04", "1.04"}, 
  {"87Sr", "7", "7.40", "---", "6.88", "0.5", "7.4", "16.(3.)"}, 
  {"88Sr", "82.58", "7.15", "0", "6.42", "0", "6.42", "0.058"}, 
  {"Y", "100", "7.75", "1.1", "7.55", "0.15", "7.7", "1.28"}, 
  {"Zr", "---", "7.16", "---", "6.44", "0.02", "6.46", "0.185"}, 
  {"90Zr", "51.45", "6.4", "0", "5.1", "0", "5.1", "0.011"}, 
  {"91Zr", "11.32", "8.7", "-1.08", "9.5", "0.15", "9.7", "1.17"}, 
  {"92Zr", "17.19", "7.4", "0", "6.9", "0", "6.9", "0.22"}, 
  {"94Zr", "17.28", "8.2", "0", "8.4", "0", "8.4", "0.0499"}, 
  {"96Zr", "2.76", "5.5", "0", "3.8", "0", "3.8", "0.0229"}, 
  {"Nb", "100", "7.054", "-0.139", "6.253", "0.0024", "6.255", "1.15"}, 
  {"Mo", "---", "6.715", "---", "5.67", "0.04", "5.71", "2.48"}, 
  {"92Mo", "14.84", "6.91", "0", "6", "0", "6", "0.019"}, 
  {"94Mo", "9.25", "6.80", "0", "5.81", "0", "5.81", "0.015"}, 
  {"95Mo", "15.92", "6.91", "---", "6", "0.5", "6.5", "13.1"}, 
  {"96Mo", "16.68", "6.20", "0", "4.83", "0", "4.83", "0.5"}, 
  {"97Mo", "9.55", "7.24", "---", "6.59", "0.5", "7.1", "2.5"}, 
  {"98Mo", "24.13", "6.58", "0", "5.44", "0", "5.44", "0.127"}, 
  {"100Mo", "9.63", "6.73", "0", "5.69", "0", "5.69", "0.4"}, 
  {"Tc", "(2.3E5a)", "6.8", "---", "5.8", "0.5", "6.3", "20.(1.)"}, 
  {"Ru", "---", "7.03", "---", "6.21", "0.4", "6.6", "2.56"}, 
  {"96Ru", "5.5", "---", "0", "---", "0", "---", "0.28"}, 
  {"98Ru", "1.9", "---", "0", "---", "0", "---", "<8."}, 
  {"99Ru", "12.7", "---", "---", "---", "---", "---", "6.9(1.0)"}, 
  {"100Ru", "12.6", "---", "0", "---", "0", "---", "4.8"}, 
  {"101Ru", "17", "---", "---", "---", "---", "---", "3.3"}, 
  {"102Ru", "31.6", "---", "0", "---", "0", "144.8", "1.17"}, 
  {"104Ru", "18.7", "---", "0", "---", "0", "4.483", "0.31"}, 
  {"Rh", "100", "5.88", "---", "4.34", "0.3", "4.6", "144.8"}, 
  {"Pd", "---", "5.91", "---", "4.39", "0.093", "4.48", "6.9"}, 
  {"102Pd", "1.02", "7.7(7)", "0", "7.5(1.4)", "0", "7.5(1.4)", "3.4"}, 
  {"104Pd", "11.14", "7.7(7)", "0", "7.5(1.4)", "0", "7.5(1.4)", "0.6"}, 
  {"105Pd", "22.33", "5.5", "-2.6(1.6)", "3.8", "0.8", "4.6(1.1)", "20.(3.)"}, 
  {"106Pd", "27.33", "6.4", "0", "5.1", "0", "5.1", "0.304"}, 
  {"108Pd", "26.46", "4.1", "0", "2.1", "0", "2.1", "8.55"}, 
  {"110Pd", "11.72", "7.7(7)", "0", "7.5(1.4)", "0", "7.5(1.4)", "0.226"}, 
  {"Ag", "---", "5.922", "---", "4.407", "0.58", "4.99", "63.3"}, 
  {"107Ag", "51.83", "7.555", "1.00", "7.17", "0.13", "7.3", "37.6(1.2)"}, 
  {"109Ag", "48.17", "4.165", "-1.60", "2.18", "0.32", "2.5", "91.0(1.0)"}, 
  {"Cd", "---", "4.87-0.70i", "---", "3.04", "3.46", "6.5", "2520.(50.)"}, 
  {"106Cd", "1.25", "5.(2.)", "0", "3.1", "0", "3.1(2.5)", "1"}, 
  {"108Cd", "0.89", "5.4", "0", "3.7", "0", "3.7", "1.1"}, 
  {"110Cd", "12.51", "5.9", "0", "4.4", "0", "4.4", "11"}, 
  {"111Cd", "12.81", "6.5", "---", "5.3", "0.3", "5.6", "24"}, 
  {"112Cd", "24.13", "6.4", "0", "5.1", "0", "5.1", "2.2"}, 
  {"113Cd", "12.22", "-8.0-5.73i", "---", "12.1", "0.3", "12.4", "20600.(400.)"}, 
  {"114Cd", "28.72", "7.5", "0", "7.1", "0", "7.1", "0.34"}, 
  {"116Cd", "7.47", "6.3", "0", "5", "0", "5", "0.075"}, 
  {"In", "---", "4.065-0.0539i", "---", "2.08", "0.54", "2.62", "193.8(1.5)"}, 
  {"113In", "4.3", "5.39", "(+/-)0.017", "3.65", "0.000037", "3.65", "12.0(1.1)"}, 
  {"115In", "95.7", "4.01-0.0562i", "-2.1", "2.02", "0.55", "2.57", "202.(2.)"}, 
  {"Sn", "---", "6.225", "---", "4.871", "0.022", "4.892", "0.626"}, 
  {"112Sn", "1", "6.(1.)", "0", "4.5(1.5)", "0", "4.5(1.5)", "1"}, 
  {"114Sn", "0.7", "6.2", "0", "4.8", "0", "4.8", "0.114"}, 
  {"115Sn", "0.4", "6.(1.)", "---", "4.5(1.5)", "0.3", "4.8(1.5)", "30.(7.)"}, 
  {"116Sn", "14.7", "5.93", "0", "4.42", "0", "4.42", "0.14"}, 
  {"117Sn", "7.7", "6.48", "---", "5.28", "0.3", "5.6", "2.3"}, 
  {"118Sn", "24.3", "6.07", "0", "4.63", "0", "4.63", "0.22"}, 
  {"119Sn", "8.6", "6.12", "---", "4.71", "0.3", "5", "2.2"}, 
  {"120Sn", "32.4", "6.49", "0", "5.29", "0", "5.29", "0.14"}, 
  {"122Sn", "4.6", "5.74", "0", "4.14", "0", "4.14", "0.18"}, 
  {"124Sn", "5.6", "5.97", "0", "4.48", "0", "4.48", "0.133"}, 
  {"Sb", "---", "5.57", "---", "3.9", "0.007", "3.9", "4.91"}, 
  {"121Sb", "57.3", "5.71", "-0.05", "4.1", "0.0003", "4.1", "5.75"}, 
  {"123Sb", "42.7", "5.38", "-0.10", "3.64", "0.001", "3.64", "3.8"}, 
  {"Te", "---", "5.80", "---", "4.23", "0.09", "4.32", "4.7"}, 
  {"120Te", "0.096", "5.3", "0", "3.5", "0", "3.5", "2.3"}, 
  {"122Te", "2.6", "3.8", "0", "1.8", "0", "1.8", "3.4"}, 
  {"123Te", "0.908", "-0.05-0.116i", "-2.04", "0.002", "0.52", "0.52", "418.(30.)"}, 
  {"124Te", "4.816", "7.96", "0", "8", "0", "8", "6.8(1.3)"}, 
  {"125Te", "7.14", "5.02", "-0.26", "3.17", "0.008", "3.18", "1.55"}, 
  {"126Te", "18.95", "5.56", "0", "3.88", "0", "3.88", "1.04"}, 
  {"128Te", "31.69", "5.89", "0", "4.36", "0", "4.36", "0.215"}, 
  {"130Te", "33.8", "6.02", "0", "4.55", "0", "4.55", "0.29"}, 
  {"I", "100", "5.28", "1.58", "3.5", "0.31", "3.81", "6.15"}, 
  {"Xe", "---", "4.92", "3.04", "2.96", "0", "---", "23.9(1.2)"}, 
  {"124Xe", "0.1", "---", "0", "---", "0", "---", "165.(20.)"}, 
  {"126Xe", "0.09", "---", "0", "---", "0", "---", "3.5"}, 
  {"128Xe", "1.91", "---", "0", "---", "0", "---", "<8"}, 
  {"129Xe", "26.4", "---", "---", "---", "---", "---", "21.(5.)"}, 
  {"130Xe", "4.1", "---", "0", "---", "0", "---", "<26."}, 
  {"131Xe", "21.2", "---", "---", "---", "---", "---", "85.(10.)"}, 
  {"132Xe", "26.9", "---", "0", "---", "0", "---", "0.45"}, 
  {"134Xe", "10.4", "---", "0", "---", "0", "---", "0.265"}, 
  {"136Xe", "8.9", "---", "0", "---", "0", "---", "0.26"}, 
  {"Cs", "100", "5.42", "1.29", "3.69", "0.21", "3.9", "29.0(1.5)"}, 
  {"Ba", "---", "5.07", "---", "3.23", "0.15", "3.38", "1.1"}, 
  {"130Ba", "0.11", "-3.6", "0", "1.6", "0", "1.6", "30.(5.)"}, 
  {"132Ba", "0.1", "7.8", "0", "7.6", "0", "7.6", "7"}, 
  {"134Ba", "2.42", "5.7", "0", "4.08", "0", "4.08", "2.0(1.6)"}, 
  {"135Ba", "6.59", "4.67", "---", "2.74", "0.5", "3.2", "5.8"}, 
  {"136Ba", "7.85", "4.91", "0", "3.03", "0", "3.03", "0.68"}, 
  {"137Ba", "11.23", "6.83", "---", "5.86", "0.5", "6.4", "3.6"}, 
  {"138Ba", "71.7", "4.84", "0", "2.94", "0", "2.94", "0.27"}, 
  {"La", "---", "8.24", "---", "8.53", "1.13", "9.66", "8.97"}, 
  {"138La", "0.09", "8.(2.)", "---", "8.(4.)", "0.5", "8.5(4.0)", "57.(6.)"}, 
  {"139La", "99.91", "8.24", "3.0", "8.53", "1.13", "9.66", "8.93"}, 
  {"Ce", "---", "4.84", "---", "2.94", "0.001", "2.94", "0.63"}, 
  {"136Ce", "0.19", "5.80", "0", "4.23", "0", "4.23", "7.3(1.5)"}, 
  {"138Ce", "0.25", "6.70", "0", "5.64", "0", "5.64", "1.1"}, 
  {"140Ce", "88.48", "4.84", "0", "2.94", "0", "2.94", "0.57"}, 
  {"142Ce", "11.08", "4.75", "0", "2.84", "0", "2.84", "0.95"}, 
  {"Pr", "100", "4.58", "-0.35", "2.64", "0.015", "2.66", "11.5"}, 
  {"Nd", "---", "7.69", "---", "7.43", "9.2", "16.6", "50.5(1.2)"}, 
  {"142Nd", "27.16", "7.7", "0", "7.5", "0", "7.5", "18.7"}, 
  {"143Nd", "12.18", "14.(2.)", "(+/-)21.(1.)", "25.(7.)", "55.(7.)", "80.(2.)", "337.(10.)"}, 
  {"144Nd", "23.8", "2.8", "0", "1", "0", "1", "3.6"}, 
  {"145Nd", "8.29", "14.(2.)", "---", "25.(7.)", "5.(5.)", "30.(9.)", "42.(2.)"}, 
  {"146Nd", "17.19", "8.7", "0", "9.5", "0", "9.5", "1.4"}, 
  {"148Nd", "5.75", "5.7", "0", "4.1", "0", "4.1", "2.5"}, 
  {"150Nd", "5.63", "5.3", "0", "3.5", "0", "3.5", "1.2"}, 
  {"Pm", "(2.62a)", "12.6", "(+/-)3.2(2.5)", "20.0(1.3)", "1.3(2.0)", "21.3(1.5)", "168.4(3.5)"}, 
  {"Sm", "---", "0.80-1.65i", "---", "0.422", "39.(3.)", "39.(3.)", "5922.(56.)"}, 
  {"144Sm", "3.1", "-3.(4.)", "0", "1.(3.)", "0", "1.(3.)", "0.7"}, 
  {"147Sm", "15.1", "14.(3.)", "(+/-)11.(7.)", "25.(11.)", "143(19.)", "39.(16.)", "57.(3.)"}, 
  {"148Sm", "11.3", "-3.(4.)", "0", "1.(3.)", "0", "1.(3.)", "2.4"}, 
  {"149Sm", "13.9", "-19.2-11.7i", "(+/-)31.4-10.3i", "63.5", "137.(5.)", "200.(5.)", "42080.(400.)"}, 
  {"150Sm", "7.4", "14.(3.)", "0", "25.(11.)", "0", "25.(11.)", "104.(4.)"}, 
  {"152Sm", "26.6", "-5.0", "0", "3.1", "0", "3.1", "206.(6.)"}, 
  {"154Sm", "22.6", "9.3", "0", "11.(2.)", "0", "11.(2.)", "8.4"}, 
  {"Eu", "---", "7.22-1.26i", "---", "6.57", "2.5", "9.2", "4530.(40.)"}, 
  {"151Eu", "47.8", "6.13-2.53i", "(+/-)4.5-2.14i", "5.5", "3.1", "8.6", "9100.(100.)"}, 
  {"153Eu", "52.2", "8.22", "(+/-)3.2", "8.5", "1.3", "9.8", "312.(7.)"}, 
  {"Gd", "---", "6.5-13.82i", "---", "29.3", "151.(2.)", "180.(2.)", "49700.(125.)"}, 
  {"152Gd", "0.2", "10.(3.)", "0", "13.(8.)", "0", "13.(8.)", "735.(20.)"}, 
  {"154Gd", "2.1", "10.(3.)", "0", "13.(8.)", "0", "13.(8.)", "85.(12.)"}, 
  {"155Gd", "14.8", "6.0-17.0i", "(+/-)5.(5.)-13.16i", "40.8", "25.(6.)", "66.(6.)", "61100.(400.)"}, 
  {"156Gd", "20.6", "6.3", "0", "5", "0", "5", "1.5(1.2)"}, 
  {"157Gd", "15.7", "-1.14-71.9i", "(+/-)5.(5.)-55.8i", "650.(4.)", "394.(7.)", "1044.(8.)", "259000.(700.)"}, 
  {"158Gd", "24.8", "9.(2.)", "0", "10.(5.)", "0", "10.(5.)", "2.2"}, 
  {"160Gd", "21.8", "9.15", "0", "10.52", "0", "10.52", "0.77"}, 
  {"Tb", "100", "7.38", "-0.17", "6.84", "0.004", "6.84", "23.4"}, 
  {"Dy", "---", "16.9-0.276i", "---", "35.9", "54.4(1.2)", "90.3", "994.(13.)"}, 
  {"156Dy", "0.06", "6.1", "0", "4.7", "0", "4.7", "33.(3.)"}, 
  {"158Dy", "0.1", "6.(4.)", "0", "5.(6.)", "0", "5.(6.)", "43.(6.)"}, 
  {"160Dy", "2.34", "6.7", "0", "5.6", "0", "5.6", "56.(5.)"}, 
  {"161Dy", "19", "10.3", "(+/-)4.9", "13.3", "3.(1.)", "16.(1.)", "600.(25.)"}, 
  {"162Dy", "25.5", "-1.4", "0", "0.25", "0", "0.25", "194.(10.)"}, 
  {"163Dy", "24.9", "5.0", "1.3", "3.1", "0.21", "3.3", "124.(7.)"}, 
  {"164Dy", "28.1", "49.4-0.79i", "0", "307.(3.)", "0", "307.(3.)", "2840.(40.)"}, 
  {"Ho", "100", "8.01", "-1.70", "8.06", "0.36", "8.42", "64.7(1.2)"}, 
  {"Er", "---", "7.79", "---", "7.63", "1.1", "8.7", "159.(4.)"}, 
  {"162Er", "0.14", "8.8", "0", "9.7", "0", "9.7", "19.(2.)"}, 
  {"164Er", "1.56", "8.2", "0", "8.4", "0", "8.4", "13.(2.)"}, 
  {"166Er", "33.4", "10.6", "0", "14.1", "0", "14.1", "19.6(1.5)"}, 
  {"167Er", "22.9", "3.0", "1.0", "1.1", "0.13", "1.2", "659.(16.)"}, 
  {"168Er", "27.1", "7.4", "0", "6.9", "0", "6.9", "2.74"}, 
  {"170Er", "14.9", "9.6", "0", "11.6", "0", "11.6(1.2)", "5.8"}, 
  {"Tm", "100", "7.07", "0.9", "6.28", "0.1", "6.38", "100.(2.)"}, 
  {"Yb", "---", "12.43", "---", "19.42", "4", "23.4", "34.8"}, 
  {"168Yb", "0.14", "-4.07-0.62i", "0", "2.13", "0", "2.13", "2230.(40.)"}, 
  {"170Yb", "3.06", "6.77", "0", "5.8", "0", "5.8", "11.4(1.0)"}, 
  {"171Yb", "14.3", "9.66", "-5.59", "11.7", "3.9", "15.6", "48.6(2.5)"}, 
  {"172Yb", "21.9", "9.43", "0", "11.2", "0", "11.2", "0.8"}, 
  {"173Yb", "16.1", "9.56", "-5.3", "11.5", "3.5", "15", "17.1(1.3)"}, 
  {"174Yb", "31.8", "19.3", "0", "46.8", "0", "46.8", "69.4(5.0)"}, 
  {"176Yb", "12.7", "8.72", "0", "9.6", "0", "9.6", "2.85"}, 
  {"Lu", "---", "7.21", "---", "6.53", "0.7", "7.2", "74.(2.)"}, 
  {"175Lu", "97.39", "7.24", "(+/-)2.2", "6.59", "0.6", "7.2", "21.(3.)"}, 
  {"176Lu", "2.61", "6.1-0.57i", "(+/-)3.0+0.61i", "4.7", "1.2", "5.9", "2065.(35.)"}, 
  {"Hf", "---", "7.7", "---", "7.6", "2.6", "10.2", "104.1"}, 
  {"174Hf", "0.2", "10.9(1.1)", "0", "15.(3.)", "0", "15.(3.)", "561.(35.)"}, 
  {"176Hf", "5.2", "6.61", "0", "5.5", "0", "5.5", "23.5(3.1)"}, 
  {"177Hf", "18.6", "0.8(1.0)", "(+/-)0.9(1.3)", "0.1", "0.1", "0.2", "373.(10.)"}, 
  {"178Hf", "27.1", "5.9", "0", "4.4", "0", "4.4", "84.(4.)"}, 
  {"179Hf", "13.7", "7.46", "(+/-)1.06", "7", "0.14", "7.1", "41.(3.)"}, 
  {"180Hf", "35.2", "13.2", "0", "21.9", "0", "21.9(1.0)", "13.04"}, 
  {"Ta", "---", "6.91", "---", "6", "0.01", "6.01", "20.6"}, 
  {"180Ta", "0.012", "7.(2.)", "---", "6.2", "0.5", "7.(4.)", "563.(60.)"}, 
  {"181Ta", "99.988", "6.91", "-0.29", "6", "0.011", "6.01", "20.5"}, 
  {"W", "---", "4.86", "---", "2.97", "1.63", "4.6", "18.3"}, 
  {"180W", "0.1", "5.(3.)", "0", "3.(4.)", "0", "3.(4.)", "30.(20.)"}, 
  {"182W", "26.3", "6.97", "0", "6.1", "0", "6.1", "20.7"}, 
  {"183W", "14.3", "6.53", "---", "5.36", "0.3", "5.7", "10.1"}, 
  {"184W", "30.7", "7.48", "0", "7.03", "0", "7.03", "1.7"}, 
  {"186W", "28.6", "-0.72", "0", "0.065", "0", "0.065", "37.9"}, 
  {"Re", "---", "9.2", "---", "10.6", "0.9", "11.5", "89.7(1.)"}, 
  {"185Re", "37.4", "9.0", "(+/-)2.0", "10.2", "0.5", "10.7", "112.(2.)"}, 
  {"187Re", "62.6", "9.3", "(+/-)2.8", "10.9", "1", "11.9", "76.4(1.)"}, 
  {"Os", "---", "10.7", "---", "14.4", "0.3", "14.7", "16"}, 
  {"184Os", "0.02", "10.(2.)", "0", "13.(5.)", "0", "13.(5.)", "3000.(150.)"}, 
  {"186Os", "1.58", "11.6(1.7)", "0", "17.(5.)", "0", "17.(5.)", "80.(13.)"}, 
  {"187Os", "1.6", "10.(2.)", "---", "13.(5.)", "0.3", "13.(5.)", "320.(10.)"}, 
  {"188Os", "13.3", "7.6", "0", "7.3", "0", "7.3", "4.7"}, 
  {"189Os", "16.1", "10.7", "---", "14.4", "0.5", "14.9", "25.(4.)"}, 
  {"190Os", "26.4", "11.0", "0", "15.2", "0", "15.2", "13.1"}, 
  {"192Os", "41", "11.5", "0", "16.6", "0", "16.6(1.2)", "2"}, 
  {"Ir", "---", "10.6", "---", "14.1", "0.(3.)", "14.(3.)", "425.(2.)"}, 
  {"191Ir", "37.3", "---", "---", "---", "---", "---", "954.(10.)"}, 
  {"193Ir", "62.7", "---", "---", "---", "---", "---", "111.(5.)"}, 
  {"Pt", "---", "9.60", "---", "11.58", "0.13", "11.71", "10.3"}, 
  {"190Pt", "0.01", "9.0", "0", "10.(2.)", "0", "10.(2.)", "152.(4.)"}, 
  {"192Pt", "0.79", "9.9", "0", "12.3(1.2)", "0", "12.3(1.2)", "10.0(2.5)"}, 
  {"194Pt", "32.9", "10.55", "0", "14", "0", "14", "1.44"}, 
  {"195Pt", "33.8", "8.83", "-1.00", "9.8", "0.13", "9.9", "27.5(1.2)"}, 
  {"196Pt", "25.3", "9.89", "0", "12.3", "0", "12.3", "0.72"}, 
  {"198Pt", "7.2", "7.8", "0", "7.6", "0", "7.6", "3.66"}, 
  {"Au", "100", "7.63", "-1.84", "7.32", "0.43", "7.75", "98.65"}, 
  {"Hg", "---", "12.692", "---", "20.24", "6.6", "26.8", "372.3(4.0)"}, 
  {"196Hg", "0.2", "30.3(1.0)", "0", "115.(8.)", "0", "115.(8.)", "3080.(180.)"}, 
  {"198Hg", "10.1", "---", "0", "---", "0", "---", "2"}, 
  {"199Hg", "17", "16.9", "(+/-)15.5", "36.(2.)", "30.(3.)", "66.(2.)", "2150.(48.)"}, 
  {"200Hg", "23.1", "---", "0", "---", "0", "---", "<60."}, 
  {"201Hg", "13.2", "---", "---", "---", "---", "---", "7.8(2.0)"}, 
  {"202Hg", "29.6", "---", "0", "---", "0", "9.828", "4.89"}, 
  {"204Hg", "6.8", "---", "0", "---", "0", "---", "0.43"}, 
  {"Tl", "---", "8.776", "---", "9.678", "0.21", "9.89", "3.43"}, 
  {"203Tl", "29.524", "6.99", "1.06", "6.14", "0.14", "6.28", "11.4"}, 
  {"205Tl", "70.476", "9.52", "-0.242", "11.39", "0.007", "11.4", "0.104"}, 
  {"Pb", "---", "9.405", "---", "11.115", "0.003", "11.118", "0.171"}, 
  {"204Pb", "1.4", "9.90", "0", "12.3", "0", "12.3", "0.65"}, 
  {"206Pb", "24.1", "9.22", "0", "10.68", "0", "10.68", "0.03"}, 
  {"207Pb", "22.1", "9.28", "0.14", "10.82", "0.002", "10.82", "0.699"}, 
  {"208Pb", "52.4", "9.50", "0", "11.34", "0", "11.34", "0.00048"}, 
  {"Bi", "100", "8.532", "---", "9.148", "0.0084", "9.156", "0.0338"}, 
  {"Po", "---", "---", "0.259", "0", "---", "---", "---"}, 
  {"At", "---", "---", "---", "0", "---", "---", "---"}, 
  {"Rn", "---", "---", "---", "0", "---", "12.6", "---"}, 
  {"Fr", "---", "---", "---", "0", "---", "---", "---"}, 
  {"Ra", "(1.60E3a)", "10.0(1.0)", "0", "13.(3.)", "0", "13.(3.)", "12.8(1.5)"}, 
  {"Ac", "---", "---", "---", "0", "---", "---", "---"}, 
  {"Th", "100", "10.31", "0", "13.36", "0", "13.36", "7.37"}, 
  {"Pa", "(3.28E4a)", "9.1", "---", "10.4", "0.1(3.3)", "10.5(3.2)", "200.6(2.3)"}, 
  {"U", "---", "8.417", "---", "8.903", "0.005", "8.908", "7.57"}, 
  {"233U", "(1.59E5a)", "10.1", "(+/-)1.(3.)", "12.8", "0.1", "12.9", "574.7(1.0)"}, 
  {"234U", "0.005", "12.4", "0", "19.3", "0", "19.3", "100.1(1.3)"}, 
  {"235U", "0.72", "10.47", "(+/-)1.3", "13.78", "0.2", "14", "680.9(1.1)"}, 
  {"238U", "99.275", "8.402", "0", "8.871", "0", "8.871", "2.68"}, 
  {"Np", "(2.14E6a)", "10.55", "---", "14", "0.5", "14.5", "175.9(2.9)"}, 
  {"Pu", "---", "---", "---", "---", "---", "7.7", "---"}, 
  {"238Pu", "(87.74a)", "14.1", "0", "25.0(1.8)", "0", "25.0(1.8)", "558.(7.)"}, 
  {"239Pu", "(2.41E4a)", "7.7", "(+/-)1.3(1.9)", "7.5", "0.2", "7.7", "1017.3(2.1)"}, 
  {"240Pu", "(6.56E3a)", "3.5", "0", "1.54", "0", "1.54", "289.6(1.4)"}, 
  {"242Pu", "(3.76E5a)", "8.1", "0", "8.2", "0", "8.2", "18.5"}, 
  {"Am", "(7.37E3a)", "8.3", "(+/-)2.(7.)", "8.7", "0.3", "9.0(2.6)", "75.3(1.8)"}, 
  {"Cm", "---", "---", "---", "0", "---", "---", "---"}, 
  {"244Cm", "(18.10a)", "9.5", "0", "11.3", "0", "11.3", "16.2(1.2)"}, 
  {"246Cm", "(4.7E3a)", "9.3", "0", "10.9", "0", "10.9", "1.36"}, 
  {"248Cm", "(3.5E5a)", "7.7", "0", "7.5", "0", "7.5", "3"}
};

/* module to return the Neutron scattering cross sections by elemental symbol e */
double getSdm(char e[3]){
	int nEl = sizeof(xcross)/sizeof(XCROSS);
	for (int i=0; i<nEl; i++) 
	    if (!strcmp(xcross[i].Element, e)) {
		double mass = getmass(e);
		double Scattxs = atof(xcross[i].Scattxs);
		if (mass <= 0.0 || Scattxs <= 0.0) break;
		return Scattxs/mass;
		//return pow(Scattxs/mass,2);
	    }
	printf("\n********FETAL ERROR! CANNOT FIND CROSS SECTION of element \"%s\"\n\n", e);
	exit(1);
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

/* equivalent to Linux getline, for compatibility with Microsoft C++ */
 
int yetline(char **line, size_t *len, FILE *born)
{
	char tmp[MAX_STR_LEN];
	int i;
	for (i=0; ; i++) 
	{
		foo = fscanf(born, "%c", tmp+i);
		if (tmp[i]=='\n') break;
		if (feof(born)) 
		{
			tmp[i]='\n';
			break;
		}
	}
	if (i==0 && feof(born)) return -1;
	*len = i;
	tmp[i+1] = 0;
	*line = new char[i+2];
	strcpy (*line, tmp);
	return 0;
}

/* routine for returning Linux time
	option - Yphon command line options
	Nopt - number of command line options
*/

void datetimehost(char **option, int Nopt) {
    char dtime[MAX_STR_LEN];
    char hsname[MAX_STR_LEN];
    struct tm *tm;
    time_t clock;
    time (&clock);
    tm = localtime (&clock);
    strftime (dtime,MAX_STR_LEN,"%c",tm);
    //gethostname (hsname,MAX_STR_LEN); //linux
	//GetComputerName(hsname,MAX_STR_LEN);
	strcpy(hsname, "PC WINDOWS USER");

    FILE *log = fopen("run.log", "a");
    fprintf (log, "time=%s host=%s\n", dtime,hsname);
    for (int i=0; i<Nopt; i++) fprintf (log, " %s", option[i]);
    fprintf (log, "\n");
    //fprintf (log, " %s\n", arg);
    fclose(log);
}

/* routine for printing a double array "a" into stream "out" with format "fmt" */

void println(FILE *out, const char* fmt, double a[3]) {
       for (int i=0; i<3; i++) fprintf(out, fmt, a[i]);
       fprintf(out, "\n");
}

/* routine for printing a double array "a" into "stdout" with format "fmt" */

void println(const char* fmt, double a[3]) {
       for (int i=0; i<3; i++) printf(fmt, a[i]);
       printf("\n");
}

/* routine for printing a m*n double matrix "a" into stream "out" with format "fmt" */

void print_matrix(FILE *out, int m, int n, dcmplx* a, int lda ) {
        for( int i = 0; i < m; i++ ) {
            for( int j = 0; j < n; j++ )
                fprintf( out, "%16.12lf %16.12lf", real(a[i+j*lda]), imag(a[i+j*lda]) );
            fprintf( out, "\n" );
        }
}

/* routine for printing the Q vector and 
   nTHz phonon frequencies contained in "f" into stream "out" with format "fmt" */

void print_line(FILE *out, double Q[3], char *fmt, double *f, int nTHz ) {
        for( int i = 0; i < 3; i++ )
            fprintf( out, fmt, Q[i] );

        gsl_vector *eval = gsl_vector_alloc (nTHz);
	for (int i=0; i<nTHz; i++) 
	    gsl_vector_set (eval, i, f[i]);
	gsl_sort_vector(eval);

        for( int i = 0; i < nTHz; i++ )
            fprintf( out, " %8.4lf", gsl_vector_get(eval,i) );
        fprintf( out, "\n" );
}

/* rountine to reduce two integers */
char *irri(int x, int xx)
{
int i, ii, jj;
char *tmp = new char[17];
	if (x==0) {
		sprintf(tmp, "%d", x);
		return tmp;
	}
	for (i=abs(x); i>=2; i--) if ((x/i)*i==x && (xx/i)*i==xx) break;
	ii = x/i;
	jj = xx/i;
	sprintf(tmp, "%d/%d", ii, jj);
	return tmp;
}

/* routines for misc matrix operations */

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

void mproduct(double a[3][3], double *bb, double c[3][3])
{
double b[3][3];
int i, k; 
	for (i=0; i <3; i++) 
	    for (k=0; k<3; k++) b[i][k] = bb[i*3+k];
	mproduct (a, b, c);
}


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

void mproductL(double a[3], double *bb, double c[3])
{
double b[3][3];
int i, k;
        for (i=0; i <3; i++)
            for (k=0; k<3; k++) b[i][k] = bb[i*3+k];
        mproduct (a, b, c);
}

void mproductR(double *bb, double a[3], double c[3])
{
double b[3][3];
int i, k;
        for (i=0; i <3; i++)
            for (k=0; k<3; k++) b[i][k] = bb[i*3+k];
        mproduct (b, a, c);
}

/* normalize a vector "a" */

void Anormal(double *a)
{
	double N = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	if (N <= THR) N = 1.e0;
        for (int i=0; i <3; i++) a[i] /= N;
}

/* mod a vector "a" by a mode "1.0" */
double normal(double a)
{
	a = fmod(a, 1.e0);
	if (fabs(a) < THR) a = 0.e0;
	else if (a < 0.e0) a += 1.e0;
	if (fabs(a-1.e0) < THR) a = 0.e0;
	return a;
}

/* return the scalar length of a vector */
double normal(double *a)
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

/* return the sqare of the scalar length of a vector */
double normal2(double *a)
{
	return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
}

/* integerizing the double vector "a" */
int *Inormal(double a[3])
{
	int *Ix = new int[3];
	for (int i=0; i<3; i++) {
	    Ix[i] = (int) a[i];
	    if (a[i] < -THR2) Ix[i] -= 1;
	    for (int j=-2; j<=2; j++) 
		if (fabs(a[i] - (Ix[i] + j))<THR2) Ix[i] += j;
	}
        return Ix;
}

/* mod "a" by PI2 */
double ModtoPI(double a) 
{
	a=fmod(a, PI2); 
	if (a<-PI-THR) a += PI2;
	if (a>PI+THR) a -= PI2;
	return a;
}

/* return the cell volume of a 3 bv 3 matrix */

double volume(double a[3][3])
{
        return a[0][0]*a[1][1]*a[2][2]
            + a[0][1]*a[1][2]*a[2][0]
            + a[0][2]*a[1][0]*a[2][1]
            - a[0][2]*a[1][1]*a[2][0]
            - a[0][0]*a[1][2]*a[2][1]
            - a[0][1]*a[1][0]*a[2][2];
}

/* return the cell volume by the three lattice vector a, b, c */

double volume(double *a, double *b, double *c)
{
	return a[0]*b[1]*c[2]
            + a[1]*b[2]*c[0]
            + a[2]*b[0]*c[1]
            - a[2]*b[1]*c[0]
            - a[0]*b[2]*c[1]
            - a[1]*b[0]*c[2];
}

/* routines for misc vector operations */

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

/* inverse "mout" of a matrix 3 by 3 "m" */

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

/* inverse matrix "mmout" of a matrix "mm" */

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

/* check if "a" is an exact Q points by the dot product with supercell "b" */

int exactK(double a[3], double b[3][3])
{
	for (int i=0; i<3; i++) {
	    double tmp = dotproduct(a, b[i]);
		double imp = floor(tmp+0.5);
	    if (fabs(imp - tmp) > THR*THR) return 0;
	}
	return 1;
}

/* class for Linux CPU time */

class CPUTIM {

    private:
	clock_t start, middle;

    public:
	CPUTIM() {
	    start = clock();
	    middle = start;
	}

	/* out - out stream to print to */
	void elptime(FILE *out) {
	    clock_t end = clock();
	    double total = ((double) (end - start)) / CLOCKS_PER_SEC;
	    double elapsed = ((double) (end - middle)) / CLOCKS_PER_SEC;
	    middle = end;
	    fprintf (out, " Section time %.3f Sec., Total time = %.3f Sec.\n",
		elapsed, total);
		fflush(out);
	}
};


/* make the accoustic sum rule using a self-consistant loop method for Yphon option
  "-tranI 2" 
  Fij - force constant matrix (3*natomS by 3*natomS)
  natomS - number of atoms in the supercell
  wFij - weight matrix to Fij
*/

double makeTI(double *Fij, int natomS, double *wFij) {
    int dimN = natomS*3;
    int nR = 3;
    int M = dimN*nR;
    double *R = new double[M];
    double *A = new double[M];
    for (int i=0; i<dimN; i++) {
	for (int s=0; s<nR; s++) {
	    if (i%nR==s) R[i+s*dimN] = 1.e0;
	    else R[i+s*dimN] = 0.e0;
	}
    }

    for (int i=0; i<dimN; i++) {
	for (int r=0; r<nR; r++) {
	    A[i+r*dimN] = 0.e0;
	    for (int j=0; j<dimN; j++)
	    	A[i+r*dimN] -= Fij[i*dimN+j]*R[j+r*dimN];
	}
    }

    double *B = new double[M*M];
    for (int i=0; i<dimN; i++) {
	for (int r=0; r<nR; r++) {
	    for (int j=0; j<dimN; j++) {
		for (int t=0; t<nR; t++) {
		    int ii=i+r*dimN;
		    int jj=j+t*dimN;
		    double tmp = wFij[i*dimN+j]*wFij[i*dimN+j]*R[ii]*R[jj];
		    if (i==j) {
		        for (int k=0; k<dimN; k++) {
			    tmp += wFij[i*dimN+k]*wFij[i*dimN+k]*R[k+r*dimN]*R[k+t*dimN];
		        }   
		    }
		    B[ii*M+jj] = .25e0*tmp;
		}
	    }
	}
    }

       gsl_matrix_view m 
         = gsl_matrix_view_array (B, M, M);
     
       gsl_vector_view b
         = gsl_vector_view_array (A, M);
     
       gsl_vector *x = gsl_vector_alloc (M);
       
       int s;
     
       gsl_permutation * p = gsl_permutation_alloc (M);
     
       gsl_linalg_LU_decomp (&m.matrix, p, &s);

       gsl_linalg_LU_det(&m.matrix, s);
       //double det = gsl_linalg_LU_det(&m.matrix, s);
       //printf ("det = %g, s=%d\n",det, s);
     
       gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
     
    for (int i=0; i<dimN; i++) {
        for (int j=0; j<dimN; j++) {
            double tmp = 0.e0;
            for (int s=0; s<nR; s++) {
                double xi = gsl_vector_get (x, i+s*dimN);
                double xj = gsl_vector_get (x, j+s*dimN);
                tmp += xi*R[j+s*dimN]+xj*R[i+s*dimN];
            }
	    Fij[i*dimN+j] += .25e0*wFij[i*dimN+j]*wFij[i*dimN+j]*tmp;
        }
    }

       gsl_permutation_free (p);
       gsl_vector_free (x);

    double dout=0.e0;
    for (int i=0; i<dimN; i++) {
	for (int r=0; r<nR; r++) {
	    double tmp= 0.e0;
	    for (int j=0; j<dimN; j++)
	    	tmp -= Fij[i*dimN+j]*R[j+r*dimN];
	    dout += pow(tmp,2);
	}
    }
    dout = sqrt(dout/M);
    return dout;
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
    char *str(int i) {return ListString[i];} //return a field
    int operator[] (int i) {return atoi(ListString[i]);} //convert a field into integer
};

/* class defining the atomic posttions and symbols */

class ATM {

    private:

    public:
        double x, y, z;
        char *sym;

	/* set up by the explicit atomic posttions and symbols */
        ATM(double _x, double _y, double _z, const char *_sym) {
            x = _x;
            y = _y;
            z = _z;
            int s = strlen(_sym);
            sym = new char[s+1];
            //strcpy(sym, _sym);
	    int i;
	    for (i=0;i<s && isLetter(_sym[i]);i++) sym[i] = _sym[i]; 
	    sym[i] = '\000';
        }

	/* set up by a text line */
        ATM(char *line) {
	    SplitItem *pos = new SplitItem(line);
            x = atof(pos->str(0));
            y = atof(pos->str(1));
            z = atof(pos->str(2));
            if (pos->GetN()>3) {
		char *_sym = pos->str(3);
                int s = strlen(_sym);
	        int i;
	        for (i=0;i<s && isLetter(_sym[i]);i++) sym[i] = _sym[i]; 
	        sym[i] = '\000';
	    } else {
		sym = new char[2];
		strcpy(sym, "X");
	    }
        }

	/* return atomic positions */
	double *getxyz() {
	    double *xyz = new double[3];
	    xyz[0] = x;
	    xyz[1] = y;
	    xyz[2] = z;
	    return xyz;
	}

	/* set up atomic positions */
	void setxyz(double *xyz) {
	    x = xyz[0];
	    y = xyz[1];
	    z = xyz[2];
	}

	/* check if var is very small */
        int isq(double var) {
            if (fabs(var)<THR2) return 1;
            else if (fabs(fabs(var)-1.e0)<THR2) return 1;
            else return 0;
        }

	/* check if "xx,yy,zz" and "x,y,z" are the same position */
        int operator==(ATM s) {
            double xx = x - s.x;
            double yy = y - s.y;
            double zz = z - s.z;
	    if (strcmp(s.sym, sym)) return 0;
            if (isq(xx)+isq(yy)+isq(zz) == 3 ) return 1;
            else return 0;
        }

	// check if f and x are the same position
        int equalRot(ATM s) {
	    if (strcmp(s.sym, sym)) return 0;
            double xx = x - s.x; xx = fmod(xx,1.0);
	    if (xx <= -0.5) xx += 1.0; if (xx >= 0.5) xx -= 1.0;
            double yy = y - s.y; yy = fmod(yy,1.0);
	    if (yy <= -0.5) yy += 1.0; if (yy >= 0.5) yy -= 1.0;
            double zz = z - s.z; zz = fmod(zz,1.0);
	    if (zz <= -0.5) zz += 1.0; if (zz >= 0.5) zz -= 1.0;
            if (sqrt(xx*xx + yy*yy + zz*zz) < THR3 ) return 1;
            else return 0;
	}

	/* check if "xx,yy,zz" and "x,y,z" are the same position without compare atomic symbol */
        int equalSQS(ATM s) {
            double xx = x - s.x;
            double yy = y - s.y;
            double zz = z - s.z;
            if (isq(xx)+isq(yy)+isq(zz) == 3 ) return 1;
            else return 0;
        }
};

/* class to save the index of the primitive unit cell in the supercell */

class CINDEX {

    private:

    public:
	/* index primitive unit cell by the integer x, y, z and double X, Y, Z */
        int x, y, z;
	double X, Y, Z;

        CINDEX(int _x, int _y, int _z) {
            x = _x;
            y = _y;
            z = _z;
	    X = (double)x;
	    Y = (double)y;
	    Z = (double)z;
        }

	/* check if "xx,yy,zz" and "x,y,z" are the same unit cell */
        int operator==(CINDEX s) {
            int xx = x - s.x;
            int yy = y - s.y;
            int zz = z - s.z;
            if ( xx==0 && yy==0 && zz==0 ) return 1;
            else return 0;
        }

	int get_cx() {return x;} //x index of a PUC in supercel
	int get_cy() {return y;} //y index of a PUC in supercel
	int get_cz() {return z;} //z index of a PUC in supercel
	int get_cell() {return x;} //sequential index of a PUC in supercel
	int get_atomP() {return y;} //atom index in PUC
	int get_atomS() {return z;} //atom index in supercel
};

/* class of calculating and saving the atomic distance infor for speeding up calculation 
   This is an essential part making Yphon efficient and fast,
   by means of pre-calculating and saving the atomis distance infor. 
*/

class RDIS {

    private:

    public:
        double *Rx, *RR; //postion difference and distance to a "central PUC"
        double *xP; //internal postion difference within PUC
        double **IxS; // Relative distance in the unit of PUC lat. vec.
        int *NxS; //number of atoms that are in the same distance. Think abour super-supercell

	/* atomic distance for fast calculation */
	/* input
	   atomP - atomic infor at the primitive unit cell (PUC)
	   natomP - number of atoms in the PUC
	   cellA - pointer to index of PUC in the supercell, index of atom in the PUC, index of atom in the supercell
	   cellI - pointer to index containing three integers index PUC positions in the supercell
	   M - index of the center/reference atom in the supercell
	   natomS - number atoms in the supercell
	   scell_inv_avec - reciprocal lat. vec. of the supercell
	   scell - lat. vec. of supercell
	   avec - lat. vec. of PUC
	   mode - =1 stable working mode
		  =2 under developing 
	   
	   output
	   Rx = relative positions of cell M to other atoms in the supercell
	   RR = atomic distance of cell M to other atoms in the supercell
	   xP = relative positions of atoms in the PUC
	   NxS - point to number of atomic positions having the same minimum distance
	   IxS - point to the relative distance vec
	*/
        RDIS (ATM **atomP, int natomP, 
	    ATM **atomS, CINDEX **cellI, CINDEX **cellA, int M, int natomS, 
	    double scell_inv_avec[3][3], double scell[3][3], 
	    double avec[3][3], int mode) {
		double Ix[3];
		IxS = (double **) malloc ( (size_t) (natomS*(sizeof(double *))) );
		NxS = new int[natomS];
		for (int i=0; i<natomS; i++) NxS[i] = 0;

		Rx = new double[natomS*3];
		RR = new double[natomS];
		xP = new double[natomP*3];

		int i = cellA[M]->get_atomP();
		for (int N=0; N<natomP; N++) {
		    xP[N*3] = atomP[N]->x - atomP[i]->x;
		    xP[N*3+1] = atomP[N]->y - atomP[i]->y;
		    xP[N*3+2] = atomP[N]->z - atomP[i]->z;
		}

		double xA = atomS[M]->x;
		double yA = atomS[M]->y;
		double zA = atomS[M]->z;
		int iM = cellA[M]->get_cell();
		int xM = cellI[iM]->get_cx();
		int yM = cellI[iM]->get_cy();
		int zM = cellI[iM]->get_cz();
                for (int N=0; N<natomS; N++) {
		    double a[3], b[3];
		    if (mode==2) getReduceDIS(Ix, a,xA,yA,zA,N, scell, atomS);
		    else getReduceDIS(Ix, a,xA,yA,zA,N, scell, atomS, NRDIS);
		    mproduct(a, scell_inv_avec, b);
		    mproduct(b, avec, a);
		    RR[N] = normal(a);

		    int jN = cellA[N]->get_cell();
		    mproduct(Ix, scell_inv_avec, b);
	 	    Rx[N*3] = b[0] + (double)(cellI[jN]->get_cx() - xM);
	 	    Rx[N*3+1] = b[1] + (double)(cellI[jN]->get_cy() - yM);
	 	    Rx[N*3+2] = b[2] + (double)(cellI[jN]->get_cz() - zM);

		    if (mode==2) getReduceDIS(RR[N],N,Ix, xA,yA,zA, scell, 
			scell_inv_avec, atomS, natomS, cellA, 
			xM, yM, zM, cellI);
		    else getReduceDIS(RR[N],N,Ix, xA,yA,zA, scell, 
			scell_inv_avec, atomS, natomS, cellA, 
			xM, yM, zM, cellI, NRDIS);
		}
        }

	double *get_Rx(int i) {return Rx+3*i;} //return pointer of relative atomic positions in supercell
	double get_RR(int i) {return RR[i];} //return relative distance in supercell
	double *get_xP(int i) {return xP+3*i;} //return pointer of relative atomic positions in PUC

	/* calculate the Reduced distance dd and index Ix */
	/* x,y,z - given position
	   Ix - relative position index of atom N to x, y, z for minimum distance
	   N - index of relative atom in the supercell
	   scell - lat. vec. of supercell
	   atomS - atomic infor in the supercell
	   dd - minimum relative position vector
	   nrd - step range of super/mirror of supercell for looking for mininum distance
		(think about side, face, corner etc)
	*/
        static void getReduceDIS(double Ix[3], double dd[3],
            double x, double y, double z, int N, 
	    double scell[3][3], ATM **atomS, int nrd=0) {
            int i, j, k;
            double R;
            double a[3], b[3];
            R = 1.e72;
            for (i=-nrd; i<=nrd; i++) { // nrd - the boundary to search for mininum distance
                a[0] = (double)i + atomS[N]->x - x;
                for (j=-nrd; j<=nrd; j++) {
                    a[1] = (double)j + atomS[N]->y - y;
                    for (k=-nrd; k<=nrd; k++) {
                        a[2] = (double)k + atomS[N]->z - z;
			mproduct (a, scell, b);
                        double dis = normal2(b);
                        if (dis < R) {
                            R = dis;
                            dd[0] = a[0];
                            dd[1] = a[1];
                            dd[2] = a[2];
                            Ix[0] = (double)i;
                            Ix[1] = (double)j;
                            Ix[2] = (double)k;
                        }
                    }
                }
            }
        }

	/* set up the relative position vector IxS for the minimum distances
		(might more than one, esspecially important for interpolation 
		in the vincinity  of Gamma point phonon DOS
		think about side, face, corner etc)
	   R - minimum distance
	   N - index of the relative atom in the supercell
	   xI - relative index of atom N to x, y, z for minimum distance
	   x,y,z - given position
	   scell - supercel lat. Vec.
	   scell_inv_avec - supercell lat. vec. vs PUC lat. vec.
	   xM, yM, zM - relative positions
	   nrd - step range of super/mirror of supercell for looking for mininum distance
	*/
        void getReduceDIS(double R, int N, double xI[3], 
	    double x, double y, double z,
            double scell[3][3], double scell_inv_avec[3][3],
	    ATM **atomS, int natomS, CINDEX **cellA,
	    int xM, int yM, int zM, CINDEX **cellI, int nrd=0) {
            int i, j, k;
            double a[3], b[3], Ix[3];
            for (i=-nrd; i<=nrd; i++) {
                for (j=-nrd; j<=nrd; j++) {
                    for (k=-nrd; k<=nrd; k++) {
                	a[0] = (double)i + atomS[N]->x - x;
                    	a[1] = (double)j + atomS[N]->y - y;
                        a[2] = (double)k + atomS[N]->z - z;
                        mproduct (a, scell, b);
                        double dis = normal(b);
                        if (fabs(dis-R) < THR3) {
//printf ("Rdis called with R = %lf, N= %d, xyz = %lf %lf %lf\n", R, N, x, y, z);
                            Ix[0] = (double)i - xI[0];
                            Ix[1] = (double)j - xI[1];
                            Ix[2] = (double)k - xI[2];
                            mproduct (Ix, scell_inv_avec, b);
			    if (NxS[N]==0) IxS[N] = (double *) malloc ( (size_t) 
				(3*(sizeof(double))) );
			    else IxS[N] = (double *) realloc (IxS[N], (size_t) 
				(3*(NxS[N]+1)*(sizeof(double))) );

			    double *p = IxS[N]+NxS[N]*3;
			    p[0] = b[0]; 
			    p[1] = b[1]; 
			    p[2] = b[2]; 
			    NxS[N]++;
                        }
                    }
                }
            }
        }

	int get_NxS(int n) {return NxS[n];} //return number of atomic positions having the same minimum distance
	double *get_IxS(int N, int n) {return IxS[N]+n*3;} //return the pointer to the relative distance vec 
};

/* base class for CijRDIS */

class SymRDIS {

    private:

    public:
        double **IxS, *RR;
        int *NxS;

        SymRDIS () {}

	int get_NxS(int n) {return NxS[n];}

	double *get_IxS(int N, int n) {return IxS[N]+n*3;}
	double get_RR(int N) {return RR[N];}
};

/* class symmetrizing and saving the atomic distance infor for speeding up calculation 
	see class RDIS for the meanings of all variables
	this class will be deleted due to low accuracy in calculating elastic constants from
	force constants
*/

class CijRDIS : public SymRDIS {

    private:

    public:
	/* atomic distance for fast calculation for elastic constant calculation */
        /* input
           atomS - atomic infor within the supercell
           natomS - number atoms in the supercell
           M - index of the center/reference atom in the supercell
           scell - lat. vec. of supercell

           output
	   Rx = relative positions of cell M to other atoms in the supercell
	   RR = atomic distance of cell M to other atoms in the supercell
	   xP = relative positions of atoms in the PUC
	   NxS - point to number of atomic positions having the same minimum distance
	   IxS - point to the relative distance vec
        */

        CijRDIS (ATM **atomS, int natomS, double scell[3][3], int M) {
		IxS = (double **) malloc ( (size_t) (natomS*(sizeof(double *))) );
		RR = new double[natomS];
		NxS = new int[natomS];
		for (int i=0; i<natomS; i++) NxS[i] = 0;

                for (int N=0; N<natomS; N++) {
		    RR[N] = getReduceDIS(N, M, scell, atomS, NRDIS/2);
		    getReduceDIS(RR[N], N, M, scell, atomS, natomS, NRDIS/2);
		}
        }

        /* atomic distance for fast calculation for DM calculation */
        /* input
           atomS - atomic infor within the supercell
           natomS - number atoms in the supercell
           M - index of the center/reference atom in the supercell
           scell - lat. vec. of supercell
           inv_avec - reciprocal lat. vec. of the PUC
        */

        CijRDIS (ATM **atomS, int natomS, double scell[3][3], double inv_avec[3][3], int M) {
		IxS = (double **) malloc ( (size_t) (natomS*(sizeof(double *))) );
		RR = new double[natomS];
		NxS = new int[natomS];
		for (int i=0; i<natomS; i++) NxS[i] = 0;

                for (int N=0; N<natomS; N++) {
		    RR[N] = getReduceDIS(N, M, scell, atomS, NRDIS);
		    getReduceDIS(RR[N], N, M, scell, atomS, natomS, NRDIS, inv_avec);
		}
        }

        /* find the minimum distance between M and N
           M - index of the center/reference atom in the supercell
           N - index of the relative atom in the supercell
           scell - lat. vec. of supercell
           atomS - atomic infor within the supercell
	   nrd - step range of super/mirror of supercell for looking for mininum distance
        */
        double getReduceDIS(int N, int M, double scell[3][3], 
	    ATM **atomS, int nrd=0) {
            int i, j, k;
            double R;
            double a[3], b[3];
            R = 1.e72;
            for (i=-nrd; i<=nrd; i++) {
                a[0] = (double)i + atomS[N]->x -  atomS[M]->x;
                for (j=-nrd; j<=nrd; j++) {
                    a[1] = (double)j + atomS[N]->y - atomS[M]->y;
                    for (k=-nrd; k<=nrd; k++) {
                        a[2] = (double)k + atomS[N]->z - atomS[M]->z;
			mproduct (a, scell, b);
                        double dis = normal(b);
                        if (dis < R) R = dis;
                    }
                }
            }
	    return R;
        }

        /* calculate the relative distance vec between M and N
	   R - the minimum distance between M and N
           M - index of the center/reference atom in the supercell
           N - index of the relative atom in the supercell
           scell - lat. vec. of supercell
           atomS - atomic infor within the supercell
	   nrd - step range of super/mirror of supercell for looking for mininum distance
	   inv_avec - reciprocal lat. vec. of the PUC
	   IxS - containing the relative coordinates
        */
        void getReduceDIS(double R, int N, int M,
            double scell[3][3], ATM **atomS, int natomS, int nrd, double inv_avec[3][3]) {
            int i, j, k;
            double a[3], b[3];
            for (i=-nrd; i<=nrd; i++) {
                for (j=-nrd; j<=nrd; j++) {
                    for (k=-nrd; k<=nrd; k++) {
                	a[0] = (double)i + atomS[N]->x - atomS[M]->x;
                    	a[1] = (double)j + atomS[N]->y - atomS[M]->y;
                        a[2] = (double)k + atomS[N]->z - atomS[M]->z;
                        mproduct (a, scell, b);
                        double dis = normal(b);
                        //if (fabs(dis-R) < THR+THR3*R) {
                        if (fabs(dis-R) < THR3) {
//printf ("CijRdis called with R = %lf,  THR+THR3*R= %lf\n", R, THR+THR3*R);
			    if (NxS[N]==0) IxS[N] = (double *) malloc ( (size_t) 
				(3*(sizeof(double))) );
			    else IxS[N] = (double *) realloc (IxS[N], (size_t) 
				(3*(NxS[N]+1)*(sizeof(double))) );

			    double *p = IxS[N]+NxS[N]*3;
                            mproduct (b, inv_avec, a);
			    p[0] = a[0]; 
			    p[1] = a[1]; 
			    p[2] = a[2]; 
			    NxS[N]++;
                        }
                    }
                }
            }
        }

        /* calculate the relative distance vec between M and N
	   R - the minimum distance between M and N
           M - index of the center/reference atom in the supercell
           N - index of the relative atom in the supercell
           scell - lat. vec. of supercell
           atomS - atomic infor within the supercell
	   nrd - step range of super/mirror of supercell for looking for mininum distance
	   IxS containing the Cartesian coordinates
        */
        void getReduceDIS(double R, int N, int M,
            double scell[3][3], ATM **atomS, int natomS, int nrd=0) {
            int i, j, k;
            double a[3], b[3];
            for (i=-nrd; i<=nrd; i++) {
                for (j=-nrd; j<=nrd; j++) {
                    for (k=-nrd; k<=nrd; k++) {
                        a[0] = (double)i + atomS[N]->x - atomS[M]->x;
                        a[1] = (double)j + atomS[N]->y - atomS[M]->y;
                        a[2] = (double)k + atomS[N]->z - atomS[M]->z;
                        mproduct (a, scell, b);
                        double dis = normal(b);
                        if (fabs(dis-R) < THR) {
                            if (NxS[N]==0) IxS[N] = (double *) malloc ( (size_t)
                                (3*(sizeof(double))) );
                            else IxS[N] = (double *) realloc (IxS[N], (size_t)
                                (3*(NxS[N]+1)*(sizeof(double))) );

                            double *p = IxS[N]+NxS[N]*3;
                            p[0] = b[0];
                            p[1] = b[1];
                            p[2] = b[2];
                            NxS[N]++;
                        }
                    }
                }
            }
        }
};

/* class for vector operations */

class VECTOR {

    private:
	double v[3];

    public:
	VECTOR () {}

	/* define the vector using atomic positions "pos" defined in ATM class */
	VECTOR (ATM *pos) {
	    v[0] = pos->x;
            v[1] = pos->y;
            v[2] = pos->z;
	}

	/* define the vector by explicitly providing atomic positions _v as vector */
	VECTOR (double _v[3]) {
	    for (int i=0; i<3; i++) v[i] = _v[i];
	}

	/* define the vector by explicitly providing atomic positions */
	VECTOR (double _x, double _y, double _z) {
	    v[0] = _x;
            v[1] = _y;
            v[2] = _z;
	}

	double operator [] (int i) {return v[i];}

	/* normalize the vector with floor and ceil */
	void Ceil(double floor, double ceil) {
	    for (int i=0; i<3; i++) 
		while (v[i] >= ceil-THR) v[i] -= 1.e0;

	    for (int i=0; i<3; i++) 
		while (v[i] <= floor-THR) v[i] += 1.e0;
	}

	void operator = (VECTOR _v) { //copy a vector
	    for (int i=0; i<3; i++) v[i] = _v.v[i];
	}  

	VECTOR operator + (VECTOR _v) { //vector plus
	    VECTOR *tmp = new VECTOR();
	    for (int i=0; i<3; i++) tmp->v[i] = v[i] + _v.v[i];
	    return *tmp;
	}

	VECTOR operator - (VECTOR _v) { //vector minus
	    VECTOR *tmp = new VECTOR();
	    for (int i=0; i<3; i++) tmp->v[i] = v[i] - _v.v[i];
	    return *tmp;
	}

	VECTOR operator += (VECTOR _v) { //adding a vector
	    VECTOR *tmp = new VECTOR();
	    for (int i=0; i<3; i++) tmp->v[i] = v[i] + _v.v[i];
	    return *tmp;
	}

	VECTOR operator -= (VECTOR _v) { //minusing a vector
	    VECTOR *tmp = new VECTOR();
	    for (int i=0; i<3; i++) tmp->v[i] = v[i] - _v.v[i];
	    return *tmp;
	}

	operator double * (void) { return v; };
};


/* will be deleted after the pos2s analysis is developped */

typedef struct {
    double MoM[10];
} RAMAN;

/* will be deleted after the pos2s analysis is developped */

class GammaPhonon {

    private:
    ATM **atomP;
    int natomP;
    double *mass;
    double *Rv;
    int nTHz;
    double *pos;

    public:
    GammaPhonon(ATM **_atomP, int _natomP, double *_Rv, double avec[3][3]) {
	atomP = _atomP;
	natomP = _natomP;
	nTHz = 3*natomP;
	Rv = _Rv;
	pos = new double[nTHz];
	mass = new double[natomP];
        
	double tmass = 0.e0;
	for (int j=0; j<natomP; j++) {
		mproduct(atomP[j]->getxyz(), avec, pos+3*j);
		mass[j] = getmass(atomP[j]->sym);
		tmass += mass[j];
	}

        double xc[3];
	for (int j=0; j<3; j++) xc[j] = 0.e0;
	double sum = 0.e0;
	for (int j=0; j<nTHz; j++) {
		xc[j%3] += mass[j/3]*pos[j];
		Rv[j] /= sqrt(mass[j/3]);
		sum += Rv[j]*Rv[j];
	}
	for (int j=0; j<3; j++) xc[j] /= tmass;

	for (int j=0; j<nTHz; j++) {
		Rv[j] /= sqrt(sum);
		pos[j] -= xc[j%3];
	}
//printf(" XC = %8.4lf %8.4lf %8.4lf\n", xc[0], xc[1], xc[2]); 
    }

    double ByDipole(double *BCharg, double qvec[3][3]) {
	double RX[3], tmp[3];
	for (int j=0; j<3; j++) RX[j] = 0.e0;
	for (int N=0; N<natomP; N++) {
		mproduct(Rv+N*3, BCharg+N*9, tmp);
		for (int j=0; j<3; j++) RX[j] += tmp[j];
	}
	mproduct(dirE, qvec, tmp);
	Anormal(tmp);
	return (dotproduct(RX, tmp));
    }

    double BySymm() {
        double X=0.e0, Y=0.e0, Z=0.e0;
        for (int N=0; N<natomP; N++) {
                X += Rv[N*3];
                Y += Rv[N*3+1];
                Z += Rv[N*3+2];
        }
        return sqrt(X*X +Y*Y+Z*Z);
    }

    void ByRaman(RAMAN *raman) {
        raman->MoM[0] = GSum(0);
        raman->MoM[1] = GSum(1);
        raman->MoM[2] = GSum(2);
        raman->MoM[3] = GSum(3);
        raman->MoM[4] = GSum(4);
        raman->MoM[5] = GSum(5);
        raman->MoM[6] = GSum(6);
        raman->MoM[7] = GSum(7);
        raman->MoM[8] = GSum(8);
        raman->MoM[9] = GSum(9);
    }

    double GSum(int key) {
        double X=0.e0, Y=0.e0, Z= 0.e0;
	double a[3];
        for (int N=0; N<natomP; N++) {
	    cross_product(Rv+N*3, pos+N*3, a);
	    double R;
	    if (key==1) {
		R = pos[N*3];
	    } else if (key==2) {
		R = pos[N*3+1];
	    } else if (key==3) {
		R = pos[N*3+2];
	    } else if (key==4) {
		R = pos[N*3]*pos[N*3] - pos[N*3+1]*pos[N*3+1];
	    } else if (key==5) {
		R = 2.e0* pos[N*3+2]*pos[N*3+2] - pos[N*3]*pos[N*3] - pos[N*3+1]*pos[N*3+1];
	    } else {
		R =1.e0;
	    }
            X += R*a[0];
            Y += R*a[1];
            Z += R*a[2];
        }
        return sqrt(X*X+Y*Y+Z*Z);
    }
};

/* Major class of the mixed space approach */

class LatDYN {

    private:
	/* gsl's are the gsl running memory for GSL */
        static gsl_vector *eval;
        static gsl_matrix_complex *evec;
        static gsl_matrix_complex *m;
	static gsl_eigen_hermv_workspace *wv;
	static gsl_eigen_herm_workspace *w;
        static gsl_matrix *evecR;
        static gsl_matrix *mR;
	//static gsl_eigen_symmv_workspace *wR;
	static gsl_eigen_symm_workspace *wR;

        static double *runa, *runb, *runc;

    public:
	static double avec[3][3], qvec[3][3], scell[3][3]; //lattice vectors of the PUC and supercell
	static double inv_avec[3][3], scell_inv_avec[3][3]; //inverse matrix of PUC and supercell
	static int nTHz; //number of phonon branches which is three time of the number of atoms in the PUC
	static double *frequencies;
	static int Nfreq;//total number of frequencies
	static int Debug;//debug control
	double vecq[3], vecQ[3]; //Q vector
	double *eigenvalue; //eigen frequencies
	static dcmplx *dmat; //total dynamical matrix
	static dcmplx *dmatNAG; //dynamical matrix from the dipole-dipole interaction
	static double *dmatG; //Gamma point dynamical matrix
	static RDIS **disR; //speeding up memory
	static SymRDIS **disNR; //speeding up memory
	static CINDEX **cellA; //primitive unit cell index
	static double *Fij; //real space force constants BY VASP.5
	static double BornCC; //Born effective charge tensor

	static int pvdosN; //partial phonon DOS calculation
	static int *xeN;
        static double *xe;
	static double *pvdos;

	/* defining Q point vector by two forms
	  x,y,z - Q point
	*/
        LatDYN (double x, double y, double z) {
	    vecq[0] = x*PI2;
	    vecq[1] = y*PI2;
	    vecq[2] = z*PI2;
            vecQ[0] = x;
            vecQ[1] = y;
            vecQ[2] = z;
        }

	/* defining Q point vector by product with a matrix lat
	  x,y,z - Q point
	  lat - transformation matrix
	 */
        LatDYN (double x, double y, double z, double *lat) {
	    vecQ[0] = x;
	    vecQ[1] = y;
	    vecQ[2] = z;
	    mproduct (vecQ, lat, vecq);
	    for (int i=0; i<3; i++) vecq[i] *= PI2;
        }

	/* prepare the workspace lattice dynamics calculation
	  avec - lat. vec. of PUC
	  scell - lat. vec. of supercell
	  atomP - atimic infor (positions and atomic symbols etc) within PUC
	  atomS - atimic infor (positions and atomic symbols etc) within supercell
	  natomS - number of atoms in the supercell
	  natomP - number of atoms in the PUC
	  cellA - pointer to index of PUC in the supercell, index of atom in the PUC, index of atom in the supercell
	  cellI - pointer to index containing three integers index PUC positions in the supercell
	  Fij - force constant matrix
	  MM - array containing indexes of primitive unit cells to be evaluated as the 0th PUC (reference PUC)
	  NX - number of primitive unit cells (PUC) to be averaged on
	  mode - calculation mode (do not set it), 
	  emin, emax - define the frequency range for partial DOS calculation
	*/
	static void SetWorkSpace(double _avec[3][3], double _scell[3][3],
	    ATM **atomS, ATM **atomP, CINDEX **cellI, CINDEX **_cellA, 
	    double *_Fij, int natomS, int natomP, int *MM, int NX, int mode,
	    double emin, double emax) {
	    cellA = _cellA;
	    Fij = _Fij;
	    for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
		    avec[i][j]=_avec[i][j];
		    scell[i][j]=_scell[i][j];
		}
	    }
	    inv_m(avec, inv_avec); //inverse of the supercell lattive vectors
	    mproduct(scell, inv_avec, scell_inv_avec);
            double fac = PI2/volume(avec);
	    /* The factor VASPtoeVA is from e**2/Bohr = Hartree 
				resulting e**2 = Hartree*Bohr */
	    double VASPtoeVA = 27.211384523e0 * 0.5291772085936e0;
	    BornCC = (fac+fac)* VASPtoeVA;

	    /* make the reciprocal PUC by vector cross product */
            cross_product(avec[1], avec[2], qvec[0]);
            cross_product(avec[2], avec[0], qvec[1]);
            cross_product(avec[0], avec[1], qvec[2]);

            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                        qvec[i][j] *= fac/PI2;
                }
            }

	    /* precalculating atomic distances for speeding up calculation */
	    disR = (RDIS **) malloc ( (size_t) (NX*sizeof(RDIS *)) );
	    if (!CTLDM) disNR = (SymRDIS **) malloc ( (size_t) (NX*sizeof(CijRDIS *)) );
	    for (int nx=0; nx<NX; nx++) {
		disR[nx] = new RDIS(atomP, natomP, atomS, cellI, cellA, 
		    MM[nx], natomS, scell_inv_avec, scell, avec, mode); //fast distance index
		if (!CTLDM) disNR[nx] = new CijRDIS(atomS, natomS, scell, inv_avec, MM[nx]);
	    }

	    nTHz = natomP*3; //this is the number of phonon branches
            eval = gsl_vector_alloc (nTHz);
            dmat = new dcmplx[nTHz*nTHz]; //dynamical matrix
            dmatNAG = new dcmplx[nTHz*nTHz]; //Cochran and Cowley matrix for the nonanalytical term
            dmatG = new double[nTHz*nTHz]; //used by high accurate G point calculation
	    runa = new double[nTHz]; //running memory
	    runb = new double[nTHz];
	    runc = new double[nTHz];
            evec = gsl_matrix_complex_alloc (nTHz, nTHz);
            m = gsl_matrix_complex_alloc (nTHz, nTHz);
	    w = gsl_eigen_herm_alloc (nTHz);
	    wv = gsl_eigen_hermv_alloc (nTHz);
            evecR = gsl_matrix_alloc (nTHz, nTHz);
            mR = gsl_matrix_alloc (nTHz, nTHz);
	    //wR = gsl_eigen_symmv_alloc (nTHz);
	    wR = gsl_eigen_symm_alloc (nTHz);

	    pvdosN = 12001;
	    xeN = new int[pvdosN];
	    xe = new double[pvdosN];
	    pvdos = new double[pvdosN*natomP];
	    double x = emin;
	    double de = (emax-emin)*1.1e0/(double)(pvdosN-1);
	    for (int i=0; i<pvdosN; i++) {
		xeN[i] = 0;
		xe[i] = x;
		x += de;
	        for (int j=0; j<natomP; j++) pvdos[i*natomP+j] = 0.e0;
	    }
	}

	/* make the dynamic memory saving the frequencies
	  nFrq - number of frequencies to be stored
	 */

	void static setNFrq(int nFrq) {
	    Nfreq = 0;
	    frequencies = new double[nFrq];
	}

	/* making the dynamical matrix (DM)
	  natomS - number of atoms in the supercell
	  MM - array containing indexes of primitive unit cells to be evaluated as the 0th PUC (reference PUC)
	  NX - number of primitive unit cells to be averaged on
	  FACTOR - natomP/NX, where natomP is number of atoms in PUC
	  RedDM - control for if using the reduced dynamical matrix
	  tranI - control for ASR
	  calcNA - control for if couting the D-D interaction
	  dmat - dynamical matrix
	*/
        void ReduceDM(int natomS, int *MM, int NX, double FACTOR, 
	    int RedDM=1, int tranI=0, int calcNA=0) {
            for (int i=0; i<nTHz*nTHz; i++) dmat[i] = 0.e0;

	    int dimN = natomS*3;
            for (int nx=0; nx<NX; nx++) {
                int M = MM[nx];
                int i = cellA[M]->get_atomP();
                for (int N=0; N<natomS; N++) {
                    int j = cellA[N]->get_atomP();
		    if (j>i) continue;
		    //if (disR[nx]->get_RR(N) > DIS)continue;
		    int Np = disR[nx]->get_NxS(N);

		    double bQ = dotproduct(disR[nx]->get_Rx(N),vecq);
		    if (RedDM==0) 
		        bQ += dotproduct(disR[nx]->get_xP(j),vecq);
                    dcmplx fac = polar(FACTOR,bQ);

		    if (Np > 1) {
			double *ff = disR[nx]->get_IxS(N,0);
			dcmplx f2 = 0.e0;
			for (int f = 0; f < Np; f++) {
			    f2 += polar(1.e0, dotproduct(ff,vecq));
			    ff += 3;
			}   
			fac *= f2/(double)(Np);
		    }

                    int kk = i*3*nTHz+j*3;
                    int MMNN = M*3*dimN+N*3;
                    for (int x=0; x<3; x++) {
                        for (int y=0; y<3; y++) {
                            if (calcNA ==0) dmat[kk+y] += Fij[MMNN+y]*fac;
                            else dmat[kk+y] += (dmatNAG[kk+y]+Fij[MMNN+y])*fac;
			    if (tranI) dmat[kk+y] += dmatG[kk+y]*fac;
/*
                            dmat[kk+y] += Fij[MMNN+y]*fac;
			    dmat[kk+y] += dmatG[kk+y]*fac;
*/
                        }
			kk += nTHz;
			MMNN += dimN;
                    }
                }
            }
        }

	/* Simplified non-reduced DM calculation
	  natomS - number of atoms in the supercell
	  MM - primitive unit cells to be evaluated as the 0th PUC (reference PUC)
	  NX - number of primitive unit cells (PUC) to be averaged on
	  FACTOR - natomP/NX, where natomP is number of atoms in PUC
	  tranI - control for ASR
	  calcNA - control for if couting the D-D interaction
	  dmat - dynamical matrix
	  sigma - Parlinski empirical parameter, will be deleted
	  sqs - SQS calculation (obselete)
	*/

        void NonReduceDM(int natomS, int *MM, int NX, double FACTOR, 
	    int tranI=0, int calcNA=0, double sigma=0.e0, int sqs=0) {
            for (int i=0; i<nTHz*nTHz; i++) dmat[i] = 0.e0;

	    int dimN = natomS*3;
            for (int nx=0; nx<NX; nx++) {
                int M = MM[nx];
                int i = cellA[M]->get_atomP();
                for (int N=0; N<natomS; N++) {
                    int j = cellA[N]->get_atomP();
		    if (j>i && sqs==0) continue;
		    int Np = disNR[nx]->get_NxS(N);

                    dcmplx fac = 0.e0;
		    double *ff = disNR[nx]->get_IxS(N,0);
		    double gg = 0.e0;
		    for (int np=0; np<Np; np++) {
			double bQ = dotproduct(ff,vecq);
			if (sigma != 0.e0 ) {
			    /* complete obselete (a test for Parliski empirical interpolation for LO-TO splitting) */
            		    double a[3], b[3];
			    mproduct (ff, avec, a);
			    mproduct (vecq, qvec, b);
			    double dd = normal(a) - disNR[nx]->get_RR(N);
			    double rf = exp(-sigma*sigma*dd*normal(b));
                            fac += polar(rf, bQ);
			    gg += rf;
			} else {
                            fac += polar(1.e0, bQ);
			    gg += 1.e0;
			}
			ff += 3;
		    }
                    fac *= FACTOR/gg;

                    int kk = i*3*nTHz+j*3;
                    int MMNN = M*3*dimN+N*3;
                    for (int x=0; x<3; x++) {
                        for (int y=0; y<3; y++) {
                            if (calcNA ==0) dmat[kk+y] += Fij[MMNN+y]*fac;
                            else dmat[kk+y] += (dmatNAG[kk+y]+Fij[MMNN+y])*fac;
			    if (tranI) dmat[kk+y] += dmatG[kk+y]*fac;
                        }
			kk += nTHz;
			MMNN += dimN;
                    }
                }
            }
	    if (sqs!=0) {
		/* for sqs only, under development or obselete */
            	for (int i=0; i<nTHz; i++) {
                    for (int j=i; j<nTHz; j++) {
			int ij = i*nTHz + j;
			int ji = j*nTHz + i;
			dmat[ij] = 0.5e0*(dmat[ij]+conj(dmat[ji]));
			dmat[ji] = conj(dmat[ij]);
		    }
	    	}
	    }
        }

	/* make dynamical matrix for D-D interaction 
	  atomP - atomic infor in primitive unit cell 
	  natomP - number of atoms in the primitive unit cell
	  dielec - dielectric constant tensor
          BCharg - Born effective charge tensor
	  mode - calculation mode (do not set it), 
          Kcell - number of the primitive unit cell in the supercell
	  Dir - direction toward the Gamma point
	  dmatNAG - Cochran and Cowley matrix for the nonanalytical term
	*/

	void BornGamma(ATM **atomP, int natomP, double *dielec, double *BCharg, 
	    int mode, int Kcell, double *Dir) {
	    double qq[3], qa[3], qb[3], a[3], tmpq[3];
	    for (int i=0; i<3; i++) tmpq[i]= vecq[i];

	    double tmp = normal(tmpq);
	    /* tmpq is to define the direction to the Gamma point */
            if (tmp == 0.e0) {
	        tmp = normal(Dir);
		if (tmp == 0.e0 ) for (int i=0; i<3; i++) tmpq[i]=.1e0;
                else for (int i=0; i<3; i++) tmpq[i]=.1e0*Dir[i];
	    } 

	    /* calculate the denomintor - the dielectric part */
	    double qEq, fac;
	    if (DMODE==0) {
	        mproduct(tmpq, qvec, qq);
	        mproduct(qq, dielec, a);
	        qEq = dotproduct(a,qq);
	      if (fabs(qEq) < 1.e-12) {
		printf("********WARNING, point %lf %lf %lf, tmpq= %lf %lf %lf not handled Born\n\n",
		    vecq[0], vecq[1], vecq[2], tmpq[0], tmpq[1], tmpq[2]);
		printf(" qq vector is :\n\n");
		println(" %10.6lf", qq);
		printf(" a vector is :\n\n");
		println(" %10.6lf", a);
		printf(" die tensor is :\n\n");
		for (int i=0;i<3;i++) println(" %10.6lf", dielec+i*3);
		return;
	      }
	    } else {
		/* testing, not USED */
	        mproduct(tmpq, qvec, qq);
	        mproduct(qq, inv_dielec, a);
	        qEq = dotproduct(a,qq);
	        tmp  = normal(qq);
		qEq /= pow(tmp,4);
	    }

	    for (int t=0; t<nTHz*nTHz; t++) dmatNAG[t] = 0.e0;
	    if (DMODE==0) {
	        fac = BornCC/qEq/(double)Kcell;
	    } else {
		/* testing, not USED */
	        //fac = BornCC*qEq/(double)Kcell;
	        fac = BornCC/(double)Kcell;
	    }

/*
*/
	    dcmplx *sij = new dcmplx[9];
	    if (Debug) for (int ii=0;ii<9;ii++) sij[ii] = 0.e0;

	    /* calculate the BEC part */
	    for (int i=0; i<natomP; i++) {
	    	double imass = getmass(atomP[i]->sym); //mass of ith atom in PUC
	        mproduct(qq, BCharg+i*9, qa);
	    	for (int j=0; j<natomP; j++) {
		    if (!Debug && mode!=3) if (j>i) continue;
	            double jmass = getmass(atomP[j]->sym); //mass of jth atom in PUC
	            dcmplx fac1 = fac*sqrt(1.e0/(imass*jmass));
		    if (DMODE==0) {
		      if (mode==0) {
                        double bQ = (atomP[j]->x-atomP[i]->x)*vecq[0]
                            + (atomP[j]->y-atomP[i]->y)*vecq[1]
                            + (atomP[j]->z-atomP[i]->z)*vecq[2];
                        fac1 *= polar(1.e0, -bQ);
		      }

	              mproduct(qq, BCharg+j*9, qb);
	              //mproductR(BCharg+j*9, qq, qb);
		      for (int x=0; x<3; x++) {
			int ix = i*3 + x;
                        for (int y=0; y<3; y++) {
			    int jy = j*3 + y;
                            dmatNAG[ix*nTHz+jy] = -fac1*qa[x]*qb[y];
			    if (Debug) sij[3*x+y] += 
				-fac1*qa[x]*qb[y]/sqrt(1.e0/(imass*jmass));
                        }
                      }
		    } else {
/*
	              mproduct(qq, BCharg+j*9, qb);
		      for (int x=0; x<3; x++) {
			int ix = i*3 + x;
                        for (int y=0; y<3; y++) {
			    int jy = j*3 + y;
                            dmatNAG[ix*nTHz+jy] = -fac1*qa[x]*qb[y];
                        }
		      }
*/
		      /* obselete, do not use */
		      double *Bi = BCharg+i*9;
		      double *Bj = BCharg+j*9;
		      for (int x=0; x<3; x++) {
			int ix = i*3 + x;
                        for (int y=0; y<3; y++) {
			    int jy = j*3 + y;
			    for (int xx=0; xx<3; xx++) {
			      for (int yy=0; yy<3; yy++) {
                                dmatNAG[ix*nTHz+jy] += 
				  -fac1*Bi[xx*3+x]*inv_dielec[xx*3+yy]*Bj[yy*3+y];
				  //-fac1*qa[xx]*inv_dielec[xx*3+yy]*qb[yy];
			      }  
			    }
                        }
                      }
		    }  
                }
	    }

	    if (Debug) {
                printf("\nSummation of the nonanalytical term\n\n");
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) printf("(%12.6f, %12.6f) ",
			real(sij[i*3+j]), imag(sij[i*3+j]));
		    printf("\n");
		}
	    }
	}


#ifdef GSL
	/* major eigenvalue solver using GSL
	  frequencies are aved in the array "frequencies" */
	void solveDM() {
	    double sR=0.e0, sI=0.e0;
	    for (int i=0; i<nTHz; i++) {
		for (int j=0; j<nTHz; j++) {
		    if (j>i) continue;
		    double R = real(dmat[i*nTHz+j]);
		    double I = imag(dmat[i*nTHz+j]);
		    //if (Debug) printf(" (%.4lf,%.4lf)\n", R, I);
		    gsl_complex z;
		    GSL_SET_COMPLEX (&z, R, I);
		    gsl_matrix_complex_set (m,i,j,z);
		    gsl_matrix_set (mR,i,j,R);
		    sR += fabs(R);
		    sI += fabs(I);
	        }
            }

	    if (1.e-8*sR>sI) {
	        //if (Debug) printf(" Solving REAL, real = %.8lf, imag = %.8lf\n", sR, sI);
	        //gsl_eigen_symmv (mR, eval, evecR, wR);
	        gsl_eigen_symm (mR, eval, wR);
	    } else {
	        //if (Debug) printf(" Solving IMAG, real = %.8lf, imag = %.8lf\n", sR, sI);
	        //gsl_eigen_hermv (m, eval, evec, w);
	        gsl_eigen_herm (m, eval, w);
	    }
            //gsl_eigen_hermv_free (w);
            //gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

            for (int i=0;i<nTHz;i++) {
                double evalue = gsl_vector_get (eval, i);
                if (evalue>=0.e0) evalue = -sqrt(evalue)*toTHz;
                else evalue = sqrt(-evalue)*toTHz;
		//if (vecq[0]==0.e0 && vecq[1]==0.e0 && vecq[2]==0.e0
		//    && evalue < THR) frequencies[Nfreq++] = 0.e0;
		//else frequencies[Nfreq++] = evalue;
		frequencies[Nfreq++] = evalue;
		double b[3];
		mproduct(vecq, qvec, b);
		for (int j=0; j<3; j++) b[j] /= PI2;
                //if (Debug) printf ("EV %d = %lf THz for Q= %7.4lf %7.4lf %7.4lf\n", 
		//    i, evalue, b[0], b[1], b[2]);
		    //i, evalue, vecQ[0], vecQ[1], vecQ[2]);
            }
	}

	/* solving Gamma point vibrational modes 
	  fexact - file name which the exact Q point phonon infor will be printed to
	  atomP - atomic infor in primitive unit cell 
	  natomP - number of atoms in the primitive unit cell
	  calcNA - key to consider D-D interaction
	  BCharg - BEC tensor
	  ILOTO -  for differing LO-TO modes
          IRed - for IR mode
	  IRaman - for Raman mode
	  RRvv - vibrational eigenvector
	*/
        double solveDMev(FILE *fexact, ATM **atomP, int natomP, int calcNA, double *BCharg,
		double *ILOTO, double *IRed, RAMAN *IRaman, double *RRvv) {
            double b[3];
            mproduct(vecq, qvec, b);
            for (int j=0; j<3; j++) b[j] /= PI2;
            fprintf (fexact, "   Qval= %9.6lf %9.6lf %9.6lf, q= %9.6lf %9.6lf %9.6lf\n\n",  b[0], b[1], b[2], vecq[0]/PI2, vecq[1]/PI2, vecq[2]/PI2);
	    fprintf(fexact, "                    ");
	    for (int i=0; i<natomP; i++) fprintf(fexact, "   %2s X %3d        %2s Y %3d        %2s Z %3d     ", 
		atomP[i]->sym, i, atomP[i]->sym, i, atomP[i]->sym, i);
	    fprintf(fexact, "\n                   ");
	    for (int i=0; i<3*natomP; i++) fprintf(fexact, " ---------------");
	    fprintf(fexact, "\n");

            for (int i=0; i<nTHz; i++) {
                for (int j=0; j<nTHz; j++) {
                    if (j>i) continue;
                    double R = real(dmat[i*nTHz+j]);
                    double I = imag(dmat[i*nTHz+j]);
                    gsl_complex z;
                    GSL_SET_COMPLEX (&z, R, I);
                    gsl_matrix_complex_set (m,i,j,z);
                }
            }

	/* solving eigen vectors */
            gsl_eigen_hermv (m, eval, evec, wv);
            gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	    double foundunstable = 0.0;
            for (int i=0;i<nTHz;i++) {
                double evalue = gsl_vector_get (eval, i);
                if (evalue<=0.e0) evalue = sqrt(-evalue)*toTHz;
                else {
		    evalue = -sqrt(evalue)*toTHz;
		    if (evalue<=-0.01) foundunstable = min(foundunstable,evalue);
		}
                frequencies[Nfreq++] = evalue;
                fprintf (fexact, " %4d %8.4lf THz ", i, evalue);
	    	double *Rv = new double[nTHz];
		for (int j=0; j<nTHz; j++) {
                    gsl_complex z = gsl_matrix_complex_get (evec,j,i);
		    fprintf (fexact, " %7.4lf %7.4lf",  GSL_REAL(z),  GSL_IMAG(z) );
		    RRvv[i*nTHz+j] = GSL_REAL(z);
		    Rv[j] = GSL_REAL(z);
		}
		fprintf (fexact, "\n");
		if (vecq[0]==0.e0 && vecq[1]==0.e0 && vecq[2]==0.e0 ) {
    			GammaPhonon GM(atomP, natomP, Rv, avec);
			if (calcNA) ILOTO[i] = GM.ByDipole(BCharg, qvec);
			else ILOTO[i] = 0.e0;
			IRed[i] = GM.BySymm();
			GM.ByRaman(IRaman+i);
		}
            }
	    return foundunstable;
        }

	/* solving partial phonon DOS to be stored in "pvdos" through eigen vectors
	  atomP - atomic infor in primitive unit cell 
	  natomP - number of atoms in the primitive unit cell
    	*/
        void solveDMPVDOS(ATM **atomP, int natomP) {

            for (int i=0; i<nTHz; i++) {
                for (int j=0; j<nTHz; j++) {
                    if (j>i) continue;
                    double R = real(dmat[i*nTHz+j]);
                    double I = imag(dmat[i*nTHz+j]);
                    gsl_complex z;
                    GSL_SET_COMPLEX (&z, R, I);
                    gsl_matrix_complex_set (m,i,j,z);
                }
            }

            gsl_eigen_hermv (m, eval, evec, wv);

	    double emin = xe[0];
	    double de = xe[1] - xe[0];

            for (int i=0;i<nTHz;i++) {
                double evalue = gsl_vector_get (eval, i);
                if (evalue>=0.e0) evalue = -sqrt(evalue)*toTHz;
                else evalue = sqrt(-evalue)*toTHz;
                frequencies[Nfreq++] = evalue;
		int j0 = (int)((evalue - emin)/de + 0.5e0);
		if (j0>=pvdosN) continue;
		if (j0<0) continue;

		for (int j=0; j<nTHz; j++) {
                    gsl_complex z = gsl_matrix_complex_get (evec,j,i);
		    double xR = GSL_REAL(z);
		    double xI = GSL_IMAG(z);
		    int N = j/3;
		    int Np = pvdosN*N+j0;
		    pvdos[Np] += (xR*xR + xI*xI);
		    xeN[j0] += 1;
		}
            }
        }

#endif

	/* not used, for IMSL eigenvalue solver */
#ifdef IMSL
	void solveDM() {
	    int N = nTHz*nTHz;
	    d_complex *m = new d_complex[N];
	    for (int i=0; i<N; i++) {
		double R = real(dmat[i]);
		double I = imag(dmat[i]);
		m[i] = imsl_zd_convert(R,I);
            }
	    double *eval = imsl_z_eig_herm (nTHz, m, 0);

            for (int i=0;i<nTHz;i++) {
                double evalue = eval[i];
                if (evalue>=0.e0) evalue = -sqrt(evalue)*toTHz;
                else evalue = sqrt(-evalue)*toTHz;
		if (vecq[0]==0.e0 && vecq[1]==0.e0 && vecq[2]==0.e0
		    && evalue < THR) frequencies[Nfreq++] = 0.e0;
		else frequencies[Nfreq++] = evalue;
		double b[3];
		mproduct(vecq, qvec, b);
		for (int j=0; j<3; j++) b[j] /= PI2;
                if (Debug) printf ("EV %d = %lf THz for Q= %7.4lf %7.4lf %7.4lf\n", 
		    i, evalue, b[0], b[1], b[2]);
            }
	}
#endif

	/* not used, for INTEL eigenvalue solver */
#ifdef INTEL
        void solveDM() {

            int N = nTHz*nTHz;
            MKL_Complex16 wkopt;
            MKL_Complex16* work;
	    double *eval = (double *)malloc((nTHz)*sizeof(double));
	    double *rwork = (double *)malloc((3*nTHz-2)*sizeof(double));
            MKL_Complex16 *matrix = (MKL_Complex16 *)malloc( N*sizeof(MKL_Complex16) );
            int info, lwork = 2*nTHz - 1;
            work = (MKL_Complex16 *)malloc( lwork*sizeof(MKL_Complex16) );

            for (int i=0; i<N; i++) {
                double R = real(dmat[i]);
                double I = imag(dmat[i]);
                matrix[i].real = R;
                matrix[i].imag = I;
		if (Debug) printf(" (%.4lf,%.4lf)\n", R, I);
            }

        /* Solve eigenproblem */ 
            zheev( "N", "U", &nTHz, matrix, &nTHz, eval, work, &lwork, rwork, &info );
            if( info != 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
            }

            for (int i=0;i<nTHz;i++) {
                double evalue = eval[i];
                if (evalue>=0.e0) evalue = -sqrt(evalue)*toTHz;
                else evalue = sqrt(-evalue)*toTHz;
                if (vecq[0]==0.e0 && vecq[1]==0.e0 && vecq[2]==0.e0
                    && evalue < THR) frequencies[Nfreq++] = 0.e0;
                else frequencies[Nfreq++] = evalue;
                double b[3];
                mproduct(vecq, qvec, b);
                for (int j=0; j<3; j++) b[j] /= PI2;
                if (Debug) printf ("EV %d = %lf THz for Q= %7.4lf %7.4lf %7.4lf\n",
                    i, evalue, b[0], b[1], b[2]);
            }
        }
#endif

	/* solving phonon frequencies in real space with the real space force constants matrix  Fij
	  eeigenfrequencies to be stored in "evalue"
	  natomS - number of atoms in the supercell
	*/
        static void SrealFij(double *Fij, int natomS, double *evalue) {
            int dimN = natomS*3;
            gsl_vector *eval = gsl_vector_alloc (dimN);
            gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc (dimN);
                        
            gsl_matrix *m = gsl_matrix_alloc (dimN, dimN);
            for (int i=0; i<dimN; i++) 
                for (int j=0; j<dimN; j++) gsl_matrix_set (m,i,j,Fij[i*dimN+j]);

            gsl_eigen_symm (m, eval, w);

            gsl_eigen_symm_free (w);
	    gsl_sort_vector(eval);

            for (int i=0;i<dimN;i++) {
                evalue[i] = gsl_vector_get (eval, i);
                if (evalue[i]>=0.e0) evalue[i] = -sqrt(evalue[i])*toTHz;
                else evalue[i] = sqrt(-evalue[i])*toTHz;
//printf ("eigen i=%d\n", i);
            }
            gsl_vector_free (eval);
        }       

        /* solving phonon frequencies & eigenvec in real space with the real space force constants matrix  Fij
          eeigenfrequencies to be stored in "evalue"
          natomS - number of atoms in the supercell
        */
        static void SrealFijV(double *Fij, int natomS, double *evalue, 
		double scell[3][3], double inv_scell[3][3], ATM **atomS, double Amode, const char *fh) {
            int dimN = natomS*3;
            gsl_vector *eval = gsl_vector_alloc (dimN);
            gsl_matrix *evec = gsl_matrix_alloc (dimN, dimN);
            gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (dimN);
            //gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc (dimN);

            gsl_matrix *m = gsl_matrix_alloc (dimN, dimN);
            for (int i=0; i<dimN; i++)
                for (int j=0; j<dimN; j++) gsl_matrix_set (m,i,j,Fij[i*dimN+j]);

            gsl_eigen_symmv (m, eval, evec, w);
            //gsl_eigen_symm (m, eval, w);

            gsl_eigen_symmv_free (w);
            //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
            gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
            //gsl_sort_vector(eval);

            for (int i=0;i<dimN;i++) {
                evalue[i] = gsl_vector_get (eval, i);
                if (evalue[i]>=0.e0) evalue[i] = -sqrt(evalue[i])*toTHz;
                else evalue[i] = sqrt(-evalue[i])*toTHz;
//printf ("i=%d, %lf\n", i, evalue[i]);
                if (evalue[i] < -0.001) {
        /* make the POSCAR file for further relax for the imagninary phonon mode*/
		  for (double bmode=-1.0; bmode<=1.0001; bmode += 2.0) {
        	    char file[256];
		    if (bmode<0.0) sprintf(file,"%s%d-",fh,i);
		    else sprintf(file,"%s%d+",fh,i);
        	    FILE *symmfile = fopen(file, "w");
        	    fprintf(symmfile, "Imaginary phonon No.%d = %.3lfTHz\n 1.00\n", i, evalue[i]);
        	    for (int ii=0; ii<3; ii++) println(symmfile, "%18.12f  ", scell[ii]);

        	    int *n = new int[natomS];
        	    int N = 1;
        	    n[0] = 1;

        	    fprintf(symmfile, " %4s", atomS[0]->sym);
        	    for (int ii=1; ii<natomS; ii++) {
            	        if (strcmp(atomS[ii-1]->sym, atomS[ii]->sym)) {
                	    fprintf(symmfile, " %4s", atomS[ii]->sym);
                	    n[N] = 1;
                	    N++;
            	        } else {
                	    n[N-1]++;
            	        }
        	    }
        	    fprintf(symmfile, "\n");

        	    for (int ii=0; ii<N; ii++) fprintf(symmfile, " %4d", n[ii]);
        	    fprintf(symmfile, "\nD\n");

		    double *vec = new double[dimN];
		    double umax = 0.0;
        	    for (int j=0; j<natomS; j++) {
            	        double jmass = 1./sqrt(getmass(atomS[j]->sym));
			//for (int ii=0; ii<3; ii++) vec[3*j+ii] = gsl_matrix_get(evec,i,j*3+ii)*jmass;
			for (int ii=0; ii<3; ii++) vec[3*j+ii] = gsl_matrix_get(evec,j*3+ii,i)*jmass;
			umax = max(umax, normal(vec+3*j));
		    }
        	    for (int j=0; j<natomS; j++) {
		        double v[3], c[3];
			for (int ii=0; ii<3; ii++) v[ii] = bmode*Amode/umax*vec[j*3+ii];
			mproduct(v, inv_scell, c);
                	fprintf(symmfile, " %15.10lf %15.10lf %15.10lf %2s %9.6lf %9.6lf %9.6lf %9.6lf\n",
                            atomS[j]->x+c[0], atomS[j]->y+c[1], atomS[j]->z+c[2], atomS[j]->sym, v[0], v[1], v[2], normal(v)/Amode);
        	    }	
        	    fclose(symmfile);
		  }
                }
            }
            gsl_vector_free (eval);
            gsl_matrix_free (evec);
        }

	/* printing the phonon frequencies at the specific Q points
	 print - key to decide if print more text comments */
	static void printQ(int print) {
            gsl_vector *eval = gsl_vector_alloc (Nfreq);
	    for (int i=0; i<Nfreq; i++) 
		gsl_vector_set (eval, i, frequencies[i]);
	    gsl_sort_vector(eval);

	    if (print>0) printf("Phonon frequencies:\n");
	    for (int i=0; i<Nfreq; i++) {
		frequencies[i] = gsl_vector_get (eval, i);
		if (print>0)printf("%20.3lf\n", frequencies[i]*1.e12);
	    }
	    if(print>0)printf("end\n");
	    gsl_vector_free(eval);
	}

	/* printing the phonon frequencies into temp file strem "cpq" for PDOS calculation */
	static void printPQ(FILE *cpq) {
            gsl_vector *eval = gsl_vector_alloc (nTHz);
	    for (int i=0; i<nTHz; i++) 
		gsl_vector_set (eval, i, frequencies[Nfreq-nTHz+i]);
	    gsl_sort_vector(eval);

	    for (int i=0; i<nTHz; i++) {
		fprintf(cpq, " %20.3lf", 1.e12*gsl_vector_get(eval, i));
	    }
	    fprintf(cpq,"\n");

	    gsl_vector_free(eval);
	}

	/* print/sort the phonon dispersions, in particular for the crossing of two dispersions
	  out - file stream
	  index - location of frequencies
	  Dir - Direction
	*/
        double printDIS(FILE *out, int index, double *Dir, double *fmin) { 
                gsl_vector *eval = gsl_vector_alloc (nTHz);
	        for (int i=0; i<nTHz; i++) 
		    gsl_vector_set (eval, i, frequencies[index*nTHz+i]);
	        gsl_sort_vector(eval);

            	for (int j=0; j<nTHz; j++) 
		    frequencies[index*nTHz+j] = gsl_vector_get (eval, j);

		if (index == 0) {
		    for (int i=0; i<nTHz; i++) 
			runa[i] = frequencies[index*nTHz+i];
		    for (int i=0; i<nTHz; i++) runc[i] = runa[i];
		} else if (index == 1) {
		    for (int i=0; i<nTHz; i++) 
			runb[i] = frequencies[index*nTHz+i];
		    for (int i=0; i<nTHz; i++) runc[i] = runb[i];
		} else {
		    for (int i=0; i<nTHz; i++) runa[i] = runb[i];
		    for (int i=0; i<nTHz; i++) runb[i] = runc[i];
		}

		if (index >1) {
		    int *kk = new int[nTHz];
		    for (int i=0; i<nTHz; i++) kk[i] = 0;
		    for (int i=0; i<nTHz; i++) {
			double fx = runb[i] + runb[i] - runa[i];
			double dd = 1.e72;
			int kused = 0;
			for (int k=0; k<nTHz; k++) {
			    if (kk[k]) continue;
			    double ff = frequencies[index*nTHz+k];
			    if (fabs(ff - fx) < dd) {
				dd = fabs(ff - fx);
				runc[i] = ff;
				kused = k;
			    }
			}
			kk[kused] = 1;
		    }
		    delete kk;
		}

		double tmp[3];
		mproduct (Dir, qvec, tmp);
		double qdis = normal(tmp); 

                fprintf(out, "%10.6lf %9.6lf %9.6lf %9.6lf ", qdis, 
		    vecQ[0], vecQ[1], vecQ[2]);

                for (int i=0; i<nTHz; i++) fprintf(out, " %9.6lf", runc[i]);
            	fprintf(out, "\n");

		double fmax = gsl_vector_get (eval, nTHz-1);
		*fmin = gsl_vector_get (eval, 0);
		gsl_vector_free(eval);
		return fmax;
        }

	/* calculate phonon DOS and print 
	  pdos - if not zero, also pring partial phonon DOS
	  atomP - atomic infor in primitive unit cell 
	  natomP - number of atoms in the primitive unit cell
	  nsym, nd - smoothening parameter
	  ne - number of meshes
	  DebCut - Instruct Yphon to fit the PDOS using the Debye expression in the form of C*f**2 for the phonon frequency range lower than the given frequency f. The typical value of f should be around 0.1 ~ 1.0 (in the unit of THz), depending your system and the used q mesh. One should always plot the calculated PDOS and inspect the region at around the frequency f for the smoothness of the curve to decide the proper value of f. This option is useful for the case if one wants to use the PDOS to calculate the heat capacity or Debye temperature at temperature lower than ~10 K, in particular, for superconductor etc.
	*/

	static void printDOS (int pdos, int natomP, ATM **atomP, 
	    int nsym=2, int nd=10, int ne=10001, double DebCut = -1.e0) {

	    double emin = 1.e36;
	    double emax = -1.e36;
	    for (int i=0; i<Nfreq; i++) {
		emin = min(emin, frequencies[i]);
		emax = max(emax, frequencies[i]);
	    }
	    ne = min(ne, Nfreq/10);
	    double de = (emax-emin)/(double)(ne-1);
	    if (emin > -de) {
		emin = 0.e0;
		de = (emax-emin)/(double)(ne-1);
	    }

	    double *e = new double[ne];
	    double *tdos = new double[ne];
	    double *dos = new double[ne];
	    double x = emin;

	    for (int i=0; i<ne; i++) {
		e[i] = x;
		x += de;
	    }

	    tdos[0] = 0.e0;
	    for (int i=0; i<Nfreq; i++) {
		int j = (int)((frequencies[i] - emin)/de);
		int tj = min(j+1, ne-1);
		tj = max(tj, 0);
		tdos[tj] = (double)(i);
		//tdos[tj] = (double)(i-2);
		tdos[tj] = max(0.e0, tdos[tj]);
	    }

            for (int i=1; i<ne; i++)
                tdos[i] = max(tdos[i-1], tdos[i]);

	    double fac = (double)(nTHz)/(double)Nfreq;
	    int isym = nsym*nsym*nd;
	    dos[0] = 0.e0;
	    dos[ne-1] = 0.e0;
	    for (int i=1; i<ne-1; i++) {
		int ssym = min(isym,i);
		ssym = min(ssym, ne-i-1);
		double tmp =0.e0;
		for (int j=0; j<=ssym; j++) 
		    tmp = tmp + tdos[i+j] - tdos[i-ssym+j];
		dos[i] = tmp/de/(double)(ssym*(ssym+1))*fac;
	    }

	    double DebDOS = 0.e0;
	    for (int i=0; i<ne; i++) {
		if (e[i]>= DebCut) {
		     DebDOS = dos[i]/DebCut/DebCut;
		     break;
		}
	    }

	/* phonon DOS */

	    FILE *fdos = fopen("vdos.out", "w");
//printf("DebDOS=%lf, DebCut=%lf\n", DebDOS, DebCut);
	    for (int i=0; i<ne; i++) {
		if (DebCut>0.e0 && e[i]/de >-50. && e[i] <= DebCut) 
		    fprintf(fdos, "%.12lg %.12lg\n", e[i]*1.e12, e[i]*e[i]*DebDOS*1.e-12);
		else 
		    fprintf(fdos, "%.12lg %.12lg\n", e[i]*1.e12, dos[i]*1.e-12);
	    }
	    fclose(fdos);

	    if (pdos==0) return;
	/* further print partial phonon DOS */

	    double *xdos = new double[ne];
	    int ngap = 6;
	    isym = ngap*ngap*nd;
	    xdos[0] = 0.e0;
	    xdos[ne-1] = 0.e0;
	    for (int i=1; i<ne-1; i++) {
		int ssym = min(isym,i);
		ssym = min(ssym, ne-i-1);
		double tmp =0.e0;
		for (int j=0; j<=ssym; j++) 
		    tmp = tmp + tdos[i+j] - tdos[i-ssym+j];
		xdos[i] = tmp/de/(double)(ssym*(ssym+1));
	    }

	    ngap = 0;
	    for (int i=ne/20; i<ne-1; i++) {
		if ( e[i] < 0.e0) continue;
		ngap = i+1;
		if ( xdos[i+1] < xdos[i] ) break;
	    }
	    double egapmin = e[ngap];

	    printf("\n Searching for phonon gap from f = %g THz\n", egapmin);

            emin = xe[0];
            de = xe[1] - xe[0];

	    for (int i=0; i<pvdosN; i++)
		if (xeN[i]>0) for (int j=0; j<natomP; j++) 
		    pvdos[pvdosN*j+i] /= xeN[i]/nTHz;

	    double *tpv = new double[pvdosN];
	    int Ns = 99; 
	    for (int j=0; j<natomP; j++) {
	        for (int i=0; i<pvdosN; i++) tpv[i] = pvdos[pvdosN*j+i];
	        for (int i=0; i<pvdosN; i++) {
		    int I0 = i;
		    int N = 0;
		    for (;I0>=0;I0--) {
			if (xeN[I0]>0) N++;
			if (N>=Ns) break;
			//if (xeN[I0]==0 && xe[I0] > egapmin) break;
		    }
		    int I1 = i;
		    N = 0;
		    for (;I1<pvdosN;I1++) {
			if (xeN[I1]>0) N++;
			if (N>=Ns) break;
			//if (xeN[I1]==0 && xe[I1] > egapmin) break;
		    }
		    double tmp =0.e0;
		    N = 0;
		    for (int k=I0; k<=I1; k++) {
			if (xeN[k]==0) continue;
			N += xeN[k];
			tmp += tpv[k]*xeN[k];
		    }
		    if (N>0) pvdos[pvdosN*j+i] = tmp/N;
		    else pvdos[pvdosN*j+i] = 0.e0;
		}
	    }

	/* GPDOS calculations */
	    double *S = new double[natomP];
	    double *pNdos = new double[ne];
	    double *xNdos = new double[ne];
	    printf ("\nUsed Neutron scattering cross section devided by atomic Mass\n\n");
	    for (int i=0; i<natomP; i++) {
		S[i] = getSdm(atomP[i]->sym);
		printf ("    sigma/mass of %2s = %.6lf\n", atomP[i]->sym, S[i]);
	    }

	    for (int i=0; i<ne; i++) {
		double evalue = e[i];
                int j0 = (int)((evalue - emin)/de);
                if (j0+1>=pvdosN) continue;
                if (j0<0) continue;
                double f1 = (evalue - xe[j0])/de;
                double f0 = (xe[j0+1] - evalue)/de;
		xNdos[i] = 1.e0;
		if (evalue > egapmin) 
		   if (xeN[j0]==0 && xeN[j0+1]==0) xNdos[i] = 0.e0;
		double tmp = 0.e0;
		for (int j=0; j<natomP; j++) {
		    int N = pvdosN*j+j0;
		    double v0 = pvdos[N];
		    double v1 = pvdos[N+1];
		    double val = v0*f0+v1*f1;
		    tmp += val*S[j];
		}
		pNdos[i] = tmp;
	    }

	    ngap = ne/200;
	    for (int i=0; i<ne; i++) {
		if (e[i] < egapmin) continue;
		if (xNdos[i]>0.e0) continue;
		int I0 = i;
		for (;I0>0;I0--) if (xNdos[I0] > 0.e0) break;
		int I1 = i;
		for (;I1<ne;I1++) if (xNdos[I1] > 0.e0) break;
		if (I1-I0 < ngap) xNdos[i] = 1.0;
	    }

            double tmpd = 0.e0;
            for (int i=0; i<ne; i++) tmpd +=  xNdos[i]*dos[i];
            tmpd *= e[1]-e[0];
            //for (int i=0; i<ne; i++) dos[i] = xNdos[i]*dos[i]/tmpd;

            double tmpy = 0.e0;
	    double *pYdos = new double[ne];
	    int N = 0;
            for (int i=0; i<ne; i++) 
		if (dos[i] > 0.e0 && e[i] > 0.e0 ) {
		    tmpy += pNdos[i];
		    N++;
		}
            for (int i=0; i<ne; i++) {
		if (e[i] > 0.e0 ) pYdos[i] = pNdos[i]/tmpy*N;
		else pYdos[i] = 0.e0;
	    }

            double tmpx = 0.e0;
            for (int i=0; i<ne; i++) if (e[i] > 0.e0 ) tmpx += pNdos[i]*xNdos[i]*dos[i]/tmpd;
            tmpx *= e[1]-e[0];
            for (int i=0; i<ne; i++) {
		if (e[i] > 0.e0 ) pNdos[i] = pNdos[i]/tmpx;
		else pNdos[i] = 0.e0;
	    }

	    fdos = fopen("pvdos.out", "w");
	    for (int i=0; i<ne; i++) {
		double evalue = e[i];
                int j0 = (int)((evalue - emin)/de);
                if (j0+1>=pvdosN) continue;
                if (j0<0) continue;
                double f1 = (evalue - xe[j0])/de;
                double f0 = (xe[j0+1] - evalue)/de;
                if (DebCut>0.e0 && e[i]/de >-50. && e[i] <= DebCut)
		    fprintf(fdos, "%g %g %g %g", e[i],  e[i]*e[i]*DebDOS, pNdos[i], pYdos[i]);
                else
		    fprintf(fdos, "%g %g %g %g", e[i], dos[i], pNdos[i], pYdos[i]);

		//fprintf(fdos, "%g %g %g %g", e[i], dos[i], pNdos[i], pYdos[i]);
		for (int j=0; j<natomP; j++) {
		    int N = pvdosN*j+j0;
		    double v0 = pvdos[N];
		    double v1 = pvdos[N+1];
		    double val = v0*f0+v1*f1;
		    fprintf(fdos, " %g", val);
		}
		fprintf(fdos, "\n");
	    }
	    fclose(fdos);
	}
};

/* class for Rotation matrix for restore symmetry */

class RMATRAX {
    private:
        double (*m)[3][3];
        double (*RR)[3][3];
        double (*R_R)[3][3];
        double (*t)[3];
        int Nrot;
        double vectorabc[3][3];
	double Shift[3];
	double inv_scell[3][3], inv_vectorabc[3][3];
	int natomS;
	ATM **atomS;

    public:
	/* read in rotation operations infor from file "rfile" */
        RMATRAX (char *rfile) { //rfile contains the symmetry operations, the "Rotation.sym" file by pos2s
	    FILE *cfp = fopen(rfile, "r");
            if (!cfp) {
                printf("\n********FETAL ERROR, CANNOT FIND file: %s\n\n", rfile);
                exit(1);
            }
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++) foo = fscanf(cfp, "%lf", &vectorabc[i][j]);
            for (int j=0; j<3; j++) foo = fscanf(cfp, "%lf", &Shift[j]);
	    foo = fscanf(cfp, "%d", &Nrot);

	    inv_m(vectorabc, inv_vectorabc);

	    double tmp[3][3];

            m = new double[Nrot][3][3];
            RR = new double[Nrot][3][3];
            R_R = new double[Nrot][3][3];
            t = new double[Nrot][3];

            for (int n=0; n<Nrot; n++) {
                for (int i=0; i<3; i++) 
                    for (int j=0; j<3; j++) foo = fscanf(cfp, " %lf", &m[n][i][j]);
                for (int j=0; j<3; j++) foo = fscanf(cfp, "%lf", &t[n][j]);
		mproduct(inv_vectorabc, m[n], tmp);
		mproduct(tmp, vectorabc, RR[n]);
		inv_m(RR[n], R_R[n]);
            }
	    fclose(cfp);
	}

	/* restore the Rotational symmetry of the force constant matirx 
	  Fij_Orig - original force constatns calculated by VASP.5
	  scell - supercell lattice vectors
	  atomS - atimic infor (positions and atomic symbols etc) within supercell
	  natomS - number of atoms in the supercell
	*/
        void ResRotSym(double *Fij_Orig, int _natomS, ATM **_atomS, double scell[3][3], int kvec) {
	    natomS = _natomS;
	    atomS = _atomS;
	    inv_m(scell, inv_scell);

	    int *iiijjj = new int[natomS*Nrot];
	    for (int i=0; i<natomS; i++) {
		double *pi = atomS[i]->getxyz();
		for (int r=0; r<Nrot; r++) iiijjj[i*Nrot+r] = atomRotMap(pi, atomS[i]->sym, r, kvec);
	    }

	    double fij[3][3], nfij[3][3];
	    double *nFij = new double[natomS*natomS*9];
	    for (int i=0; i<natomS*natomS*9; i++) nFij[i] = 0.e0;
	    for (int i=0; i<natomS; i++) {
		
		for (int j=0; j<natomS; j++) {
		    for (int r=0; r<Nrot; r++) {
			int ii = iiijjj[i*Nrot+r];
			int jj = iiijjj[j*Nrot+r];
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) fij[x][y] = Fij_Orig[(i*3+x)*natomS*3+j*3+y];
			Add(fij, nfij, r);
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) nFij[(ii*3+x)*natomS*3+jj*3+y] += nfij[x][y];
		    }
		}
	    }

	    for (int i=0; i<natomS*natomS*9; i++) Fij_Orig[i] = nFij[i]/(double)Nrot;
	    delete[] nFij;
	    delete[] iiijjj;
	}

	/* restore the Translational symmetry of the force constant matirx 
	  Fij_Orig - original force constatns calculated by VASP.5
	  avec - lattice vectors of PUC
	  scell - supercell lattice vectors
	  atomS - atimic infor (positions and atomic symbols etc) within supercell
	  natomS - number of atoms in the supercell
	  cellA - pointer to index of PUC in the supercell, index of atom in the PUC, index of atom in the supercell
	  cellI - pointer to index containing three integers index PUC positions in the supercell
	  Ntrans - number of translational operation
	*/
        void ResTransSym(double *Fij_Orig, int _natomS, ATM **_atomS, 
		double avec[3][3], double scell[3][3], CINDEX **cellI, CINDEX **cellA, int Ntrans) {

	    int *iiijjj = new int[natomS*Ntrans];
	    for (int i=0; i<natomS; i++) {
		double *pi = atomS[i]->getxyz();
		int M0 = cellA[0]->get_atomP();
	    	int Nt = 0;
		for (int r=0; r<natomS; r++) {
	    		int M = cellA[r]->get_atomP();
	    		if (M!=M0) continue;
			int iM = cellA[r]->get_cell();
            		//int iM = cellA[r]->get_atomP();
			int xM = cellI[iM]->get_cx();
			int yM = cellI[iM]->get_cy();
			int zM = cellI[iM]->get_cz();
			iiijjj[i*Ntrans+Nt] = atomTransMap(pi,  atomS[i]->sym, xM, yM, zM, avec, scell);
			Nt++;
		}
	    }

	    double *nFij = new double[natomS*natomS*9];
	    for (int i=0; i<natomS*natomS*9; i++) nFij[i] = 0.e0;
	    for (int i=0; i<natomS; i++) {
		for (int j=0; j<natomS; j++) {
		    for (int r=0; r<Ntrans; r++) {
			int ii = iiijjj[i*Ntrans+r];
			int jj = iiijjj[j*Ntrans+r];
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) 
				nFij[(ii*3+x)*natomS*3+jj*3+y] += Fij_Orig[(i*3+x)*natomS*3+j*3+y];
		    }
		}
	    }

	    for (int i=0; i<natomS*natomS*9; i++) Fij_Orig[i] = nFij[i]/(double)Ntrans;
	    delete[] nFij;
	    delete[] iiijjj;
	}

	/* accumulate the fij 
	   R - rotation matrix
	   R_R - inverse of R
	   Fij - force constat
	   r - index of rotational operation
	   nFij - result of Fij after ratation transformation, nFij = R**-2*Fij*R
	*/
	void Add(double Fij[3][3], double nFij[3][3], int r) {
	    double x[3][3];
	    mproduct(R_R[r], Fij, x);
	    mproduct(x, RR[r], nFij);
	}

	/* restoring the symmetry of the BEC tensor
	  dielec - dielectric constant tensor
	  BCharg - BEC tensor
	  atomS - atimic infor (positions and atomic symbols etc) within PUC
	  natomS - number of atoms in the PUC
	  scell - lat. vec. of PUC
	*/
        void ResRotSymBorn(double *dielec, double *BCharg, int _natomS, ATM **_atomS, double scell[3][3]) {
	    natomS = _natomS;
	    atomS = _atomS;
	    inv_m(scell, inv_scell);
	    double fij[3][3], nfij[3][3];
	    double *nFij = new double[natomS*natomS*9];
				for (int i=0; i<9; i++) nFij[i] = 0.e0;
		for (int r=0; r<Nrot; r++) {
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) fij[x][y] = dielec[x*3+y];
			Add(fij, nfij, r);
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) nFij[x*3+y] += nfij[x][y];
		}

	    for (int i=0; i<9; i++) dielec[i] = nFij[i]/(double)Nrot;

		printf ("\n\n Symmetrimized dielectric tensor \n");
		double *p = dielec;
		//printf (" %10.5lf %10.5lf %10.5lf\n", *p++, *p++, *p++);
		printf (" %10.5lf %10.5lf %10.5lf\n", p[0], p[1], p[2]);
		printf (" %10.5lf %10.5lf %10.5lf\n", p[3], p[4], p[5]);
		printf (" %10.5lf %10.5lf %10.5lf\n", p[6], p[7], p[8]);

	    for (int i=0; i<natomS*natomS*9; i++) nFij[i] = 0.e0;
	    for (int i=0; i<natomS; i++) {
		double *pi = atomS[i]->getxyz();
		    for (int r=0; r<Nrot; r++) {
	    if (LatDYN::Debug) printf(" atom %lf %lf %lf ", pi[0], pi[1], pi[2]);
			int ii = atomRotMap(pi,  atomS[i]->sym, r);
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) fij[x][y] = BCharg[i*9+x*3+y];
			Add(fij, nfij, r);
			for (int x=0; x<3; x++)
			    for (int y=0; y<3; y++) nFij[ii*9+x*3+y] += nfij[x][y];
		    }
	    }

	    for (int i=0; i<natomS*9; i++) BCharg[i] = nFij[i]/(double)Nrot;
	    for (int i=0; i<natomS; i++) {
	        p = BCharg+i*9;
		printf ("\n\n Symmetrimized BEC for ion %d\n", i+1);
		printf (" %10.5lf %10.5lf %10.5lf\n", p[0], p[1], p[2]);
		printf (" %10.5lf %10.5lf %10.5lf\n", p[3], p[4], p[5]);
		printf (" %10.5lf %10.5lf %10.5lf\n", p[6], p[7], p[8]);

	    	double tmp[3], tmp3[3];
		printf ("\n    Projecttion to PUC lat axis\n");
	    	for (int ii=0; ii<3; ii++) {
		    for (int j=0; j<3; j++) tmp3[j] = scell[ii][j];
		    Anormal (tmp3);
		    mproductR(p, tmp3, tmp);
		    printf ( " %10.5lf", dotproduct( tmp3, tmp));
		}
		printf ("\n");
	    }
	}

	/* accumulate the BEC tensor for averaging
	   R - rotation matrix
	   Fij - BEC tensor
	   r - index of rotational operation
	   nFij - result of Fij after rotation transformation, nFij = Fij*R
	*/
	void AddBorn(double Fij[3][3], double nFij[3][3], int r) {
	    mproduct(Fij, RR[r], nFij);
	}

	/* rotating atoms
	  _xyz - atomic position
	   r - index of rotational operation
	*/

	int atomRotMap(double *_xyz, char *sym, int r, int kvec=0) {
	    double xyz[3], pxyz[3], vxyz[3], sxyz[3], scell[3][3];
	    inv_m(inv_scell, scell);
	    mproduct (_xyz, scell, pxyz);
	    for (int i=0; i<3; i++) pxyz[i] += Shift[i];
	    mproduct (pxyz, inv_vectorabc, xyz);
	    mproduct(xyz, m[r], pxyz);
	    //mproduct(m[r], xyz, pxyz);
	    for (int i=0; i<3; i++) pxyz[i] += t[r][i];
	    //fprintf(stderr, "\n\n atom %lf %lf %lf for r=%d\n\n", _xyz[0], _xyz[1], _xyz[2], r);
	    //fprintf(stderr, "\n\n atom %lf %lf %lf for r=%d\n\n", pxyz[0], pxyz[1], pxyz[2], r);

	    mproduct(pxyz, vectorabc, vxyz);
	    for (int i=0; i<3; i++) vxyz[i] -= Shift[i];
	    //fprintf(stderr, "\n\n mapping atom %lf %lf %lf\n\n", vxyz[0], vxyz[1], vxyz[2]);
	    mproduct(vxyz, inv_scell, sxyz);

	    for (int i=0; i<3; i++) sxyz[i] = normal(sxyz[i]);
	    ATM atm(sxyz[0], sxyz[1], sxyz[2], sym);
	    if (LatDYN::Debug) printf(" mapping atom %lf %lf %lf  %s\n", _xyz[0], _xyz[1], _xyz[2], sym);
	    for (int i=0; i<natomS; i++) {
		int ii = 0;
		if (kvec==2) ii = atm.equalRot(*atomS[i]);
		else ii = (atm==*atomS[i]);
		//if (atm==*atomS[i]) {
		if (ii) {
	            if (LatDYN::Debug) printf(" mapped into atom %lf %lf %lf  %s\n", sxyz[0], sxyz[1], sxyz[2], atomS[i]->sym);
		    return i;
		}
	    }
	    printf("\n\n********FETAL ERROR, CANNOT mapping atom %lf %lf %lf %s for r=%d\n\n", _xyz[0], _xyz[1], _xyz[2], sym, r);
	    printf("\n\n********FETAL ERROR, CANNOT mapping atom %lf %lf %lf\n\n", pxyz[0], pxyz[1], pxyz[2]);
	    printf("\n\n********FETAL ERROR, CANNOT mapping atom %lf %lf %lf\n\n", sxyz[0], sxyz[1], sxyz[2]);
	    exit(1);
	    return -1;
	}

	/* translating atoms
	  _xyz - atomic position
	   x, y, z - translational operation
	   avec - lat vec of PUC
	   scell - lat vec of supercell 
	*/

	int atomTransMap(double *_xyz, char *sym, int x, int y, int z, double avec[3][3], double scell[3][3]) {
	    double xyz[3], pxyz[3], sxyz[3];
	    xyz[0] = (double)x;
	    xyz[1] = (double)y;
	    xyz[2] = (double)z;
	    mproduct (xyz, avec, pxyz);
	    mproduct (_xyz, scell, sxyz);
	    for (int i=0; i<3; i++) sxyz[i] += pxyz[i];
	    mproduct(sxyz, inv_scell, pxyz);

	    for (int i=0; i<3; i++) pxyz[i] = normal(pxyz[i]);
	    ATM atm(pxyz[0], pxyz[1], pxyz[2], sym);
	    for (int i=0; i<natomS; i++) 
		if (atm==*atomS[i]) {
	        if (LatDYN::Debug) printf(" mapped into atom %lf %lf %lf\n", sxyz[0], sxyz[1], sxyz[2]);
		return i;
	    }
	    printf("\n\n********FETAL ERROR, CANNOT mapping atom %lf %lf %lf\n\n", _xyz[0], _xyz[1], _xyz[2]);
	    printf("\n\n********FETAL ERROR, CANNOT mapping atom %lf %lf %lf\n\n", pxyz[0], pxyz[1], pxyz[2]);
	    printf("\n\n********FETAL ERROR, CANNOT mapping atom %lf %lf %lf\n\n", sxyz[0], sxyz[1], sxyz[2]);
	    exit(1);
	    return -1;
	}
};


/* class for vibrational mode analysis at Gamma point */
class PhononG {

    private:

	int Raman;
	int translational;
	int Nmode;
	double **coef;
	double *evalue;
	double *evalueB;
	double *dmat;
	dcmplx *dmatC;
	int nTHz;
	double *symmA;

    public:
	char *irrep;
	int IR;
	double *dirP;
	double *dirCart;
	int end_of_file;

	/* read in Gamma point infor from the symmetry file 
	  _nTHz = number of phonon frequencies at Gamma pionts (dispersions)
	*/
	PhononG (FILE *file, int _nTHz) {
	    end_of_file = 0;

	    char line[80];
	    foo = fscanf (file, "%s %d", line, &Nmode); //read in the Irrep and number of modes (Nmode)
	    if (feof(file)) {
		end_of_file =1;
		return;
	    }
	    irrep = new char(strlen(line) + 1);
	    strcpy(irrep, line);

	    foo = fscanf (file, "%s %d", line, &IR); 
	    foo = fscanf (file, "%s %d", line, &Raman); 
	    foo = fscanf (file, "%s %d", line, &translational); 
	    nTHz = _nTHz;

	    dirP = 0;
	    evalue = new double[Nmode];
	    evalueB = new double[Nmode];
	    dmat = new double[Nmode*Nmode];
	    dmatC = new dcmplx[Nmode*Nmode];
	    coef = new double*[Nmode];
	    symmA = new double[nTHz*Nmode];

	    for (int i=0; i<Nmode; i++) {
		coef[i] = new double[nTHz];
		double *p = coef[i];
		for (int j=0; j<nTHz; j++) foo = fscanf(file, "%lf", p++);

		double tmp = 0.e0;
		p = coef[i];
		for (int j=0; j<nTHz; j++) tmp += p[j]*p[j];

		tmp = 1.e0/sqrt(tmp);
		for (int j=0; j<nTHz; j++) p[j] *= tmp;
	    }
	}

	/*  find the overlapping between the eigenvec RRvv with the Nmode of individual mode*/
    	void Analysis (double *RRvv) {
	    for (int i=0; i<nTHz; i++) {
		for (int r=0; r<Nmode; r++) {
		    double *pr = coef[r];
		    double tmp =0.e0;
		    for (int j=0; j<nTHz; j++) tmp += RRvv[i*nTHz+j]*pr[j];
		    symmA[i*Nmode+r] = tmp;
		}
	    }
	}

	/* printing the results of irr representation 
	  EEvv - eigenvector
	  EEss - eigenvalue
	  IRRep - Irrep with numbering
	*/
    	void print_symmA (double *EEvv, double *EEss, char **IRRep) {
	    for (int i=0; i<nTHz; i++) {
		double tmp =0.e0;
		for (int r=0; r<Nmode; r++) tmp += fabs(symmA[i*Nmode+r]);
	    	if (tmp>THR) {
		    if (tmp < EEss[i]) continue;
		    EEss[i] = tmp;
		    IRRep[i] = new char[strlen(irrep)+1];
		    strcpy(IRRep[i], irrep);
		    //printf (" State %4d with energy %8.4lf, is Projected into %s with weight %8.4lf\n", i, EEvv[i], IRRep[i], tmp);
		}
	    }
	}

	/* set up the electric field direction */
    	void SetDirP () {
	    dirCart = new double[3];
	    double inv_qvec[3][3];
	    inv_m(LatDYN::qvec, inv_qvec);
            for (int j=0; j<3; j++) dirCart[j] = 0.e0;
	    for (int i=0; i<Nmode; i++) {
		double *pi = coef[i];
                for (int N=0; N<nTHz/3; N++) {
                    for (int j=0; j<3; j++) dirCart[j] += pi[3*N+j];
		}
            }  

            if (normal (dirCart) > THR) {
		dirP = new double[3];
		Anormal(dirCart); 
		mproduct(dirCart, inv_qvec, dirP);
		Anormal(dirP); 
	    }
        }

	/* making the dynamical matrix 
	  Fij - force constatn matrix
	*/
        void makeDM(double *Fij) {
	    for (int i=0; i<Nmode; i++) {
		double *pi = coef[i];
		for (int j=0; j<Nmode; j++) {
		    double *pj = coef[j];
		    double tmp = 0.e0;
		    for (int k=0; k<nTHz; k++) {
			for (int l=0; l<nTHz; l++) {
			    if (l<=k) tmp += pi[k]*Fij[k*nTHz+l]*pj[l];
			    else tmp += pi[k]*Fij[l*nTHz+k]*pj[l];
			}
		    }
		    dmat[i*Nmode+j] = tmp;
		}
	    }
	}

	/* return the Text representations for different modes */
	const char *getmode() {
            if (IR) return "ir_active"; 
            else if (Raman) return "raman_active";
            else return "silent_mode";
	}

	/* real space sover at Gamma point 
	  Tij - real space force constant
	  Kcell - number of PUC in the supercell
	  Dij - D-D interaction matrix induced by vibration
	*/
        void SrealFij(double *Tij, int Kcell=-1, dcmplx *Dij=0) {
	    double *Fij = new double[nTHz*nTHz];
	    if (Dij) for (int i=0; i<nTHz*nTHz; i++) Fij[i] = Tij[i] + ((double)Kcell)*real(Dij[i]);
	    else for (int i=0; i<nTHz*nTHz; i++) Fij[i] = Tij[i];
            makeDM(Fij);

            gsl_vector *eval = gsl_vector_alloc (Nmode);
            gsl_matrix *evec = gsl_matrix_alloc (Nmode, Nmode);
            gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (Nmode);
            
            gsl_matrix *m = gsl_matrix_alloc (Nmode, Nmode);
            for (int i=0; i<Nmode; i++)
                for (int j=0; j<Nmode; j++) gsl_matrix_set (m,i,j,dmat[i*Nmode+j]);
                
            gsl_eigen_symmv (m, eval, evec, w);
            //gsl_eigen_symm (m, eval, w);
                
            gsl_eigen_symmv_free (w);
            //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
            gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
            //gsl_sort_vector(eval);
                
	    if (Kcell==0 || Dij) {
                printf("\n%3d %4s Modes of %s ", Nmode, irrep, getmode());
	        if (dirP) printf(": P= %6.3lf,%6.3lf,%6.3lf (%6.3lf,%6.3lf,%6.3lf) ", 
		dirCart[0], dirCart[1], dirCart[2],
	 	dirP[0], dirP[1], dirP[2]);
                if (translational) printf(" with one tranlational mode ");
                printf("\n");
	    }
        
            for (int i=0;i<Nmode;i++) {
                double tmp = gsl_vector_get (eval, i);
		if (tmp >=0) tmp = -sqrt(tmp)*toTHz;
		else tmp = sqrt(-tmp)*toTHz;
		if (Dij) evalueB[i] = tmp;
		else evalue[i] = tmp;

		if (LatDYN::Debug) {
                    printf (", Mode=");
                    for (int j=0; j<Nmode; j++) printf (" %7.4lf",  gsl_matrix_get (evec,j,i));
		}

		if (Kcell==0) {
		    printf("%3d %4s %8.4lf  (%8.2lf)", i, irrep, evalue[i],  evalue[i]*1.e10/LIGHTSPEED);
                    printf ("\n");
                } else if (Dij) {
		    printf("%3d %4s %8.4lf  %8.4lf  (%8.2lf  %8.2lf)", i, irrep, evalue[i], evalueB[i],
			evalue[i]*1.e10/LIGHTSPEED, evalueB[i]*1.e10/LIGHTSPEED);
                    printf ("\n");
		}
            }       
        }

	/* complex space solver at Gamma point
	  Tij - real space force constant
	  Kcell - number of PUC in the supercell
	  Dij - D-D interaction matrix induced by vibration
	  atomP - atomic infor at the primitive unit cell (PUC)
	  natomP - number of atoms in the PUC
	  BCharg - BEC tensor
	  tZmu - for the mode dependent charge - underdevelopping
	*/
        void SrealDij(double *Tij, int Kcell, dcmplx *Dij, ATM **atomP, int natomP, double *BCharg, double *tZmu) {
			double *Fij = new double[nTHz*nTHz];
			for (int i=0; i<nTHz*nTHz; i++) Fij[i] = Tij[i] + ((double)Kcell)*real(Dij[i]);
            makeDM(Fij);

            gsl_vector *eval = gsl_vector_alloc (Nmode);
            gsl_matrix *evec = gsl_matrix_alloc (Nmode, Nmode);
            gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (Nmode);
            
            gsl_matrix *m = gsl_matrix_alloc (Nmode, Nmode);
            for (int i=0; i<Nmode; i++)
                for (int j=0; j<Nmode; j++) gsl_matrix_set (m,i,j,dmat[i*Nmode+j]);
                
            gsl_eigen_symmv (m, eval, evec, w);
            //gsl_eigen_symm (m, eval, w);
                
            gsl_eigen_symmv_free (w);
            //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
            gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
            //gsl_sort_vector(eval);

			double *Rv = new double[nTHz];

            printf("\n%3d %4s Modes of %s ", Nmode, irrep, getmode());
	        if (dirP) printf(": P= %6.3lf,%6.3lf,%6.3lf (%6.3lf,%6.3lf,%6.3lf) ", 
				dirCart[0], dirCart[1], dirCart[2],
	 			dirP[0], dirP[1], dirP[2]);
            if (translational) printf(" with one tranlational mode ");
            printf("\n No irrep        THz                 (cm-1)           Z*(x)     Z*(y)     Zz(z)\n");
       
            for (int i=0;i<Nmode;i++) {
				double tmp = gsl_vector_get (eval, i);
				if (tmp >=0) tmp = -sqrt(tmp)*toTHz;
				else tmp = sqrt(-tmp)*toTHz;
				evalueB[i] = tmp;

				if (LatDYN::Debug) {
                    printf (", Mode=");
                    for (int j=0; j<Nmode; j++) printf (" %7.4lf",  gsl_matrix_get (evec,j,i));
				}

				for (int j=0; j<nTHz; j++)	Rv[j] = 0.e0;

				for (int j=0; j<Nmode; j++) {
					double z = gsl_matrix_get (evec,j,i);
				    double *pj = coef[j];
					for (int l=0; l<nTHz; l++) Rv[l] += z*pj[l];
				}

 				double sum = 0.e0;
				for (int j=0; j<nTHz; j++) {
					sum += Rv[j]*Rv[j];
					Rv[j] /= sqrt(getmass(atomP[j/3]->sym));
				}

				for (int j=0; j<nTHz; j++) Rv[j] /= sqrt(sum);

				double Zmu[3];
				for (int j=0; j<3; j++) Zmu[j] = 0.e0;

				for (int j=0; j<nTHz; j++) {
					int jA = j/3;
					int jy = j%3;
					double *pB = BCharg+jA*9;
					Zmu[0] += Rv[j]*pB[jy];
					Zmu[1] += Rv[j]*pB[3+jy];
					Zmu[2] += Rv[j]*pB[6+jy];
				}

				printf("%3d %4s %8.4lf  %8.4lf  (%8.2lf  %8.2lf)  %8.4lf  %8.4lf  %8.4lf", i, irrep, evalue[i], evalueB[i],
				evalue[i]*1.e10/LIGHTSPEED, evalueB[i]*1.e10/LIGHTSPEED, Zmu[0], Zmu[1], Zmu[2]);
				printf ("\n");

				if (evalueB[i]>=0.01e0) for (int j=0; j<3; j++) tZmu[j] += Zmu[j]/evalueB[i];
            }       
        } 
};

/* Deriving elastic constants from Fij (GPa). Not accurate, under developping
	the following codes will be deleted due to low accuracy in calculating elastic constants from
	force constants
*/
void Cij(double *Gij, int natomP, double avec[3][3], ATM **atomP,
	    double *stress) {
	double Cijkl[3][3][3][3], Xijkl[3][3][3][3], Cij[6][6];
	double ss[3][3];
	int i,j,k,l;

	for (i=0; i<3; i++) ss[i][i] = -stress[i]*0.1e0;
	ss[0][1] =  stress[3]*0.1e0;
	ss[1][0] =  stress[3]*0.1e0;
	ss[1][2] =  stress[4]*0.1e0;
	ss[2][1] =  stress[4]*0.1e0;
	ss[0][2] =  stress[5]*0.1e0;
	ss[2][0] =  stress[5]*0.1e0;

        CijRDIS **disR = (CijRDIS **) malloc ( (size_t) (natomP*sizeof(CijRDIS *)) );
        for (int nx=0; nx<natomP; nx++)
            disR[nx] = new CijRDIS(atomP, natomP, avec, nx);

	double V = volume(avec);
	for (i=0; i<3; i++) {
	    for (j=0; j<3; j++) {
		for (k=0; k<3; k++) {
		    for (l=0; l<3; l++) {
			double tmp = 0.e0;
			for (int m=0; m<natomP; m++) {
			for (int n=0; n<natomP; n++) {
			    int MX = disR[m]->get_NxS(n);
			    double fac = 1.e0/MX;
			    for (int mx=0; mx<MX; mx++) {
		    	    	double *aa = disR[m]->get_IxS(n, mx);
			    	int ii = m*3 + i;
			    	int jj = n*3 + j;
			    	int ij = ii*natomP*3+jj;
			    	tmp += Gij[ij]*aa[k]*aa[l]*fac;
			    }
			}
			}
			Xijkl[i][j][k][l] = 0.5e0*tmp*UtoGPa/V;
		    }
		}
	    }
	}

        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                for (k=0; k<3; k++) {
                    for (l=0; l<3; l++) {
			Cijkl[i][k][j][l] = Xijkl[i][j][k][l] + Xijkl[k][j][i][l] - Xijkl[i][k][j][l];
			if (i==j) Cijkl[i][j][k][l] -= ss[k][l];
			if (k==j) Cijkl[i][j][k][l] -= ss[i][l];
			if (k==i) Cijkl[i][j][k][l] += ss[j][l];
		    }
		}
	    }
	}

        for (i=0; i<3; i++) {
            for (j=i; j<3; j++) {
		int ij;
		if (i==1 && j==2) ij = 3;
		else if (i==0 && j==2) ij = 4;
		else if (i==0 && j==1) ij = 5;
		else ij = i;
                for (k=0; k<3; k++) {
                    for (l=k; l<3; l++) {
			int kl;
			if (k==1 && l==2) kl = 3;
			else if (k==0 && l==2) kl = 4;
			else if (k==0 && l==1) kl = 5;
			else kl = k;
			Cij[ij][kl] = (Cijkl[i][j][k][l]
				     + Cijkl[i][j][l][k]
				     + Cijkl[j][i][k][l]
				     + Cijkl[j][i][l][k]
				     + Cijkl[k][l][i][j]
				     + Cijkl[k][l][j][i]
				     + Cijkl[l][k][i][j]
				     + Cijkl[l][k][j][i])*0.125e0;
		    }
		}
	    }
	}

        for (i=0; i<6; i++) {
            for (j=0; j<6; j++) printf(" %10.2lf", Cij[i][j]);
	    printf("\n");
	}
}

double LatDYN::avec[3][3];
double LatDYN::inv_avec[3][3];
double LatDYN::scell_inv_avec[3][3];
double LatDYN::qvec[3][3];
double LatDYN::scell[3][3];
RDIS **LatDYN::disR;
SymRDIS **LatDYN::disNR;
CINDEX **LatDYN::cellA;
double *LatDYN::Fij;
int LatDYN::nTHz;
int LatDYN::Nfreq=0;
int LatDYN::Debug=0;

double *LatDYN::runa;
double *LatDYN::runb;
double *LatDYN::runc;
dcmplx *LatDYN::dmat;
dcmplx *LatDYN::dmatNAG;
double *LatDYN::dmatG;

int LatDYN::pvdosN;
int *LatDYN::xeN;
double *LatDYN::xe;
double *LatDYN::pvdos;

double LatDYN::BornCC;
double *LatDYN::frequencies;
gsl_vector *LatDYN::eval;
gsl_matrix_complex *LatDYN::evec;
gsl_matrix_complex *LatDYN::m;
gsl_eigen_herm_workspace *LatDYN::w;
gsl_eigen_hermv_workspace *LatDYN::wv;
gsl_matrix*LatDYN::evecR;
gsl_matrix*LatDYN::mR;
gsl_eigen_symm_workspace *LatDYN::wR;

int main(int argc, char* argv[])
{
FILE * fp_in = stdin;
char * line = NULL;
size_t len = 0;
int natomP, natomS, Ncell, Kcell;
double avec[3][3], inv_avec[3][3];
double bvec[3][3];
double gvec[3][3];
double scell[3][3], inv_scell[3][3];
double Rcell[3][3], inv_Rcell[3][3];
double a[3], b[3], c[3];
int i, j, k;
ATM **atomS, **atomP;
CINDEX **cellI, **cellA;
char *pdisfile = 0, *Bornfile = 0, *Gammafile=0, *Gfile = 0, *Rfile = 0;
int tranRS = 0;
const char *symmetry=0;
const char *expt=".dat";
int mall = 0, RedDM=0;
int nsym = 3, NdisQ=501;
int print = 0, silent = 0, tranI=1, plot=0, mode=1;
int nqx=-1, nqy=-1, nqz=-1, kvec=0, calcNA=1;
double funit=1.e0, eunit=1.e0;
double sigma=0.e0, fmaxT=-1.e0;
int sqs = 0;
int pvdos = 0;
FILE *born=0;
char **option;
int Nopt, nofetalerror=0;
int NEDOS=10001;
int PQ=0;
int eps = 0;
double DebCut = -1.0e0;
double Amode = 0.0;

	CPUTIM *cputim = new CPUTIM(); //starting CPU counting
	ostype = getenv("OSTYPE"); //return the environmental operating system
	shell = getenv("SHELL"); //rerurn the environmental Linux SHELL
	if (!ostype) ostype = "Unknown system"; //sometimes, bash does not return OSTYPE
	if (!shell) shell = "/bin/bash"; //if no SHELL environment, assuming bash
	
	setamass();

	Nopt = argc;
	option = new char*[argc];
	for (int i=0; i<argc; i++) option[i] = argv[i];
	FILE *optionfile = fopen("run", "w"); //save the options in a file "run" so that it remember them
						// then it will be easied for next run
	fprintf(optionfile, "#!/bin/csh -f\n");
	for (int i=0; i<argc; i++) fprintf(optionfile, "%s ", option[i]);
	fprintf(optionfile, "\n");
	fclose(optionfile);
        if (!strcmp(ostype,"linux")) foo = system("chmod u+x run");

	datetimehost(option, Nopt);

	/* handling command line options */
        for (int i=1; i<Nopt; i++) {
            if (!strcmp(option[i], "-Debug")) {
                LatDYN::Debug = 1; //turn on debug mode
            } else if (!strcmp(option[i], "-mode")) {
                mode = atoi(option[++i]); //noe used
            } else if (!strcmp(option[i], "-pvdos")) {
                pvdos = 1; //GPDOS calculation
            } else if (!strcmp(option[i], "-eps")) {
                eps = 1; //output figure in postscript file
            } else if (!strcmp(option[i], "-PQ")) {
                PQ = 1; //for used by thermal conductivity by Huazhi
            } else if (!strncmp(option[i], "-nof",4)) {
                nofetalerror = 1; //do not check fetal error
            } else if (!strncmp(option[i], "-sqs",4)) {
                sqs = 1; //SQS phonon calculation
            } else if (!strcmp(option[i], "-sigma")) {
                sigma = atof(option[++i]); //not USE, Parlinski parameter
	        RedDM = 0;
            } else if (!strcmp(option[i], "-DebCut")) {
                DebCut = atof(option[++i]); //Debye cutoff frequecies for long wave A*w**2 fitting
            } else if (!strcmp(option[i], "-Amode")) {
                Amode = atof(option[++i]); //Amplitude for unstable phonon mode
            } else if (!strcmp(option[i], "-thr")) {
                THR = atof(option[++i]); //general threshold
            } else if (!strcmp(option[i], "-thrE")) {
                THRE = atof(option[++i]);
            } else if (!strcmp(option[i], "-thr2")) {
                THR2 = atof(option[++i]); //threshold for restring symmtry
            } else if (!strcmp(option[i], "-fmax")) {
                fmaxT = atof(option[++i]); //threshold for restring symmtry
            } else if (!strcmp(option[i], "-nedos")) {
                NEDOS = atoi(option[++i]); //number of points for Phonon DOS calculation
            } else if (!strcmp(option[i], "-thr3")) {
                THR3 = atof(option[++i]); //not USED
            } else if (!strcmp(option[i], "-noNA")) {
                calcNA = 0; //do not consider D-D interaction
            } else if (!strcmp(option[i], "-bvec")) {
                kvec = 1; //using the primitive unit cell from the BEC file
            } else if (!strcmp(option[i], "-bvec1")) {
                kvec = 2;
            } else if (!strcmp(option[i], "-nrdis")) {
                NRDIS = atoi(option[++i]); //for speeding up calculations
            } else if (!strcmp(option[i], "-tranI")) {
                tranI = atoi(option[++i]); //for methods of make accoustic rule
            } else if (!strcmp(option[i], "-plot")) {
                plot = 1; //plot using THz unit
		funit =1.e0;
            } else if (!strcmp(option[i], "-plotTHz")) {
                plot = 1;
		funit =1.e0;
            } else if (!strcmp(option[i], "-plotcm")) {
                plot = 2; //plot using the cm-1 unit
		funit =1.e10/LIGHTSPEED;
		eunit =1.e10/LIGHTSPEED;
            } else if (!strcmp(option[i], "-plotmeV")) {
                plot = 3; //plot using the meV unit
		funit = 1.e0/meVtoTHz;
		eunit = 1.e0/meVtoTHz;
            } else if (!strcmp(option[i], "-emeVtoTHz")) {
		eunit = meVtoTHz; //conversion factor for experimental data
            } else if (!strcmp(option[i], "-print")) {
                print = atoi(option[++i]); //print option
            } else if (!strcmp(option[i], "-expt")) {
                expt = option[++i]; //name of the experimental data file
            } else if (!strcmp(option[i], "-silent")) {
                silent = 1; //do not print too much in the screen
            } else if (!strcmp(option[i], "-symm")) {
                symmetry = option[++i]; //define the symmetry of the crystall
            } else if (!strcmp(option[i], "-reddm")) {
                RedDM = 1-RedDM; //using the nonreduced DM or the reduced DM
            } else if (!strcmp(option[i], "-dmode")) {
                DMODE = 1-DMODE; //not USED
            } else if (!strcmp(option[i], "-ctldm")) {
                CTLDM = 1-CTLDM; //not USED
            } else if (!strcmp(option[i], "-tranRS")) {
                tranRS = 1;
            } else if (!strcmp(option[i], "-pdis")) {
                pdisfile = option[++i]; //name of the dispersion definition file
		if (access(pdisfile, R_OK)) {
		    printf("\n*******CANNOT access phonon dispersion path file \"%s\"\n\n", pdisfile);
		    exit(1);
		}
            } else if (!strcmp(option[i], "-Mass") || !strcmp(option[i], "-mass")) {
                char *Massfile = option[++i]; //name of the mass definition file
		if (access(Massfile, R_OK)) {
		    printf("\n*******CANNOT access mass redefinition file \"%s\"\n\n", Massfile);
		    exit(1);
		}
		addamass(Massfile);
            } else if (!strcmp(option[i], "-Born") || !strcmp(option[i], "-born")) {
                Bornfile = option[++i]; //name of the BEC definition file
		if (access(Bornfile, R_OK)) {
		    printf("\n*******CANNOT access Born effective tensor file \"%s\"\n\n", Bornfile);
		    exit(1);
		}
            } else if (!strcmp(option[i], "-Gfile") || !strcmp(option[i], "-gfile")) {
                Gfile = option[++i]; //name of the group representation definition file of phonon at Gamma point
		if (access(Gfile, R_OK)) {
		    printf("\n*******CANNOT access Gamma point irreducible representation file \"%s\"\n\n", Gfile);
		    exit(1);
		}
            } else if (!strcmp(option[i], "-Rfile") || !strcmp(option[i], "-rfile")) {
                Rfile = option[++i]; // //name of the file containing the rotation operations
		if (access(Rfile, R_OK)) {
		    printf("\n*******CANNOT access ratation symmetry restroe file \"%s\"\n\n", Rfile);
		    exit(1);
		}
            } else if (!strcmp(option[i], "-GC")) {
                Gammafile = option[++i]; //not USED
		if (access(Gammafile, R_OK)) {
		    printf("\n*******CANNOT access file \"%s\"\n\n", Gammafile);
		    exit(1);
		}
            } else if (!strcmp(option[i], "-mall")) {
                mall = 1; //average over the whole cell for translational symmetry
            } else if (!strcmp(option[i], "-nsym")) {
                nsym = atoi(option[++i]); //key for phonon DOS smoothness calculation
            } else if (!strcmp(option[i], "-ndisq")) {
                NdisQ = atoi(option[++i]); //number of points in each dispersion
            } else if (!strcmp(option[i], "-dirE")) {
                dirE[0] = atof(option[++i]); //direction of electric field
                dirE[1] = atof(option[++i]);
                dirE[2] = atof(option[++i]);
		Anormal(dirE);
		tranI = 2;
            } else if (!strcmp(option[i], "-nq")) {
                nqx = atoi(option[++i]); //define q mesh
                nqy = atoi(option[++i]);
                nqz = atoi(option[++i]);
            } else {
		printf("\n*******Fetal UNKNOWN OPTION \"%s\"\n\n", option[i]);
		exit(1);
	    }
        }

	FILE *fexact = fopen("exactQ.out", "w"); //file saving the eigenvalue and vectors for the exact Q

	if (Bornfile!=0 && calcNA==1) calcNA = 1;
	else calcNA = 0;

    /* for dispersion calculation, set up matix for converting into the conventional dirction */
    if (pdisfile) {
        if (!symmetry) {
			symmetry = ".scc";
            if (strlen(pdisfile) >= 4)
                symmetry = pdisfile + strlen(pdisfile)-4;
        } else symmetry = ".scc";
    }
	
	/* 
	internal conversion for direction of phonon dispersion 
	*/
	double *lattice;
        if (!symmetry) lattice = scc;
	else if (!strcmp(symmetry, ".fcc")) lattice = fcc;
        else if (!strcmp(symmetry, ".bcc")) lattice = bcc;
        else if (!strcmp(symmetry, ".hcp")) lattice = hcp;
        else if (!strcmp(symmetry, "tet2")) lattice = tet2;
        else lattice = scc;

	for (i=0; i<3; i++) a[i] = dirE[i];
        mproduct (a, lattice, dirE);

	//skip comment line leading by '#'
        for (i=0; ;i++) {
	    long fpos = ftell(fp_in);
	    yetline(&line, &len, fp_in);
	    if (line[0] != '#') {
		fseek(fp_in, fpos, SEEK_SET);
                break;
	    }
	}

	/* read in the primitive unit cell lattice vectors */
	for (i=0; i<3; i++)
	    for (j=0; j<3; j++) foo = fscanf(fp_in, "%lf", &avec[i][j]);

        printf("\nPrimitive lattice vector from Supercell\n\n");
        for (i=0; i<3; i++) println("%18.12f  ",avec[i]);

	/* read in the supercell lattice vectors */
        printf("\nSupercell lattice vector\n\n");
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) foo = fscanf(fp_in, "%lf", &scell[i][j]);
            println("%18.12f  ",scell[i]);
	}

	/* read in the number of atoms in the supercell */
        foo = fscanf(fp_in, "%d %d", &natomS, &Kcell);

	/* read in BEC tensor */
	if (Bornfile!=0) {
	    born = fopen(Bornfile, "r");
		/* read in the primitive unit cell lattice vectors */
	    for (i=0; i<3; i++) 
		for (j=0; j<3; j++) foo = fscanf(born, "%lf", &bvec[i][j]);

            printf("\nPrimitive lattice vector from Bornfile\n\n");
            for (i=0; i<3; i++) println("%18.12f  ",bvec[i]);


	    /* for -bvec option */
	    if (kvec==1) {
		double ratio = volume(avec)/volume(bvec);
		if (ratio>=1.e0) {
		    ratio /= (int)(ratio+0.5);
		    ratio = pow(ratio, 1.e0/3.e0);
		} else ratio=1.e0;
//printf ("ratio=%lf abce=%lf, %lf\n", ratio, volume(avec), volume(bvec));
		for (i=0; i<3; i++) 
		    for(j=0;j<3;j++) {
			bvec[i][j] *= ratio;
			avec[i][j] = bvec[i][j];
		    }
	    }
	    else if (kvec==2)
		for (i=0; i<3; i++) 
		    for(j=0;j<3;j++)
			avec[i][j] = bvec[i][j];

	    inv_m(avec, inv_avec);

            yetline(&line, &len, born);

	    Kcell = (int)(0.5e0 + volume(scell)/volume(bvec));
	    natomP = natomS/Kcell;

	    atomP = new ATM *[natomP];

	    /* read in the atomic positions in the BEC file */
            printf("\n Atom position(s) in the Born effective charge file:\n\n");
	    for (int i=0; i<natomP; i++) {
                yetline(&line, &len, born);
		printf("%s", line);
		SplitItem *pos = new SplitItem(line);
		for (int k=0; k<3; k++) a[k] = normal(fmod(atof(pos->str(k)),1.e0));
		mproduct(a, bvec, b);
		mproduct(b, inv_avec, a);
		for (int k=0; k<3; k++) a[k] = normal(fmod(a[k],1.e0));
		atomP[i] = new ATM(a[0], a[1], a[2], pos->str(3));
	    }
	} else {
	    natomP = 0;
	    atomP = (ATM **) malloc ( (size_t) (sizeof(ATM *)) );
	}

	inv_m(avec, inv_avec);
	//July 1, 2013 major change
	//inv_m(scell, inv_scell);
	mproduct(scell, inv_avec, Rcell);

        printf("\nSupercell is a transformation of the Primitive cell by\n\n");
        for (i=0; i<3; i++) println("%18.12f  ",Rcell[i]);

	int error = 0;
        printf("\nThe tranforming Matrix is truncated into\n\n");
        for (int i=0; i<3; i++) { 
                for (j=0; j<3; j++) {
		    double tmp = floor(0.5+Rcell[i][j]);
		    //if (abs(Rcell[i][j]-tmp) > THR) error = 1;
		    if (fabs(Rcell[i][j]-tmp) > 0.4e0) error = 1;
		    Rcell[i][j] = tmp;
		}
                println("%18.12f  ", Rcell[i]);
        }

	if (error == 1) {
	    printf ("\n********Fetal ERROR, the truncating is a little larger ! May be due to the Born file\n\n");
	    exit(1);
	}

	inv_m(Rcell, inv_Rcell);
	if (kvec==1 || kvec==2) mproduct (Rcell, avec, scell);
	else mproduct (inv_Rcell, scell, avec);

	inv_m(avec, inv_avec);
	//July 1, 2013 major change
	inv_m(scell, inv_scell);

        yetline(&line, &len, fp_in);
        yetline(&line, &len, fp_in);
	int Direct = 1;
	if (tolower(line[0])=='c') Direct = 0;

	atomS = (ATM **) malloc ( (size_t) (natomS*sizeof(ATM *)) ); //locate memory for supercell atomic positions
        printf("\n%4d Atom(s) in the Supercell\n\n", natomS);

	    /* read in the atomic positions in the supercell file */
        for (int j=0; j<natomS; j++) {
        	yetline(&line, &len, fp_in);
		SplitItem *pos = new SplitItem(line);
		for (int k=0; k<3; k++) {
			b[k] = atof(pos->str(k));
			a[k] = b[k];
		}
		if (Direct==0) mproduct(b,inv_scell, a);
		for (int k=0; k<3; k++) a[k] = normal(fmod(a[k],1.e0));

		atomS[j] = new ATM(a[0], a[1], a[2], pos->str(3));
		if (LatDYN::Debug) printf(" %.10lf %.10lf %.10lf %2s %4d\n", 
		    atomS[j]->x, atomS[j]->y, atomS[j]->z, atomS[j]->sym, j);
	}

	/* indexing the primitive unit cell in the supercell */
	Ncell = 0;
	cellA = (CINDEX **) malloc ( (size_t) natomS*(sizeof(CINDEX *)) );
	cellI = (CINDEX **) malloc ( (size_t) natomS*(sizeof(CINDEX *)) );

        if (LatDYN::Debug) 
	    printf("\nPrimitive atom index in the Supercell\n\n");

        for (int n=0; n<natomS;n++) { 
		a[0] = atomS[n]->x;
		a[1] = atomS[n]->y;
		a[2] = atomS[n]->z;
		mproduct(a, scell, c);
		mproduct(c, inv_avec, b);
		int *Ix = Inormal(b);
                CINDEX *tx = new CINDEX(Ix[0], Ix[1], Ix[2]);
                int Icell;
                for (Icell=0; Icell<Ncell;Icell++) {
                    if ((*tx) == (*cellI[Icell])) break;
                }
                if (Icell==Ncell) {
                    if (LatDYN::Debug) printf (" %3d %3d %3d %3d %lf %lf %lf\n", 
			Ncell, Ix[0], Ix[1], Ix[2], b[0], b[1], b[2]);
                    cellI[Ncell++] = tx;
		}

		for (int jj=0; jj<3; jj++) c[jj] = normal(b[jj]);
		/* double check the data slefconsistant */
		double check[3];
		for (int jj=0; jj<3; jj++) check[jj] = b[jj] - (double)Ix[jj];
		for (int jj=0; jj<3; jj++) {
		    if (fabs(check[jj] - c[jj])<THR) continue;
		    printf("********FETAL ERROR! non-selfconsistant cordinate!\n\n");
		    printf(" primitive coordinate from normal function is :\n\n");
		    println(stdout, " %10.6lf", c);
		    printf(" primitive coordinate from Inormal function is :\n\n");
		    println(stdout, " %10.6lf", check);
		}

		//ATM *tmp = new ATM(c[0], c[1], c[2], atomS[n]->sym);
		ATM *tmp = new ATM(check[0], check[1], check[2], atomS[n]->sym);

                int ii;
		for (ii=0; ii<natomP; ii++) {
			if (sqs!=0) {
			    if (tmp->equalSQS(*atomP[ii])) break;
			} else {
			    if ((*tmp) == (*atomP[ii])) break;
			}
		}
		if (ii!=natomP) {
		    cellA[n] = new CINDEX(Icell,ii,n);
		    double bb[3], cc[3], dd[3];
		    double *aa = atomP[ii]->getxyz();
		    mproduct(aa, avec, bb);
		    for (int jj=0; jj<3; jj++) dd[jj] = Ix[jj];
		    mproduct(dd, avec, cc);
		    for (int jj=0; jj<3; jj++) bb[jj] += cc[jj];
		    mproduct(bb, inv_scell, cc);
		    if (RedDM==1 || kvec==2) atomS[n]->setxyz(cc);
		} else if (Bornfile!=0) {
printf("why sqs compare = %d, ii=%d, natomP=%d\n", sqs, ii, natomP);
	            printf("\n******** FETAL ERROR: cannot find atom %s %d %lf %lf %lf from %lf %lf %lf\n"
		        , atomS[n]->sym, n, c[0], c[1], c[2], atomS[n]->x, atomS[n]->y, atomS[n]->z);
	            printf("\n******** try -thr2 different from %lf to solve it\n", THR2);
        	    for (int j=0; j<natomP; j++) {
                	printf(" %.10lf %.10lf %.10lf %2s %4d\n",
                        atomP[j]->x, atomP[j]->y, atomP[j]->z, atomP[j]->sym, j);
        	    }	
	    	    exit(1);
		} else {
		    cellA[n] = new CINDEX(Icell,natomP,n);
		    atomP = (ATM **) realloc (atomP, (size_t) ((++natomP)*sizeof(ATM *)) );
		    atomP[natomP-1] = tmp;
		}
	}

	if (LatDYN::Debug) {
	    printf("\nTranformed Supercell\n\n");
            for (int n=0; n<natomS;n++)
	        printf(" %.10lf %.10lf %.10lf %2s %4d\n", 
		    atomS[n]->x, atomS[n]->y, atomS[n]->z, atomS[n]->sym, n);
	}

        printf("\nUsed Primitive lattice vector\n\n");
        fprintf(fexact, "\nUsed Primitive lattice vector\n\n");

	for (i=0; i<3; i++) println("%18.12f  ", avec[i]);
	for (i=0; i<3; i++) println(fexact, "%18.12f  ", avec[i]);

        printf("\n%4d Atom(s) in the Primitive cell\n\n", natomP);
        fprintf(fexact, "\n%4d Atom(s) in the Primitive cell\n\n", natomP);

	/* make the POSCAR file to be used by pos2s for symmetry analysis */
	FILE *symmfile = fopen("Symmetry.pos", "w");
        fprintf(symmfile, "Used Primitive lattice vector\n 1.00\n");
	for (i=0; i<3; i++) println(symmfile, "%18.12f  ", avec[i]);

        int *n = new int[natomP];
        int N = 1;
        n[0] = 1;

        fprintf(symmfile, " %4s", atomP[0]->sym);
        for (int i=1; i<natomP; i++) {
            if (strcmp(atomP[i-1]->sym, atomP[i]->sym)) {
                fprintf(symmfile, " %4s", atomP[i]->sym);
                n[N] = 1;
                N++;
            } else {
                n[N-1]++;
            }
        }
        fprintf(symmfile, "\n");

        for (int i=0; i<N; i++) fprintf(symmfile, " %4d", n[i]);
        fprintf(symmfile, "\nD\n");

        for (int j=0; j<natomP; j++) {
                printf(" %.10lf %.10lf %.10lf %2s %4d\n",
                        atomP[j]->x, atomP[j]->y, atomP[j]->z, atomP[j]->sym, j);
                fprintf(symmfile, " %.10lf %.10lf %.10lf %2s %4d\n",
                        atomP[j]->x, atomP[j]->y, atomP[j]->z, atomP[j]->sym, j);
                fprintf(fexact, " %.10lf %.10lf %.10lf %2s %4d\n",
                        atomP[j]->x, atomP[j]->y, atomP[j]->z, atomP[j]->sym, j);
        }
	fclose(symmfile);

        printf("\n Mapping of Supercell of Primitive cell\n\n");

	/* MM labels the reference central PUC at which the DM (Nx is its number) is calculated or averaged on */
        int *MM = new int[natomS];
        double R = 0.e0;
        for (int i=0;i<natomP; i++) {
	  printf ("Atom %2d %2.2s", i, atomP[i]->sym);
          double RX = 1.e72;
          for (int M=0; M<natomS;M++) {
            if (cellA[M]->get_atomP()!=i) continue;
	    printf (" %3d", M);
            double tmp=0.e0;
            for (int N=0; N<natomS;N++) {
                a[0] = atomS[N]->x - atomS[M]->x;
                a[1] = atomS[N]->y - atomS[M]->y;
                a[2] = atomS[N]->z - atomS[M]->z;
                mproduct(a, scell, c);
                if (normal(c)>tmp) tmp = normal(c);
            }
            if (tmp < RX) {
                RX = tmp;
                MM[i] = M;
            }
          }
	  printf ("\n");
	  R = max(R,RX);
        }

        int NX = natomP;
	if (mall==1) {
            NX = natomS;
	    for (int M=0;M<natomS;M++) MM[M] = M;
	}

        printf("\nMAXimum interaction DIS = %lf, Dynamical Matrix will be averaged over %d cell(s) built by %d centered atom(s)\n\n", R, NX/natomP, NX);

        for (int i=0; i<NX; i++) {
	    int M = MM[i];
            printf(" %2s PI= %3d SI= %3d at %10.6lf %10.6lf %10.6lf Index= %3d %3d %3d\n",
		atomS[M]->sym, cellA[M]->get_atomP(), M, 
		atomS[M]->x, atomS[M]->y, atomS[M]->z, 
		cellI[cellA[M]->get_cell()]->get_cx(),
		cellI[cellA[M]->get_cell()]->get_cy(),
		cellI[cellA[M]->get_cell()]->get_cz());
	}

        if (LatDYN::Debug) printf("\nVerifying anslational operations and they are:\n\n");

	int M0 = cellA[0]->get_atomP();
	int Ntrans = 0;
        for (int i=0; i<natomS; i++) {
	    int M = cellA[i]->get_atomP();
	    if (M!=M0) continue;
            if (LatDYN::Debug) printf("  Tran %3d  Oper= %3d %3d %3d\n",
		Ntrans,
		cellI[cellA[i]->get_cell()]->get_cx(),
		cellI[cellA[i]->get_cell()]->get_cy(),
		cellI[cellA[i]->get_cell()]->get_cz());
	    Ntrans++;
	}
        if (LatDYN::Debug) printf("\nFound %d translational operations\n\n", Ntrans);

	if (!silent) printf("\nMapping supercell :");
	if (!silent) cputim->elptime(stdout);

        int dimN = natomS*3; //number of phonons at the exact Q point
        int Rfreq = dimN*dimN;
        double *Fij_Orig = new double[Rfreq];
        double *Fij = new double[Rfreq];
        double *wFij = new double[Rfreq];
        for (i=0;i<Rfreq;i++) foo = fscanf(fp_in, "%lf", &Fij_Orig[i]); //read in VASP.5 force constants

	RMATRAX *rmatrix=0;

	if (Rfile) {
	    /* symmetrimize the force constants (fij) using the input rotation operations */
	    rmatrix = new RMATRAX(Rfile);
	    rmatrix->ResRotSym(Fij_Orig, natomS, atomS, scell, kvec);
            if (tranRS) rmatrix->ResTransSym(Fij_Orig, natomS, atomS, avec, scell, cellI, cellA, Ntrans);
	}
	    
        for (i=0;i<Rfreq;i++) Fij[i] = Fij_Orig[i];

        double *stress = new double[6];
	
	/* not useful. read in the remained stresses from the superfij file */
        for (i=0;i<6;i++) {
	    foo = fscanf(fp_in, "%lf", &stress[i]);
	    if (feof(fp_in)) {
		for (j=0; j<6; j++) stress[j] = 0.e0;
		break;
	    }
	}

	/* Symmetrization fij matrix and apply ASR */

        for (int i=0; i<dimN; i++) {
            for (int j=i; j<dimN; j++) {
		int ij = i*dimN+j;
		int ji = j*dimN+i;
                if (tranI==0) {
                    Fij[ij] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                    Fij[ji] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                } else if (tranI==1) {
                    Fij[ij] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                    Fij[ji] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                } else if (tranI==2) {
                    Fij[ij] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                    Fij[ji] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
		    wFij[ij] = 1.e0;
		    wFij[ji] = 1.e0;
                } else if (tranI==3) {
                    Fij[ij] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                    Fij[ji] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
		    wFij[ij] = Fij[ij];
		    wFij[ji] = wFij[ij];
                } else if (tranI==4) {
                    Fij[ij] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                    Fij[ji] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
		    wFij[ij] = Fij_Orig[ij];
		    wFij[ji] = Fij_Orig[ji];
                } else if (tranI==7) {
                    Fij[ij] = Fij_Orig[ij];
                    Fij[ji] = Fij_Orig[ij];
                } else if (tranI==8) {
                    Fij[ij] = Fij_Orig[ji];
                    Fij[ji] = Fij_Orig[ji];
                } else {
                    Fij[ij] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                    Fij[ji] = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
                }
	    }
	}

	if (tranI==2 || tranI==3 || tranI==4) {
	  double dout=1.e72, dnew=0.e0;
    	  for (int iter=0; iter<99 && fabs(dnew-dout)>1.e-12; iter++) {
	    dout = dnew;
            dnew = makeTI(Fij, natomS, wFij);
            if (LatDYN::Debug)
	      printf ("Fij dout = %g\n", dout);
          }   
	} 

	/* apllying different accoustic summation rule, no much effects */
	/* applying ASR */
        double *linS = new double[9];
        double *VlinS = new double[9];
        for (int i=0; i<natomS; i++) {
            //int I = cellA[i]->get_atomP();
            for (int k=0;k<9;k++) linS[k] = 0.e0;
            for (int k=0;k<9;k++) VlinS[k] = 0.e0;
            for (int x=0; x<3; x++) {
                for (int j=0; j<natomS; j++) {
                    for (int y=0; y<3; y++) {
                        int ii = i*3 + x;
                        int jj = j*3 + y;
                        int iijj = ii*dimN + jj;
                        int jjii = jj*dimN + ii;
                        linS[x*3+y] += Fij[iijj];
                        VlinS[x*3+y] += Fij[jjii];
                    }
                }
	    }

            for (int x=0; x<3; x++) {
                for (int y=0; y<3; y++) {
                    int ii = i*3 + x;
                    int jj = i*3 + y;
                    int iijj = ii*dimN + jj;
            	    //if (LatDYN::Debug)
		    //	printf(" linS=%g, VlinS=%g for %d,%d\n", linS[x*3+y], VlinS[x*3+y], ii,jj);
		    if (tranI) Fij[iijj] -= linS[x*3+y];
		}
	    }

	    if (tranI) {
              for (int x=0; x<3; x++) {
                //for (int y=0; y<x; y++) {
                for (int y=x+1; y<3; y++) {
                    int ii = i*3 + x;
                    int jj = i*3 + y;
                    int iijj = ii*dimN + jj;
                    int jjii = jj*dimN + ii;
		    Fij[iijj] = Fij[jjii];
		}
	      }
	    }
	}
	delete linS;
	delete VlinS;

        for (int i=0; i<dimN; i++) {
            for (int j=i; j<dimN; j++) {
		int ij = i*dimN+j;
		int ji = j*dimN+i;
                double tmp = 0.5e0*(Fij_Orig[ij]+Fij_Orig[ji]);
		if (tranI==2) tmp = 0.5e0*(Fij[ij]+Fij[ji]);
		Fij_Orig[ij] = tmp;
		Fij_Orig[ji] = tmp;
	    }
	}

        double *sij = new double[9];

	int nTHz = natomP*3;
        int Gfreq = nTHz*nTHz;
        double *Sij = new double[Gfreq];
        double *Mij = new double[Gfreq];
        for (int i=0; i<Gfreq; i++) Sij[i] = 0.e0;
        for (int i=0; i<Gfreq; i++) Mij[i] = 0.e0;

	/* Symmetrize the force constant matrix */
        for (int i=0; i<natomS; i++) {
            int I = cellA[i]->get_atomP();
            for (int x=0; x<3; x++) {
                for (int j=0; j<natomS; j++) {
                    int J = cellA[j]->get_atomP();
                    for (int y=0; y<3; y++) {
                        int ii = i*3 + x;
                        int jj = j*3 + y;
                        int iijj = ii*dimN + jj;
                        //int jjii = jj*dimN + ii;
                        int II = I*3 + x;
                        int JJ = J*3 + y;
                        int IIJJ = II*nTHz + JJ;

                        Sij[IIJJ] += Fij[iijj];

                        int N;
                        for (N=0;N<NX;N++) if(MM[N]==i) break;
                        if (N==NX) continue;
                        Mij[IIJJ] += Fij[iijj];
                    }
                }
            }
        }

	if (LatDYN::Debug) {
	    printf("\nDerived elastic constants from Fij (GPa)\n\n"); //not accurate, will be deleted
	    Cij(Fij, natomS, scell, atomS, stress);
	    if (!silent) printf("\nDeriving Elastic constants :");
	    if (!silent) cputim->elptime(stdout);
	}

        for (int i=0; i<Gfreq; i++) Sij[i] /= (double)(natomS/natomP);
        for (int i=0; i<Gfreq; i++) Mij[i] /= (double)(NX/natomP);

	if (LatDYN::Debug) {
	    //not accurate, will be deleted
	    printf("\nDerived elastic constants from Sij(GPa)\n\n");
	    Cij(Sij, natomP, avec, atomP, stress);
	}

	/* checking the mappings between the primitive unit cell and supercell */
        for (int i=0; i<natomS; i++) {
            int I = cellA[i]->get_atomP();
            //double imass = getmass(atomP[I]->sym);
            double imass = getmass(atomS[i]->sym);
	    if (strcmp(atomP[I]->sym, atomS[i]->sym)) {
	        printf ("******Fetal atom mappings not correct from Supercell %d %s to prim %d %s\n",i, atomS[i]->sym, I, atomP[I]->sym);
		if (!nofetalerror) exit(1);
	    }
            for (int x=0; x<3; x++) {
                int ii = i*3 + x;
                for (int j=0; j<natomS; j++) {
                    int J = cellA[j]->get_atomP();
                    //double jmass = getmass(atomP[J]->sym);
                    double jmass = getmass(atomS[j]->sym);
                    double fac = 1.e0/sqrt(imass*jmass);
                    for (int y=0; y<3; y++) {
                        int jj = j*3 + y;
                        int iijj = ii*dimN + jj;
                        //int jjii = jj*dimN + ii;
			if (tranI) {
                            int II = I*3 + x;
                            int JJ = J*3 + y;
                            int IIJJ = II*nTHz + JJ;
                            Fij[iijj] -= Mij[IIJJ]/(double)Kcell;
          		}
                        Fij[iijj] *= fac;
                        Fij_Orig[iijj] *= fac;
                    }
                }
            }
	}

	/* for the purpose of high accurate calculation at Gamma point 
	  under developping
	*/
        double *Gij = new double[Gfreq];
	if (Gammafile!=0) {
	    FILE *gammafile = fopen(Gammafile, "r");
	    for (i=0; i<3; i++)
                for (j=0; j<3; j++) foo = fscanf(gammafile, "%lf", &gvec[i][j]);
	    for (i=0; i<3; i++)
                for (j=0; j<3; j++) foo = fscanf(gammafile, "%lf", &gvec[i][j]);
            for (i=0; i<=natomP+2;i++) yetline(&line, &len, gammafile);

	    for (i=0;i<Gfreq;i++) foo = fscanf(gammafile, "%lf", &Gij[i]);
            double *stressGC = new double[6];
            for (i=0;i<6;i++) {
	        foo = fscanf(fp_in, "%lf", &stressGC[i]);
	        if (feof(fp_in)) {
		    for (j=0; j<6; j++) stressGC[j] = 0.e0;
		    break;
	        }  
	    }
	    fclose(gammafile);

	    if (tranI) {
        	for (int i=0; i<nTHz; i++) {
     	           for (int j=i; j<nTHz; j++) {
                	int ij = i*nTHz+j;
                	int ji = j*nTHz+i;
                	double tmp = 0.5e0*(Gij[ij]+Gij[ji]);
                	Gij[ij] = tmp;
                	Gij[ji] = tmp;
            	   }	
        	}
	    }

	    if (LatDYN::Debug) {
	        //not accurate, will be deleted
	        printf("\nDerived elastic constants from Gij(GPa)\n\n");
	        Cij(Gij, natomP, avec, atomP, stressGC);
	    }
	}

	/* devide the fij by atomic mass, i.e. make Hessian */
        for (int i=0; i<natomP; i++) {
            double imass = getmass(atomP[i]->sym);
            for (int ii=0;ii<9;ii++) sij[ii] = 0.e0;

            for (int j=0; j<natomP; j++) {
                double jmass = getmass(atomP[j]->sym);
                double fac = 1.e0/sqrt(imass*jmass);
                for (int x=0; x<3; x++) {
                  for (int y=0; y<3; y++) {
                        int ii = i*3 + x;
                        int jj = j*3 + y;
                        int iijj = ii*nTHz + jj;
                        sij[x*3+y] += Sij[iijj];
                        Sij[iijj] *= fac;
                        if (Gammafile!=0) Gij[iijj] *= fac;
		  }
                }
            }
            //printf("\nHorizonal Summation of Hessian at Gamma point for Atom %d\n\n", i);
            //for (int ii=0; ii<3; ii++) println("%18.12f  ",sij+ii*3);
        }

	if (!silent) printf("\nRestore the translational invariance :");
	if (!silent) cputim->elptime(stdout);

	/* read in the BEC and dielectric constat tensor */
	double *dielec = new double[9];
	double *BCharg = new double[9*natomP];
	if (calcNA!=0) {
	    printf ("\nMacroscopic high freq dielectric field tensor\n\n");
	    double *pd = dielec;
	    for (int j=0; j<3; j++) {
		yetline(&line, &len, born);
		sscanf (line, "%lf %lf %lf", pd,  pd+1,  pd+2);
		println (" %11.5lf ", pd);
		pd += 3;
	    }

            for (int ii=0;ii<9;ii++) sij[ii] = 0.e0;
	    if (LatDYN::Debug) printf ("\nBorn effective charge tensor\n");
	    for (int i=0; i<natomP; i++) {
		if (LatDYN::Debug) printf ("\nion %3d at %11.6lf %11.6lf %11.6lf %s\n\n", 
		    i, atomP[i]->x, atomP[i]->y, atomP[i]->z, atomP[i]->sym);
		yetline(&line, &len, born);
	        pd = BCharg+i*9;
		int tmp;
		for (int j=0; j<3; j++) {
		    yetline(&line, &len, born);
		    sscanf (line, "%d %lf %lf %lf", &tmp, pd,  pd+1,  pd+2);
		    if (LatDYN::Debug) println (" %11.5lf ", pd);
		    pd += 3;
		}
		//symmetrimization of BEC
/*
		for (int j=0; j<3; j++)  {
		    for (int k=j; k<3; k++) {
			double tmp = (BCharg[i*9+j*3+k] + BCharg[i*9+k*3+j])*0.5e0;
			BCharg[i*9+j*3+k] = tmp;
			BCharg[i*9+k*3+j] = tmp;
		    }
		}
*/
		for (int j=0; j<9; j++)  sij[j] += BCharg[i*9+j];
	    }

	    printf ("\nElectric neutral condition\n\n");
	    for (int j=0; j<3; j++)  println (" %11.5lf ", sij+j*3);
	    for (int ii=0; ii<natomP; ii++) 
	        for (int j=0; j<9; j++)
		    BCharg[ii*9+j] -= sij[j]/(double)natomP;

	    fclose(born);
	    if (rmatrix) rmatrix->ResRotSymBorn(dielec, BCharg, natomP, atomP, avec);
	    inv_m(dielec, inv_dielec);
	}

	delete sij;

	/* calculate the phonon frequencies at the real space, i.e., the exact Q */
	double *energy_in_real_space = new double[dimN];
	if (Amode <= 0.0) LatDYN::SrealFij(Fij_Orig, natomS, energy_in_real_space);
	else LatDYN::SrealFijV(Fij_Orig, natomS, energy_in_real_space, 
	    scell, inv_scell, atomS, Amode,"iPOSCAR.S");

	if (!silent) printf("\nSolving frequencies in real space :");
	if (!silent) cputim->elptime(stdout);

        double *energyG = new double[nTHz]; //phonons at Gamma point
        if (Amode <= 0.0) LatDYN::SrealFij(Sij, natomP, energyG);
	else LatDYN::SrealFijV(Sij, natomP, energyG,
	    avec, inv_avec, atomP, Amode,"iPOSCAR.G");


        double FACTOR = (double)natomP/(double)NX;

	if (nqx < 0) {
	    nqx = (int) pow(1.e6/(double)natomP, 1.e0/3.e0)+2;
	    nqy = nqx;
	    nqz = nqx;
	}
        printf ("\n Used q-mesh %dx%dx%d\n", nqx, nqy, nqz);

	int kx = nqx/2;
	int ky = nqy/2;
	int kz = nqz/2;
	double fx = .5e0/(double)kx;
	double fy = .5e0/(double)ky;
	double fz = .5e0/(double)kz;

	//set up workspace for GSL
	LatDYN::SetWorkSpace(avec,scell,atomS, atomP, 
	    cellI, cellA, Fij, natomS, natomP, MM, NX, mode, 
	    -0.1*energy_in_real_space[0], 1.01*energy_in_real_space[0]); //set up workspace for GSL

	if (!silent) printf("\nSetting workspace & pre-optimizing :");
	if (!silent) cputim->elptime(stdout);

        if (tranI) {
	    if (Gammafile==0) for (int i=0; i<nTHz*nTHz;i++)
	        LatDYN::dmatG[i] = Sij[i]/(double)Kcell;
	    else for (int i=0; i<nTHz*nTHz;i++)
	        LatDYN::dmatG[i] = Gij[i]/(double)Kcell;
	}

	int Nrep = 0;
	PhononG **gp=0;

	//read in group represetnation infor for Gamma point phonons
	if (Gfile) {
	    FILE *gfile = fopen(Gfile, "r"); //readin group represetnation infor for Gamma point phonons
	    gp = (PhononG **) malloc( (size_t) (sizeof(PhononG*)) );
	    while (1) {
	    	gp = (PhononG **) realloc(gp, (size_t) ((Nrep+1)*sizeof(PhononG *)) );
	    	gp[Nrep] = new PhononG(gfile, nTHz);
		if (gp[Nrep]->end_of_file) break;

		if (gp[Nrep]->IR) gp[Nrep]->SetDirP();
		gp[Nrep]->SrealFij(Sij, calcNA);

		Nrep++;
	    }
	    fclose(gfile);
	}

	if (!silent) printf("\nHandling symmetry :");
	if (!silent) cputim->elptime(stdout);

	/* prepare the menory for IR/Raman as well as LO/TO analysis */
	int exactP = 0;
	double *energyGB = new double[nTHz];
	double *ILOTO = new double[nTHz];
	double *IRed = new double[nTHz];
	RAMAN *IRaman = new RAMAN[nTHz];
	double *RRvv = new double [nTHz*nTHz];
	double *EEvv = new double[nTHz];
	double *EEss = new double[nTHz];
	char **irrep = new char*[nTHz];
	for (int i=0; i<nTHz; i++) irrep[i] = 0;

	LatDYN::setNFrq(natomS*3);
	/* loop over the q mesh, finding the exact Q within the Q mesh 
	   may not find all exact Q due to the Q mesh is not completly compatible with supercell
	*/
	printf ("\n");
	for (i=-kx; i<kx; i++) {
	    for (j=-ky; j<ky; j++) {
	    for (k=-kz; k<kz; k++) {
		a[0] = (double)i*fx;
		a[1] = (double)j*fy;
		a[2] = (double)k*fz;
		if (exactK(a,Rcell)) {
            	    fprintf (fexact, "\n\n exact Q at (%s", irri(i, 2*kx)); 
            	    fprintf (fexact, ", %s",  irri(j, 2*ky));
            	    fprintf (fexact, ", %s)",  irri(k, 2*kz));
                    LatDYN *DM = new LatDYN(a[0], a[1], a[2]);
                    if (calcNA!=0)
                        DM->BornGamma(atomP, natomP, dielec, BCharg, mode, Kcell, dirE);
                    if (CTLDM) DM->ReduceDM(natomS, MM, NX, FACTOR, RedDM, tranI, calcNA);
                    else DM->NonReduceDM(natomS, MM, NX, FACTOR, tranI, calcNA, sigma, sqs);
                    double foundunstable = DM->solveDMev(fexact, atomP, natomP, calcNA, BCharg, ILOTO, IRed, IRaman, RRvv);
		    if (foundunstable<0.0) 
			printf ("******** Found unstable mode with frequency %.3lf THz at exact Q (%s, %s, %s)\n", foundunstable, irri(i, 2*kx), irri(j, 2*ky), irri(k, 2*kz));

		    exactP++;
		    if (i==0 && j==0 && k==0) {
			for (int x=0; x <nTHz; x++) {
			    energyGB[x] = LatDYN::frequencies[LatDYN::Nfreq-nTHz+x];
			    EEvv[x] = energyGB[x];
			    EEss[x] = 0.e0;
			}
/*
            		gsl_vector *eval = gsl_vector_alloc (nTHz);
            		for (int x=0; x<nTHz; x++)
            		    gsl_vector_set (eval, x, energyGB[x]);
            		gsl_sort_vector(eval);

            		for (int x=0; x<nTHz; x++)
                		energyGB[x] = gsl_vector_get (eval, x);
            		gsl_vector_free(eval);
*/
			if (Gfile) {
			    if (calcNA!=0) {
					double tZmu[3], dion[9];
					for (int x=0; x<3; x++) tZmu[x] = 0.e0;
					for (int x=0; x<9; x++) dion[x] = dielec[x];
			        if (!silent) printf("\n\nSolving frequencies considering LO-TO splitting:\n");
			        for (int r=0; r<Nrep; r++) {
                        	    if ( gp[r]->dirP) DM->BornGamma(atomP, natomP, dielec, BCharg, mode, Kcell, gp[r]->dirP);
				    else DM->BornGamma(atomP, natomP, dielec, BCharg, mode, Kcell, dirE);
				    //gp[r]->SrealFij(Sij, Kcell, LatDYN::dmatNAG);
				    gp[r]->SrealDij(Sij, Kcell, LatDYN::dmatNAG, atomP, natomP, BCharg, tZmu);
				    gp[r]->Analysis(RRvv);
				    gp[r]->print_symmA(EEvv, EEss, irrep);
					for (int x=0; x<3; x++) {
						for (int y=0; y<3; y++) {
							dion[x*3+y] += tZmu[x]*tZmu[y]/PI2/PI2/volume(bvec)*1.e6;
						}
					}
				}

			    printf ("\nMacroscopic high freq dielectric field tensor with ionic contribution\n\n");
			    double *pd = dion;
			    for (int x=0; x<3; x++) {
					println (" %11.5lf ", pd);
					pd += 3;
			    }
			}
			}
			}
			}
            }
            }
        }
	fclose(fexact);
        printf("\n %d Exact Q points found\n\n", exactP);

	if (!silent) printf("\nFinding exact Q points :");
	if (!silent) cputim->elptime(stdout);

    if (LatDYN::Debug) {
	/* for debug use only */
        printf("\nCompare exact Q and Real Space eigenvalues:\n\n");
	LatDYN::printQ(0);
	if (LatDYN::Nfreq != dimN) {
	    printf("\n********Warning: The number of exact Q = %d notequal to %d,\n*******use -nq the change it\n\n", LatDYN::Nfreq/nTHz, dimN/nTHz);
	} else {  
	    double tmp = 0.e0;
            printf (" No. Real Space    Exact Q       Diff         meV       cm-1\n");
	    for (int i=0; i<dimN; i++) {
	        tmp += pow(energy_in_real_space[i]-LatDYN::frequencies[dimN-i-1],2);
                printf ("%4d %10.6lf %10.6lf %10.6lf, %10.5lf %10.4lf\n",
		    i,
                    energy_in_real_space[i], 
		    LatDYN::frequencies[dimN-i-1],
		    LatDYN::frequencies[dimN-i-1]-energy_in_real_space[i], 
		    LatDYN::frequencies[dimN-i-1]/meVtoTHz, 
		    LatDYN::frequencies[dimN-i-1]*1.e10/LIGHTSPEED);
	    }
	    tmp = sqrt(tmp/(double)dimN);
	    printf("\nThe averaged deviation between Real Space and Exact Q is %.6lf THz\n\n", tmp);
	}
    }

        printf ("\nFrequencies in Gamma point without & with NA term) \n\n");
        if (Gammafile!=0) {
            double *energyH = new double[nTHz];
            LatDYN::SrealFij(Gij, natomP, energyH);
            printf ("     Regular    High Accuracy with NA\n");
            for (int i=0; i<nTHz; i++)
                printf ("%4d %10.6lf  %10.6lf  %10.6lf THz\n", i, energyG[i], energyH[i], energyGB[i]);
        } else {
        printf ("      without         with         without    with     LOTO        IR    irrep\n");
        printf ("              (THz)                     (cm-1)\n");
	    int *LOs = new int[nTHz];
	    int *LOt = new int[nTHz];
	    for (int i=0; i<nTHz; i++) {
		LOs[i] = 0;
		LOt[i] = 0;
	    }
	    for (int i=0; i<nTHz; i++) {
	        for (j=0; j<nTHz; j++) {
		    if (LOt[j]) continue;
	    	    if (fabs(energyG[i]-energyGB[j]) >THRE) continue;
		    LOs[i] = 1; 
		    LOt[j] = 1;
		    break;
		}
	    }

	    int jTO = nTHz-1, j;
	    int jLO = nTHz-1;

	    /* for rough analysis (IR/Raman, LO/To) of phonons at Gamma point */
	    for (int i=0; i<nTHz; i++) {
		const char *LOi = "    ";
		const char *LOj = "    ";
		if (!LOs[i]) {
		    LOi = "?TO?";
		    LOj = "?LO?";
            	    for (; jLO >=0; jLO--) if (!LOt[jLO]) break;
		    j = jLO;
		    jLO--;
		} else {
            	    for (; jTO >=0; jTO--) if (LOt[jTO]) break;
		    j = jTO;
		    jTO--;
		}
                printf ("%4d %8.4lf %s %8.4lf %s  %8.2lf %8.2lf  %8.4lf %8.4lf", 
		    i, energyG[i], LOi, energyGB[j], LOj, 
                    energyG[i]*1.e10/LIGHTSPEED, energyGB[j]*1.e10/LIGHTSPEED,
		    ILOTO[j], IRed[j]);
		    if (irrep[j]) {
			for (int r=0; r<Nrep; r++) 
			    if (!strcmp(irrep[j], gp[r]->irrep)) printf("   %4s  %s (%.4lf)\n", irrep[j],  gp[r]->getmode(), EEss[j]);
		    } else  printf("\n");
	    }
	}
        printf ("\n");
        fflush(stdout);

	FILE *cpq = 0;
	if (PQ==1) {
	   cpq = fopen("frequency.dat", "w"); //file saving the calculated phonon frequencies
           fprintf(cpq, "%d %d %d\n", kx*2, ky*2, kz*2);
           for (i=0; i<3; i++) println(cpq, "%18.12f  ", LatDYN::qvec[i]);
	}
	if (pdisfile==0) {
	    // if no dispersion file defined, calculate the phonon DOS
            int freqN = kx*ky*kz*8*natomP*3;
	    LatDYN::setNFrq(kx*ky*kz*8*nTHz);
	    double fx = .5e0/(double)kx;
	    double fy = .5e0/(double)ky;
	    double fz = .5e0/(double)kz;
	    for (i=-kx; i<kx; i++) {
	    for (j=-ky; j<ky; j++) {
	    for (k=-kz; k<kz; k++) {
	        LatDYN *DM = new LatDYN((double)i*fx, (double)j*fy, (double)k*fz);
                if (calcNA!=0)
                    DM->BornGamma(atomP, natomP, dielec, BCharg, mode, Kcell, dirE);
                if (CTLDM) DM->ReduceDM(natomS, MM, NX, FACTOR, RedDM, tranI, calcNA);
                else DM->NonReduceDM(natomS, MM, NX, FACTOR, tranI, calcNA, sigma, sqs);
	    	if (pvdos==0) DM->solveDM();
	    	else DM->solveDMPVDOS(atomP, natomP);
		if (PQ==1)  DM->printPQ(cpq);
	    }
	    }
	    	if (!silent) {
		    printf("%12d out of %d Frequencied calculated! ", 
			LatDYN::Nfreq, freqN);
	            cputim->elptime(stdout);
		}
	    }

	    if (!silent) printf("\nCalculating PDOS from the %d frequencies\n\n", LatDYN::Nfreq);

	    printf ("\n%d frequencies\n\n", LatDYN::Nfreq);
	    LatDYN::printQ(print);
	    LatDYN::printDOS(pvdos, natomP, atomP, nsym, 10, NEDOS, DebCut);
	    FILE *plt;
		if (pvdos) plt = fopen("pvdos.plt", "w"); //partial phonon DOS
		else plt = fopen("vdos.plt", "w"); //phonon DOS
	    
		fprintf(plt, "reset\n"); //making gnuplot script
		if (eps) fprintf(plt, "set terminal postscript landscape enhanced color \"Times_Roman\" 20\n");
		else fprintf(plt, "#set terminal postscript landscape enhanced color \"Times_Roman\" 20\n");
	    fprintf(plt, "set encoding iso_8859_1\n");
	    fprintf(plt, "set pointsize 1.2\n");
	    fprintf(plt, "set size 0.95,0.95\n");
	    if (pvdos) 
		{
			if (eps) fprintf(plt, "\nset output \"pvdos.eps\"\n\n");
			else fprintf(plt, "\n#set output \"pvdos.eps\"\n\n");
		} 
		else
		{
			if (eps) fprintf(plt, "\nset output \"vdos.eps\"\n\n");
			else fprintf(plt, "\n#set output \"vdos.eps\"\n\n");
		}

		fprintf(plt, "funit=%lf\n", funit);
	    if (plot==1) fprintf(plt, "set xlabel \"Frequency (THz)\"\n");
	    else if (plot==2) fprintf(plt, "set xlabel \"Frequency (cm)\"\n");
	    else if (plot==3) fprintf(plt, "set xlabel \"Frequency (meV)\"\n");
	    if (plot==1) fprintf(plt, "set ylabel \"(THz^{-1})\"\n");
	    else if (plot==2) fprintf(plt, "set ylabel \"(cm^{-1})\"\n");
	    else if (plot==3) fprintf(plt, "set ylabel \"(meV^{-1})\"\n");
	    if (!strcmp(expt, ".dat"))
		{
		/* if no experimental data provided */
			if (!pvdos) 
				fprintf(plt, "plot 'vdos.out' using ($1*funit*1.e-12):($2/funit*1.e12) notitle w l lt -1\n\n");
			else
				fprintf(plt, "plot 'pvdos.out' using ($1*funit):($2*$3/funit) notitle w l lt -1\n\n");
		}
	    else 
		{
		/* if experimental data are provided */
			if (!pvdos) {
				fprintf(plt, "plot 'vdos.out' using ($1*funit*1.e-12):($2/funit*1.e12) notitle w l lt -1, \\\n");
				fprintf(plt, "     '%s' using ($1*funit):($2/funit) notitle w lp pt 6 lt 1\n", expt);
			} else {
				fprintf(plt, "plot 'pvdos.out' using ($1*funit):($2*$3/funit) notitle w l lt -1, \\\n");
				fprintf(plt, "     '%s' using ($1*funit):($2*funit) notitle w lp pt 6 lt 1\n", expt);
			}
	    }

		//fprintf(plt, "\nexit\n");
	    fclose(plt);
	    if (plot) ynuplot ("vdos.plt");
	} else {
		int exptUsed = 0;
	    /* for dispersion calculation
               making the gnuplot script */
	    FILE *fdis = fopen(pdisfile, "r");
	    FILE *fout = fopen("vdis.out", "w");
	    FILE *headplt = fopen("headplt.tmp", "w+"); //gnuplot head file
	    FILE *midplt = fopen("midplt.tmp", "w+"); //gnuplot body file
	    FILE *tailplt = fopen("tailplt.tmp", "w+"); //gnuplot tail file
	    FILE *exptplt = fopen("exptplt.tmp", "w+"); //expt file
	    FILE *vline = fopen("vline.dat", "w+"); // vertical lines for dispersion plot
	    double Dir[3], ToGamma[3], tmp[3], db[3], disp = 0.e0, fmax=0.e0, fmin=0.e0;
	    int dcount = 0;
	    fprintf(headplt, "qunit=1.0\n");
	    fprintf(headplt, "eunit=%lf\n", eunit);
	    fprintf(headplt, "funit=%lf\n", funit);

	    int lineN=0, *exptI=0, exptN, *exptC=0;
	    char *str6=0, *str7=0;
	    char **exptS=0;

	    double fup = fmaxT;

	    /* processing data fields in the dispersion file */
	    while (1) {
              yetline(&line, &len, fdis);
	      if (feof(fdis)) break;
	      sscanf (line, "%lf %lf %lf %lf %lf %lf", a, a+1, a+2, b, b+1, b+2); //starting and end Q point
	      SplitItem *pos = new SplitItem(line);
	      lineN = pos->GetN();
	      if (lineN < 6 ) continue;
	      if (lineN>6) {
		str6 = pos->str(6); //the label of the starting high symmetry point
		if (!strcmp(str6,"Gamma")) str6 = (char *)"{/Symbol G}";
	      }
	      if (lineN>7) {
		str7 = pos->str(7); //the label of the end high symmetry point
		if (!strcmp(str7,"Gamma")) str7 = (char *)"{/Symbol G}";
	      }

	      exptN = (lineN - 8)/3;
	      if (exptN > 0) {
			  exptUsed = 1;
		/* when experimental data are provided */
		exptI = new int[exptN];
		exptC = new int[exptN];
		exptS = new char* [exptN];
		for (int i=0; i<exptN; i++) {
		  exptI[i] = atoi(pos->str(7+i*3+1)); //index of data in the expt file
		  exptS[i] = pos->str(7+i*3+2); //column for Q point in the expt file
		  exptC[i] = atoi(pos->str(7+i*3+3)); //column frequency value
		}
	      }

	      //prepare the dipersion calculation

	      LatDYN::setNFrq(NdisQ*nTHz);
//printf ("NdisQ*nTHz=%d\n", NdisQ*nTHz);

              for (int ii=0; ii<3; ii++) Dir[ii] = b[ii] - a[ii]; //dispersion direction

	      for (i=0; i<NdisQ; i++) {
		double fac = (double)i/((double)(NdisQ-1));
		double tG[3], tG1[3];
		for (int ii=0; ii<3; ii++) {
		   db[ii] = Dir[ii]*fac;
		   tG[ii] = a[ii]+db[ii]; //direction is always toward to the Gamma point
		}

                mproduct (tG, lattice, ToGamma);
	        LatDYN *DM = new LatDYN(tG[0], tG[1], tG[2], lattice);

                mproduct (Dir, lattice, tG1);
		//if (normal(ToGamma)<THR) 
		    for (int ii=0; ii<3; ii++) ToGamma[ii] = tG1[ii];
		if (calcNA!=0) DM->BornGamma(
		  atomP, natomP, dielec, BCharg, mode, Kcell, ToGamma);
                if (CTLDM) DM->ReduceDM(natomS, MM, NX, FACTOR, RedDM, tranI, calcNA);
                else DM->NonReduceDM(natomS, MM, NX, FACTOR, tranI, calcNA, sigma, sqs);
	    	DM->solveDM();
                mproduct (db, lattice, tmp);
		double ftmp;
	        fmax = max(fmax,DM->printDIS(fout, i, tmp, &ftmp));
		fmin = min(fmin, ftmp);
	      }

	      if (fup <= 0.e0) fup = fmax;

	      if (!silent) printf("Branch %.3f %.3f %.3f to %.3f %.3f %.3f done!\n\n",
		  a[0], a[1], a[2], b[0], b[1], b[2]);
	      fprintf(fout, "\n\n");

              mproduct (Dir, lattice, ToGamma);
              mproduct (ToGamma, LatDYN::qvec, tmp);
	      double dd = normal(tmp);
		/* write to the gnuplot script */
	      fprintf(headplt, "p%d = %lf\npp%d = %lf\n", dcount, disp, dcount, dd);
	      if (dcount==0) {
		if (lineN==6) 
		  fprintf(midplt, "set xtics ( 'K%d' qunit*%lf", dcount,  disp);
		if (lineN>6) 
		  fprintf(midplt, "set xtics ( '%s' qunit*%lf", str6, disp);
	      } else {
		if (lineN==6) 
		  fprintf(midplt, ", 'K%d' qunit*%lf", dcount, disp);
		if (lineN>6) 
		  fprintf(midplt, ", '%s' qunit*%lf", str6, disp);
	      }
	      disp += dd;

		/* determine the vertical lines in the phonon dispersion plot */
	      for (double f=-fup; f<1.9*fup; f+=0.1e0*fup) 
	      	fprintf(vline, "%lf %lf\n", disp, f);
	      fprintf(vline, "\n\n");

	      //fprintf(vline, "%lf %lf\n%lf 0\n%lf %lf\n\n\n", disp, -fmax, disp, disp, fmax*1.99);

	      if (dcount==0) {
		/* first writing, include the file name */
		fprintf(tailplt, " 'vdis.out' index %d using (qunit*p%d+qunit*$1):(funit*$5) notitle w l lt -1", dcount, dcount);
		if (exptN>0) {
		  fprintf(exptplt, " '%s' index %d using (qunit*p%d+qunit*pp%d*%s):(eunit*$%d) notitle w p pt 6 lt 1", expt, exptI[0], dcount, dcount, exptS[0], exptC[0]);
		  for (int ii=1; ii<exptN; ii++) 
		    fprintf(exptplt, ", \\\n '' index %d using (qunit*p%d+qunit*pp%d*%s):(eunit*$%d) notitle w p pt %d lt %d", 
		      exptI[ii], dcount, dcount, exptS[ii], exptC[ii], ii+6, ii+1);
		}
	      } else {
		/* continue writing, file name is blanked follow gnuplot */
		fprintf(tailplt, ", \\\n '' index %d using (qunit*p%d+qunit*$1):(funit*$5) notitle w l lt -1", dcount, dcount);
		if (ftell(exptplt)==0) {
//printf("ftell = %d\n", ftell(exptplt));
		    if (exptN>0) fprintf(exptplt, " '%s' index %d using (qunit*p%d+qunit*pp%d*%s):(eunit*$%d) notitle w p pt 6 lt 1", expt, exptI[0], dcount, dcount, exptS[0], exptC[0]);
		for (int ii=1; ii<exptN; ii++) 
		  fprintf(exptplt, ", \\\n '' index %d using (qunit*p%d+qunit*pp%d*%s):(eunit*$%d) notitle w p pt %d lt %d", 
		    exptI[ii], dcount, dcount, exptS[ii], exptC[ii], ii+6, ii+1);
		} else {
		  for (int ii=0; ii<exptN; ii++) 
		    fprintf(exptplt, ", \\\n '' index %d using (qunit*p%d+qunit*pp%d*%s):(eunit*$%d) notitle w p pt %d lt %d", 
		    exptI[ii], dcount, dcount, exptS[ii], exptC[ii], ii+6, ii+1);
		}
	      }

	      for (int ii=1; ii<nTHz; ii++) 
		fprintf(tailplt, ", \\\n '' index %d using (qunit*p%d+qunit*$1):(funit*$%d) notitle w l lt -1", dcount, dcount, ii+5);
	      dcount++;
	    }
	    fprintf(headplt, "p%d = %lf\n\n", dcount, disp);
	    if (lineN>6 && lineN<8) { fprintf(midplt, ", 'K%d' qunit*%lf)\n\n", dcount, disp);
		/* if not label given, I give it */
		printf("lineN=%d, 'K%d' qunit*%lf)\n\n", lineN, dcount, disp);
	    } else fprintf(midplt, ", '%s' qunit*%lf)\n\n", str7, disp);
	    //if (!strcmp(expt, ".dat")) fprintf(tailplt, "\n\n");
	    if (!exptUsed || !strcmp(expt, ".dat")) fprintf(tailplt, "\n\n");
	    else fprintf(tailplt, ", \\\n");
	    fprintf(exptplt, "\n\n");
	/* rewind file */
	    fseek ( headplt , 0L , SEEK_SET );
	    fseek ( midplt , 0L , SEEK_SET );
	    fseek ( tailplt , 0L , SEEK_SET );
	    fseek ( exptplt , 0L , SEEK_SET );

	/* finalizing the gnuplot script */
	    FILE *plt = fopen("vdis.plt", "w");
	    fprintf(plt, "reset\n");
		if (eps) fprintf(plt, "set terminal postscript landscape enhanced color \"Times_Roman\" 20\n");
		else fprintf(plt, "#set terminal postscript landscape enhanced color \"Times_Roman\" 20\n");
	    fprintf(plt, "set encoding iso_8859_1\n");
	    fprintf(plt, "set pointsize 1.2\n");
	    fprintf(plt, "set size 0.95,0.95\n");
	    if (eps) fprintf(plt, "\nset output \"vdis.eps\"\n\n");
		else fprintf(plt, "\n#set output \"vdis.eps\"\n\n");
            while (yetline(&line, &len, headplt)!=-1) fprintf(plt, "%s", line);
            while (yetline(&line, &len, midplt)!=-1) fprintf(plt, "%s", line);
	    fprintf(plt, "set key left top\n");
	    if (plot==1) fprintf(plt, "set ylabel \"Frequency (THz)\"\n\n");
	    else if (plot==2) fprintf(plt, "set ylabel \"Frequency (cm)\"\n\n");
	    else if (plot==3) fprintf(plt, "set ylabel \"Frequency (meV)\"\n\n");
	    if (fmin > -0.01) fmin = 0.e0;
	    else {
		fmin = floor(fmin*1.1);
	      	fprintf(vline, "0.0 0.0\n%lf 0.0\n", disp);
	    }
	    fprintf(plt, "plot [x=0:qunit*p%d*1.0001] [funit*%lf:funit*%lf] \\\n", dcount,
		fmin, floor(fup*1.1+1));
	    fprintf(plt, "'vline.dat' using (qunit*$1):(funit*$2) notitle w l lt 4, \\\n");
            while (yetline(&line, &len, tailplt)!=-1) fprintf(plt, "%s", line);
            if (!exptUsed || !strcmp(expt, ".dat")) while (yetline(&line, &len, exptplt)!=-1) fprintf(plt, "#%s", line);
            else while (yetline(&line, &len, exptplt)!=-1) fprintf(plt, "%s", line);

	    fclose(headplt);
	    fclose(midplt);
	    fclose(tailplt);
	    fclose(exptplt);
	    fclose(vline);
	    fclose(plt);
	    remove("headplt.tmp");
	    remove("midplt.tmp");
	    remove("tailplt.tmp");
	    remove("exptplt.tmp");
	    fclose(fout);
	    if (plot) ynuplot("vdis.plt");
	}

	if (PQ==1) fclose(cpq);

	/* print out the reciprocal lettice vectors, helping preparing the the dispersion file */
        printf("\nPrimmitive Q/PI2 vectors in the Inverse space\n\n");
        for (i=0; i<3; i++) println("%18.12f  ",LatDYN::qvec[i]);
}
