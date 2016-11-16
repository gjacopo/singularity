/*      Program Dh_evaluation.c. Version: 2 de Diciembre, 2005   */
 
/*      Directorio de grabacion: 10.MULTIFRACTAL      */
/*      Librerias: LIB-1.2                            */
 
 
#include <stdio.h>
#include <math.h>
#include <string.h>
 
                               

//      My libraries

#include <multifractal.c>  // Version del  24 de Agosto, 2004


/*      Global variables               */

/*      General purpose variables      */

int VERBOSE=0;     // Flag controlling the verbose information

/*        Variables defining the series to be processed   */

int NSERIES=1;     // Number of output series
int D_space=1;     // Dimension of the space

int LEFF=512;      // Length/size of the series/images
int FROM_DH=0;     // Flag. If activated, the program retrieves the D(h)
                   // from previously recorded files, then calculates errors

int GEO_MAP=0;     // Flag. If activated, a map of qualities for series of the
                   // defined type and varying geometries is created
                   // (see parameters defining its geometry below)

int TYPE_MAP=0;    // Flag. If activated, a map of fixed geometry (1x16384)
                   // and varying parameters, according to the type selected,
                   // is generated (see the parameter choice below)
int TYPE=0;        // Type of multifractal to be processed.
                   // 0: Log-Poisson
                   // 1: Log-Normal
                   // 2: Log-Levi
                   // 3: Binomial
float HINF=-0.5;   // Most singular exponent (Types 0 and 3)
float CODINF=1.;   // Most singular co-dimension Type 0)
float H1=0.5;      // Upper singularity bound in binomial model
float MU=0.5;      // Average singularity (for log-Normal and log-Levi)
float SIGMA=1.;    // Singularity dispersion in log-Normal and log-Levi 
float ALPHA=1.5;   // Exponent defining log-Levi


/*      Variables concerning the analyzing techniques   */

int METHOD=15; // Composite variable, defining the method/s to be used 
              // in the analysis
              // 1*GH+2*GMWP+4*MOMENTS+8*WTMM
              // WTMM is called via LastWave, if LASTWAVE flag is set


/*          Histogram based methods        */

int NBOX=128;       // Initial number of bins in the histogram

/*          WTMM method                    */

/* Define the method used for estimating the Legendre transform */
int  DIRECT=1;  // Flag. If enabled, calculation is done by Legendre transform
#define CANONICAL (1-DIRECT)
int SUP_WTMM=0; // Flag. If given, WTMM is evaluated over suprema

int wav_wtmm=1;      /* Default choice for wavelet: gaussian */
int ord_der_wtmm=3;  // Derivation order for the wavelet in WTMM

float MINSCALE=1.;   

float ScRatio=(1./16.); 
float ScStep=2.; //Scale step

#define SIGN(x) (((x)>0.)?(1):(-1))

/*         Last Wave parameters           */

int LASTWAVE=0; // Flag. If set, WTMM is done via external calls to LastWave
int LW_AMIN=1.5;
int LW_OMAX=3;

#define LWCALL "LastWave"  // Put here the command name for the LastWave in
                           // your system


/*           Moment and WTMM methods    */

#define Nmom 65

double moms[Nmom]={-4.0,-3.6,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.3,2.6,3.0,3.5,4.0,5.0,6.0,7.0,8.0};

#define Ndist 10
int dist[Ndist]={4,5,6,8,10,12,15,18,22,30};


/*          Defining geometry for GEO_MAP    */

#define Nnums 5
#define Nleffs 4

int nums[Nnums]={1,10,100,1000,10000};
int leffs[Nleffs]={1024,4096,16384,65536};

/*          Defining parameters for TYPE_MAP    */

#define NLPs 10

const double hinfs[NLPs]=   {-.8, -.6, -.4, -.2, -.6, -.4, -.2, -.4, -.2, -.2};
const double codinfs[NLPs]={ 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.6, 0.6, 0.4};

#define Nmeans 5
#define Nsigmas 4

const double means[Nmeans]= {-1.,-0.5,0,0.5,1.};
const double sigmas[Nsigmas]={0.33,0.5,1.,2.};


/*       Internal global variables                  */

double shiftw;
double shiftg;

/*       Other parameters (non-externally accessed)  */

float DH=0.1;      // Accepted conventional error



/*       Definitions

/*         Quality struct   */

typedef struct
{
    double g;  // GH quality
    double w;  // GMWP quality
    double m;  // Moments quality
    double l;  // WTMM  quality
} Quality;

/*      Function prototypes         */
 
int main(int argc, char *argv[]);
void parse_arguments(int argc, char *argv[]);


void genera_base( int verbose, int leff, char *base);
Quality analiza_series( char *path);
double media_exp( int dimx, double sc0, double **expon);
void estima_Dh_g( char *path, char *base, int *Nh_g, double **h_g,
		  double **Dh_g, double **errDh_g);
void estima_Dh_w( char *path, char *base, int *Nh_w, double **h_w,
		  double **Dh_w, double **errDh_w);
void estima_Dh_mom( char *path, char *base, int *Nh_m, double **h_m,
		    double **Dh_m, double **errDh_m);
void estima_Dh_wtmm( char *path, char *base, int *Nh_wtmm, double **h_wtmm,
		     double **Dh_wtmm, double **errDh_wtmm);
void estima_Dh_lw( char *path, char *base, int *Nh_l, double **h_l,
		   double **Dh_l, double **errDh_l);


void WTMM_PF_compute( double *signal, int dimx, double expon, int order,
		      int nq, double *qArray, int nsc, double *scArray, 
		      double **ExtWTlis, double **Z, int *n_ext );

void WT_transform( double *signal, int dimx, double expon, int order,
		   double sc,  double *wtrans );

void define_wavelet( int dimx, double sc, double expon, int order, 
			double *wavesig);

void W_projection( double *signal,  double *filt, int dimx, double sc, 
		   double *wtrans );

double comp_coefficient(int n, double t);

int WTMM_find_extrema( double *wtrans, int *wtrans_ind, int dimx, double sc);
  
int PF_extrema_track( double *wtrans, int *wtrans_ind, 
		      double *maxsig, int *maxsig_ind,
		      int jj, int  n_max,
		      int nq, double *qArray, int nsc, double *scArray, 
		      double *ExtWTlis, double *Z );

int linefit( double *yValues, double *xValues, int N, double *slope );

int compare(const double *d1,const double *d2);

void directSpecCompute( double **Z, int N, 
			int nq, double *qArray, int nws, double *wsArray,
			double *tauq, double *h, double *Dh );

void canonMeanCompute( double **ExtWTlis, int N, int *n_ext,
		       int nq, double *qArray, int nws, double *wsArray,
		       double **sTq, double **sTqLogT /* double **logSTq */ );

void canonSpecCompute( double **sTq, double **sTqLogT, /* double **logSTq */
		       int N, int nq, double *qArray, int nws, double *wsArray,
		       double *Tauq, double *H, double *Dh);


void acumula_momentos( double **signal, double **moments);
void calcula_taup( double **moments, double *taup);
int simple_Legendre_transform( double *taup, double **h_m, double **Dh_m,
				double **errDh_m);

void genera_expon( int leff, double **signal, double **expon);
void genera_expon_1D( int leff, double *serie, double *expon);
void genera_expon_2D( int leff, double **image, double **expon);

void acumula_histograma( int leff, double *mh, double **expon, 
			 double *histo);
int calcula_Dh_histo( double sc0, double *mh, double *histo, double **h, 
		      double **Dh, double **errDh);
int calcula_Dh_histo_no_ponderado( double sc0, double *mh, double *histo, 
				   double **h, double **Dh, double **errDh);
int carga_Dh( char *nombre, int *Nh, double **h, double **Dh, double **errDh);
void lee_serie_temp( char *name_in, int *dims, int *dimt, double ***series);
int columnas_serie_temp( FILE *canal);
int lineas_serie_temp( FILE *canal);

double registra_Dh( int graba, int Nr,char *nombre, double sc0, double *h, 
		    double *Dh, double *errDh);
double theoretical_Dh( double h);
void theoretical_Deltah( double *hmin, double *hmax);
void lee_datos( int leff, char *nombre_in, double **signal);
void graba_datos( int leff, char *nombre_in, double **signal);
void lee_datos_float( int leff, char *nombre_in, double **signal);
void graba_datos_float( int leff, char *nombre_in, double **signal);
void graba_serie( int leff, char *nombre_in, double *datos);
void graba_imagen( int leff, char *nombre_in, double **datos);

void graba_geomapa( char *nombre, double **Q);
void graba_typemap( char *nombre, double **Q);

double moda(int leff, double **signal);
double moda_1D(int dimx, double *datos);
double moda_2D(int dimx, int dimy, double **datos);
double moda_por_histo( int Nhisto, double *mm, double *histo);

int main(int argc, char *argv[])
{

/*	PARAMETERS		*/


/*	DATA			*/

	char path[90],base[90];

/*	QUANTITIES TO BE COMPUTED	*/

	double **Qg,**Qw,**Qm,**Ql;

	
/*	AUXILIAR VARIABLES	*/

	Quality Qtot;
	char nombre[90];
	char GHname[90],GMWPname[90],Momname[90],WTMMname[90];
	int ileffs,inums;
	int iLPs,imeans,isigmas;

/*		Program		*/


	parse_arguments(argc,argv);

/*        Re-fitting parameters      */
	
/*
     To avoid problems, the parameter DERIVA_MODO is shortcircuited  
     to a plain rightwards finite difference
*/

	DERIVA_MODO=1;

/*     LEFF is converted to the minimum power of 2 larger than it   */

	LEFF=dimensiona(LEFF);

/*     LASTWAVE option only available for 1D signals    */

	if(D_space==2) METHOD=Mod(METHOD,8);


/*        Inputting the path to the data      */

	if(VERBOSE) printf("Please, give me the path to the data files\n");
	scanf("%s",path);


/*       Generating the base for the MF names     */

	if(VERBOSE)
	    printf("Processing for %d output %dD signals\n",NSERIES,D_space);

	

/*       Analyzing series with the methods      */

	if(GEO_MAP+TYPE_MAP==0) analiza_series(path);

	if(GEO_MAP)
	{

/*        Initialization    */

	    genera_base(0,0,base);
	    if(METHOD&0x01) 
	    {
		Qg=reservar_matriz(Nnums,Nleffs);
		sprintf(GHname,"GeoMap-GH-%s.dat",base);
	    }
	    if(METHOD&0x02)
	    {
		Qw=reservar_matriz(Nnums,Nleffs);
		sprintf(GMWPname,"GeoMap-GMWP-%s.dat",base);
	    }
	    if(METHOD&0x04)
	    {
		Qm=reservar_matriz(Nnums,Nleffs);
		sprintf(Momname,"GeoMap-Mom-%s.dat",base);
	    }
	    if(METHOD&0x08)
	    {
		Ql=reservar_matriz(Nnums,Nleffs);
		sprintf(WTMMname,"GeoMap-WTMM-%s.dat",base);
	    }

/*        Obtaining the quality matrices       */

	    for(inums=0;inums<Nnums;inums++)
	    {
		NSERIES=nums[inums];
		for(ileffs=0;ileffs<Nleffs;ileffs++)
		{
		    LEFF=leffs[ileffs];
		    Qtot=analiza_series(path);
		    if(METHOD&0x01)
		    {
			Qg[inums][ileffs]=Qtot.g;
			graba_geomapa(GHname,Qg);
		    }
		    if(METHOD&0x02)
		    {
			Qw[inums][ileffs]=Qtot.w;
			graba_geomapa(GMWPname,Qw);
		    }
		    if(METHOD&0x04)
		    {
			Qm[inums][ileffs]=Qtot.m;
			graba_geomapa(Momname,Qm);
		    }
		    if(METHOD&0x08) 
		    {
			Ql[inums][ileffs]=Qtot.l;
			graba_geomapa(WTMMname,Ql);
		    }
		}
		
	    }

/*       Recording the results         */


/*      Memory release before end    */

	    
	    if(METHOD&0x01) liberar_matriz(Qg,Nnums);
	    if(METHOD&0x02) liberar_matriz(Qw,Nnums);
	    if(METHOD&0x04) liberar_matriz(Qm,Nnums);
	    if(METHOD&0x08) liberar_matriz(Ql,Nnums);
	}
	
	if(TYPE_MAP)
	{
	    NSERIES=1;
	    LEFF=65536;
	    switch(TYPE)
	    {
		case 0:

/*        Memory reservation for the quality matrices    */

		    if(METHOD&0x01) Qg=reservar_matriz(1,NLPs);
		    if(METHOD&0x02) Qw=reservar_matriz(1,NLPs);
		    if(METHOD&0x04) Qm=reservar_matriz(1,NLPs);
		    if(METHOD&0x08) Ql=reservar_matriz(1,NLPs);

/*        Obtaining the quality matrices       */

		    for(iLPs=0;iLPs<NLPs;iLPs++)
		    {
			HINF=hinfs[iLPs];
			CODINF=codinfs[iLPs];
			Qtot=analiza_series(path);
			if(METHOD&0x01) Qg[0][iLPs]=Qtot.g;
			if(METHOD&0x02) Qw[0][iLPs]=Qtot.w;
			if(METHOD&0x04) Qm[0][iLPs]=Qtot.m;
			if(METHOD&0x08) Ql[0][iLPs]=Qtot.l;
		    }

/*       Recording the results         */
		    
		    sprintf(base,"Log-Poisson");
		    if(D_space==1) strcat(base,"_1D");
		    else strcat(base,"_2D");
		    if(METHOD&0x01)
		    {
			sprintf(nombre,"TypeMap-GH-%s.dat",base);
			graba_typemap(nombre,Qg);
		    }
		    if(METHOD&0x02)
		    {
			sprintf(nombre,"TypeMap-GMWP-%s.dat",base);
			graba_typemap(nombre,Qw);
		    }
		    if(METHOD&0x04)
		    {
			sprintf(nombre,"TypeMap-Mom-%s.dat",base);
			graba_typemap(nombre,Qm);
		    }
		    if(METHOD&0x08)
		    {
			sprintf(nombre,"TypeMap-WTMM-%s.dat",base);
			graba_typemap(nombre,Ql);
		    }


/*      Memory release before end    */

	    
		    if(METHOD&0x01) liberar_matriz(Qg,1);
		    if(METHOD&0x02) liberar_matriz(Qw,1);
		    if(METHOD&0x04) liberar_matriz(Qm,1);
		    if(METHOD&0x08) liberar_matriz(Ql,1);
		    break;

		case 1:

/*        Memory reservation for the quality matrices    */

		    if(METHOD&0x01) Qg=reservar_matriz(Nmeans,Nsigmas);
		    if(METHOD&0x02) Qw=reservar_matriz(Nmeans,Nsigmas);
		    if(METHOD&0x04) Qm=reservar_matriz(Nmeans,Nsigmas);
		    if(METHOD&0x08) Ql=reservar_matriz(Nmeans,Nsigmas);

/*        Obtaining the quality matrices       */

		    for(imeans=0;imeans<Nmeans;imeans++)
		    {
			MU=means[imeans];
			for(isigmas=0;isigmas<Nsigmas;isigmas++)
			{
			    SIGMA=sigmas[isigmas];
			    Qtot=analiza_series(path);
			    if(METHOD&0x01) Qg[imeans][isigmas]=Qtot.g;
			    if(METHOD&0x02) Qw[imeans][isigmas]=Qtot.w;
			    if(METHOD&0x04) Qm[imeans][isigmas]=Qtot.m;
			    if(METHOD&0x08) Ql[imeans][isigmas]=Qtot.l;
			}
		    }

/*       Recording the results         */

		    sprintf(base,"Log-Normal");
		    if(D_space==1) strcat(base,"_1D");
		    else strcat(base,"_2D");
		    if(METHOD&0x01)
		    {
			sprintf(nombre,"TypeMap-GH-%s.dat",base);
			graba_typemap(nombre,Qg);
		    }
		    if(METHOD&0x02)
		    {
			sprintf(nombre,"TypeMap-GMWP-%s.dat",base);
			graba_typemap(nombre,Qw);
		    }
		    if(METHOD&0x04)
		    {
			sprintf(nombre,"TypeMap-Mom-%s.dat",base);
			graba_typemap(nombre,Qm);
		    }
		    if(METHOD&0x08)
		    {
			sprintf(nombre,"TypeMap-WTMM-%s.dat",base);
			graba_typemap(nombre,Ql);
		    }

/*      Memory release before end    */

	    
		    if(METHOD&0x01) liberar_matriz(Qg,Nmeans);
		    if(METHOD&0x02) liberar_matriz(Qw,Nmeans);
		    if(METHOD&0x04) liberar_matriz(Qm,Nmeans);
		    if(METHOD&0x08) liberar_matriz(Ql,Nmeans);
		    break;

		default:
		    printf("Map not implemented for this type, sorry!\n");
		    break;
	    }


	}


/*	FREEING MEMORY BEFORE FINISHING	*/


	exit(0);

}
	
void parse_arguments(int argc, char *argv[])
{
	char **olarg,**olval,**olexp;
	float **ptrvar_f,**ptrval_f;
	int **ptrvar_i,**ptrval_i;
	int **ptrflag;
	int *type;
        int lar;
	int flagv;
	int arglen;
	int Narg0=50;  // Initialization value; it should be greater than (but not necessarily equal to) the number or arguments
	int in,Narg;

/*         Defining the expected arguments            */
	
	olarg=reservar_matriz_char(Narg0,20);   // Name of the variable 
	olval=reservar_matriz_char(Narg0,20);   // Name of the argument
	olexp=reservar_matriz_char(Narg0,300);  // Explanation
	ptrvar_f=(float **) calloc(Narg0,sizeof(float *)); // Pointer to the variable, if it is float
	ptrval_f=reservar_matriz_float(Narg0,2); // Allowed min and max
	ptrvar_i=(int **) calloc(Narg0,sizeof(int *)); // Pointer to the variable, if int
	ptrval_i=reservar_matriz_int(Narg0,2); // Allowed min and max
	ptrflag=(int **) calloc(Narg0,sizeof(int *)); // Pointer to the flag
	type=(int *) calloc(Narg0,sizeof(int)); // Type of variable

	in=0;
	
//    Argument VERBOSE. Type 0: flag

	sprintf(olarg[in],"%s","-ver");
	sprintf(olexp[in],"%s\n %s : %s\n",
		"\nGENERAL PARAMETERS\n==================",
		olarg[in],
		"Flag. If enabled, the program shows lots of (verbose) information,\nspecially for multifractal analysis. Default: DISABLED"); 
	type[in]=0;
	ptrflag[in]=&VERBOSE;
	in++;

//    Argument FROM_DH. Type 0: flag

	sprintf(olarg[in],"%s","-fromDh");
	sprintf(olexp[in]," %s : %s\n",
		olarg[in],
		"Flag. If enabled, the program takes previously computed D(h) files\nand estimates the error from them.\n Default: DISABLED"); 
	type[in]=0;
	ptrflag[in]=&FROM_DH;
	in++;

//    Argument GEO_MAP. Type 0: flag

	sprintf(olarg[in],"%s","-geomap");
	sprintf(olexp[in]," %s : %s\n",
		olarg[in],
		"Flag. If enabled, the program tries to generate quality maps\nfor each method, changing geometry but keeping the given MF type.\n Default: DISABLED"); 
	type[in]=0;
	ptrflag[in]=&GEO_MAP;
	in++;

//    Argument TYPE_MAP. Type 0: flag

	sprintf(olarg[in],"%s","-typemap");
	sprintf(olexp[in]," %s : %s\n",
		olarg[in],
		"Flag. If enabled, the program tries to generate quality maps\nfor each method, for fixed geometry (1x16384) and changing parameters in\nthe given MF type.\n Default: DISABLED"); 
	type[in]=0;
	ptrflag[in]=&TYPE_MAP;
	in++;


//    Argument NSERIES. Type 1: integer

	sprintf(olarg[in],"%s","-N");
	sprintf(olval[in],"%s","#series");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Number of series to be processed. Default:",NSERIES); 
	type[in]=1;
	ptrvar_i[in]=&NSERIES;
	ptrval_i[in][0]=1;
	ptrval_i[in][1]=10000;
	in++;

//    Argument LEFF. Type 1: integer

	sprintf(olarg[in],"%s","-dim");
	sprintf(olval[in],"%s","length");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Size of series to be processed. Default:",LEFF); 
	type[in]=1;
	ptrvar_i[in]=&LEFF;
	ptrval_i[in][0]=256;
	ptrval_i[in][1]=65536;
	in++;



//    Argument D_space. Type 1: integer

	sprintf(olarg[in],"%s","-d_space");
	sprintf(olval[in],"%s","dimension");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Dimension of the embedding space. Default:",
		D_space); 
	type[in]=1;
	ptrvar_i[in]=&D_space;
	ptrval_i[in][0]=1;
	ptrval_i[in][1]=2;
	in++;

//    Argument TYPE. Type 1: integer

	sprintf(olarg[in],"%s","-type");
	sprintf(olval[in],"%s","mult_type");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Type of multifractal to be generated.\n   0: Log-Poisson\n   1: Log-Normal\n   2: Log-Levi\n   3: Binomial\n Default:",TYPE); 
	type[in]=1;
	ptrvar_i[in]=&TYPE;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=3;
	in++;

//    Argument HINF. Type 2: float

	sprintf(olarg[in],"%s","-hinf");
	sprintf(olval[in],"%s","min_sing");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Most singular exponent. Valid for log-Poisson and binomials. Default:",HINF); 
	type[in]=2;
	ptrvar_f[in]=&HINF;
	ptrval_f[in][0]=-1.;
	ptrval_f[in][1]=0.;
	in++;

//    Argument CODINF. Type 2: float

	sprintf(olarg[in],"%s","-Codinf");
	sprintf(olval[in],"%s","min_sing_cod");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Most singular codimension. Valid for log-Poisson and binomials. Default:",CODINF); 
	type[in]=2;
	ptrvar_f[in]=&CODINF;
	ptrval_f[in][0]=0.;
	ptrval_f[in][1]=D_space;
	in++;

//    Argument H1. Type 2: float

	sprintf(olarg[in],"%s","-h1");
	sprintf(olval[in],"%s","max_sing");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Maximum singularity in binomial MFs. Default:",H1); 
	type[in]=2;
	ptrvar_f[in]=&H1;
	ptrval_f[in][0]=-1.;
	ptrval_f[in][1]=2.;
	in++;



//    Argument MU. Type 2: float


	sprintf(olarg[in],"%s","-mean");
	sprintf(olval[in],"%s","sing_av");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Singularity mean. Valid for log-Normal and log-Levi. Default:",
		MU); 
	type[in]=2;
	ptrvar_f[in]=&MU;
	ptrval_f[in][0]=-5.;
	ptrval_f[in][1]=5.;
	in++;

//    Argument SIGMA. Type 2: float

	sprintf(olarg[in],"%s","-sigma");
	sprintf(olval[in],"%s","disp_sing");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Singularity dispersion. Valid for log-Normal and log-Levi.\nDefault:",SIGMA); 
	type[in]=2;
	ptrvar_f[in]=&SIGMA;
	ptrval_f[in][0]=0.;
	ptrval_f[in][1]=5.;
	in++;

//    Argument ALPHA. Type 2: float

	sprintf(olarg[in],"%s","-alpha");
	sprintf(olval[in],"%s","exponent");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Valid for log-Levi only: log-Levi exponent. Default:",
		ALPHA); 
	type[in]=2;
	ptrvar_f[in]=&ALPHA;
	ptrval_f[in][0]=0.;
	ptrval_f[in][1]=2.;
	in++;


//    Memory usage parameters (included in <FFT1D.c>)

	in=parsing_memory(in,0,0,olarg,olval,olexp,ptrvar_f,ptrval_f, 
			  ptrvar_i,ptrval_i,ptrflag,type);

//   multifractal parameters (included in <multifractal_1D.c>)

	in=parsing_multifractal_1D(in,0,0,olarg,olval,olexp,ptrvar_f,ptrval_f,
			  ptrvar_i,ptrval_i,ptrflag,type);



//    Argument NBOX. Type 1: integer

	sprintf(olarg[in],"%s","-Nbin");
	sprintf(olval[in],"%s","#bins");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Number of histogram bins. Default:",NBOX); 
	type[in]=1;
	ptrvar_i[in]=&NBOX;
	ptrval_i[in][0]=2;
	ptrval_i[in][1]=16384;
	in++;


//    Argument METHOD. Type 1: integer

	sprintf(olarg[in],"%s","-Method");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Composite variable expressing the method or methods to be used.\nIts value is 1*GH+2*GMWP+4*MOM+8*WTMM.\n Default:",METHOD); 
	type[in]=1;
	ptrvar_i[in]=&METHOD;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=15;
	in++;

//    Argument LASTWAVE. Type 0: flag

	sprintf(olarg[in],"%s","-LastWave");
	sprintf(olexp[in],"%s\n %s : %s\n",
		"\nWTMM PARAMETERS\n==================",
		olarg[in],
		"Flag. If enabled, WTMM calculations are done by external calls to LastWave.\n Default: DISABLED"); 
	type[in]=0;
	ptrflag[in]=&LASTWAVE;
	in++;

/*             End of parameters definition      */

	Narg=in;

/*                   Parsing loop                         */

	for(lar=1,flagv=1;(lar<argc)&&(flagv==1);lar++)
	{
	  flagv=0;
	  for(in=0;(in<Narg)&&(flagv==0);in++)
	    if( (!strcmp(argv[lar],"-h")) || (!strcmp(argv[lar],"-help")) ) flagv=2;
	    else if(!strcmp(argv[lar],olarg[in])) flagv=1;
	
	  if(flagv==1)
	  {
	    in--;
	    if(type[in]%3==2)
	    {
	      lar++;
	      sscanf(argv[lar],"%f",ptrvar_f[in]);
	      if((*ptrvar_f[in]<ptrval_f[in][0])
		 ||(*ptrvar_f[in]>ptrval_f[in][1]))
	      {
		printf("\nValue out of range (%f - %f) for argument %s\n\n",
		       ptrval_f[in][0],ptrval_f[in][1],olarg[in]);
		flagv=0;
	      }
	      if(type[in]==5) *ptrflag[in]=1;
	    }
	    else if(type[in]%3==1)
	    {
	      lar++;
	      sscanf(argv[lar],"%d",ptrvar_i[in]);
	      if((*ptrvar_i[in]<ptrval_i[in][0])
		 ||(*ptrvar_i[in]>ptrval_i[in][1]))
	      {
		printf("\nValue out of range (%d - %d) for argument %s\n\n",
		       ptrval_i[in][0],ptrval_i[in][1],olarg[in]);
		flagv=0;
	      }
	      if(type[in]==4) *ptrflag[in]=1;
	    }
	    else
	    {
	      *ptrflag[in]=1;
	      if(type[in]==3) *ptrvar_i[in]=1;
	    }
	  }
	  
/*            End of parsing loop             */
	}	

	if(flagv!=1)
	{
	  printf("Usage: %s",argv[0]);
	  for(in=0;in<Narg;in++)
	    if(type[in]%3==0) printf(" [%s]",olarg[in]);
	    else printf(" [%s %s]",olarg[in],olval[in]);
	  printf("\n");
	}

	if(flagv==2)
	  for(in=0;in<Narg;in++) printf("%s",olexp[in]);


/*            Freeing memory before terminating        */

	liberar_matriz_char(olarg,Narg0);
	liberar_matriz_char(olval,Narg0);
	free(ptrvar_f);
	liberar_matriz_float(ptrval_f,Narg0);
	free(ptrvar_i);
	liberar_matriz_int(ptrval_i,Narg0);
	free(ptrflag);
	free(type);

	
/*                Termination              */

	if(flagv!=1) exit(-1);

}

void genera_base( int verbose, int leff, char *base)
{
    char aux[20];

    switch(TYPE)
    { 
	case 0:
	    sprintf(base,"Log-Poisson.h%0.2f-coD%0.2f",
		    HINF,CODINF);
	    if(verbose)
		printf("Analyzing %d series of the type Log-Poisson; hinf= %0.2f, Dinf= %0.2f\n",
		       NSERIES,HINF,(double)D_space-CODINF);
	    break;
	case 1:
	    sprintf(base,"Log-Normal.mean%0.2f-sigma%0.2f",
		    MU,SIGMA);
	    if(verbose)
		printf("Analyzing %d series of the type Log-normal; mean= %0.2f, sigma= %0.2f\n",
		       NSERIES,MU,SIGMA);
	    break;
	case 2:
	    sprintf(base,"Log-Levi_mean%0.2f-sigma%0.2f-alpha%0.2f",
		    MU,SIGMA,ALPHA);
	    if(verbose)
		printf("Analyzing %d series of the type Log-Levi; mean= %0.2f, sigma= %0.2f, alpha= %0.2f\n",
		       NSERIES,MU,SIGMA,ALPHA);
	    break;
	case 3:
	    sprintf(base,"Binomial_h0%0.2f-h1%0.2f",HINF,H1);
	    if(verbose)
		printf("Analyzing %d series of the type binomial; h0= %0.2f, h1= %0.2f\n",
		       NSERIES,HINF,H1);
	    break;
	default:
	    sprintf(base,"bugged");
	    break;
    }
    if(leff>0)
    {
	sprintf(aux,"-size%d",leff);
	strcat(base,aux);
	if(verbose) printf("and size of %d points\n",leff);
	
    }

    if(D_space==1) strcat(base,"_1D");
    else strcat(base,"_2D");
}

Quality analiza_series( char *path)
{
    Quality output={-1.,-1.,-1.,-1.};

    char nombre[90],base[90];
    double *h_g,*h_w,*h_m,*h_wtmm;
    double *Dh_g,*Dh_w,*Dh_m,*Dh_wtmm;
    double *errDh_g,*errDh_w,*errDh_m,*errDh_wtmm;
    double sc;
    double shift;
 
    int graba;
    int dimy;
    int Nh_g,Nh_w,Nh_m,Nh_wtmm;
    int ir;


    dimy=(D_space==1)?1:LEFF;
    genera_base(VERBOSE,LEFF,base);

    if(FROM_DH)
    {
	graba=0;
	if(METHOD&0x01)
	{
	    sprintf(nombre,"Dh_GH_%s-N%d",base,NSERIES);
	    carga_Dh(nombre,&Nh_g,&h_g,&Dh_g,&errDh_g);
	}
 	if(METHOD&0x02)
	{
	    sprintf(nombre,"Dh_GMWP_%s-N%d",base,NSERIES);
	    carga_Dh(nombre,&Nh_w,&h_w,&Dh_w,&errDh_w);
	}
	if(METHOD&0x04)
	{
	    sprintf(nombre,"Dh_Mom_%s-N%d",base,NSERIES);
	    carga_Dh(nombre,&Nh_m,&h_m,&Dh_m,&errDh_m);
	}
	if(METHOD&0x08)
	{
	    sprintf(nombre,"Dh_WTMM_%s-N%d",base,NSERIES);
	    carga_Dh(nombre,&Nh_wtmm,&h_wtmm,&Dh_wtmm,&errDh_wtmm);
	}
   }
    else
    {
	graba=1;
	if(METHOD&0x01) estima_Dh_g(path,base,&Nh_g,&h_g,&Dh_g,&errDh_g);
	if(METHOD&0x02) estima_Dh_w(path,base,&Nh_w,&h_w,&Dh_w,&errDh_w);
	if(METHOD&0x04) estima_Dh_mom(path,base,&Nh_m,&h_m,&Dh_m,&errDh_m);
	if(METHOD&0x08)
	{
	    if(LASTWAVE) 
		estima_Dh_lw(path,base,&Nh_wtmm,&h_wtmm,&Dh_wtmm,&errDh_wtmm);
	    else estima_Dh_wtmm(path,base,&Nh_wtmm,&h_wtmm,&Dh_wtmm,&errDh_wtmm);
	}
	
/*       Correcting wavelet shift if required       */

	if(!(HOLDER)&&(Mod(METHOD,4)>=2))
	{
	    shift=shiftg-shiftw;
	    if(VERBOSE)
		printf("Correction exponent due to lack of translational invariance: %f\n",shift);
	    
	    for(ir=0;ir<Nh_w;ir++) h_w[ir]+=shift;
	}
    }
    

/*     Summarizing results        */

    if(METHOD&0x01)
    {
	if(VERBOSE)
	    printf("Quality estimation for derivative method\n========================================\n");
	sprintf(nombre,"Dh_GH_%s-N%d",base,NSERIES);
	output.g=registra_Dh(graba,Nh_g,nombre,1./((double)LEFF),h_g,Dh_g,
			     errDh_g);
    }


    if(METHOD&0x02)
    {
 	if(D_space==1) sc=S0*escala_wavelet_1D(LEFF);
	else sc=S0*escala_wavelet_2D(LEFF,LEFF);
	if(VERBOSE)
	    printf("Quality estimation for wavelet method\n=====================================\n");
	sprintf(nombre,"Dh_GMWP_%s-N%d",base,NSERIES);
	output.w=registra_Dh(graba,Nh_w,nombre,sc,h_w,Dh_w,errDh_w);
    }


    if(METHOD&0x04)
    {
	if(VERBOSE)
	    printf("Quality estimation for moment method\n=====================================\n");
	sprintf(nombre,"Dh_Mom_%s-N%d",base,NSERIES);
	output.m=registra_Dh(graba,Nh_m,nombre,sc,h_m,Dh_m,errDh_m);
    }
	
    if(METHOD&0x08)
    {
	if(VERBOSE)
	    printf("Quality estimation for WTMM method\n=====================================\n");
	sprintf(nombre,"Dh_WTMM_%s-N%d",base,NSERIES);
	output.l=registra_Dh(graba,Nh_wtmm,nombre,sc,h_wtmm,Dh_wtmm,
			     errDh_wtmm);
    }


/*	FREEING MEMORY BEFORE FINISHING */

    if(METHOD&0x01)
    {
	free(h_g);
	free(Dh_g);
	free(errDh_g);
    }
    if(METHOD&0x02)
    {
	free(h_w);
	free(Dh_w); 
	free(errDh_w);
    } 
    if(METHOD&0x04)
    {
	free(h_m);
	free(Dh_m); 
	free(errDh_m);
    }
    if(METHOD&0x08)
    {
	free(h_wtmm);
	free(Dh_wtmm); 
	free(errDh_wtmm);
    }
 
    return(output);

}

double media_exp( int dimx, double sc0, double **expon)
{
    double mean,norm,weight;
    int dimy;
    int ix,iy;

    dimy=(D_space==1)?1:dimx;

    mean=0.;
    norm=0.;
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	weight=pow(sc0,expon[iy][ix]);
	norm+=weight;
	mean+=expon[iy][ix]*weight;
    }
    }
    mean/=norm;

    return(mean);
}

void estima_Dh_g( char *path, char *base, int *Nh_g, double **h_g, 
		  double **Dh_g, double **errDh_g)
{
    char nombre[90];
    double **signal;
    double **expon_g;
    double *histo_g;

    double m_g[2],m_g0[2];
    double sc;
    
    int dimy;
    int in,ir;

    dimy=(D_space==1)?1:LEFF;

    signal=reservar_matriz(dimy,LEFF);
    expon_g=reservar_matriz(dimy,LEFF);
    histo_g=(double *) calloc(NBOX,sizeof(double));

/*     First passage is to compute the maxima      */

    m_g[0]=1e30;
    m_g[1]=-1e30;
    for(in=0;in<NSERIES;in++)
    {
	sprintf(nombre,"%s%s-N%05d",path,base,in);
	lee_datos(LEFF,nombre,signal);
	genera_expon(LEFF,signal,expon_g);
	if(D_space==1) extrema_lista(LEFF,expon_g[0],&m_g0[0]);
	else extrema(LEFF,LEFF,expon_g,&m_g0[0]);
	m_g[0]=fMin(m_g[0],m_g0[0]);
	m_g[1]=fMax(m_g[1],m_g0[1]);
    }


/*      Second passage: the histograms are calculated    */

     if(D_space==1) sc=1./(double)LEFF;
     else sc=escala_lineal_2D(LEFF,LEFF);

     for(in=0;in<NSERIES;in++)
     {
	sprintf(nombre,"%s%s-N%05d",path,base,in);
	lee_datos(LEFF,nombre,signal);
	genera_expon(LEFF,signal,expon_g);

/* 
   If this method is explicitly invoked, we take profit of that to
      save time in the estimation of the shift for the GMWP
*/
	acumula_histograma(LEFF,&m_g[0],expon_g,histo_g);
     }

     *Nh_g=calcula_Dh_histo(1./((double)LEFF),&m_g[0],histo_g,
				  h_g,Dh_g,errDh_g);

     shiftg=moda_por_histo(NBOX,m_g,histo_g);
	  

/*	FREEING MEMORY BEFORE FINISHING */

   
     free(histo_g);
     liberar_matriz(expon_g,dimy);
     liberar_matriz(signal,dimy);

}

void estima_Dh_w( char *path, char *base, int *Nh_w, double **h_w,
		  double **Dh_w, double **errDh_w)
{
    char nombre[90];
    double **signal;
    double **expon_w;
    double *histo_w;
    double *h_g,*Dh_g,*errDh_g;

    double m_w[2],m_w0[2];
    double sc,scg;

    int Nh_g;
    int dimy;
    int in,ir;

    dimy=(D_space==1)?1:LEFF;

    signal=reservar_matriz(dimy,LEFF);
    expon_w=reservar_matriz(dimy,LEFF);
    histo_w=(double *) calloc(NBOX,sizeof(double));
	
    if((METHOD&0x01)==0) // Implicit call to GH to obtain shiftg
    {
	estima_Dh_g(path,base,&Nh_g,&h_g,&Dh_g,&errDh_g);
	free(h_g);
	free(Dh_g);
	free(errDh_g);
    }

 
/*     First passage is to compute the maxima      */

    m_w[0]=1e30;
    m_w[1]=-1e30;
    for(in=0;in<NSERIES;in++)
    {
	sprintf(nombre,"%s%s-N%05d",path,base,in);
	lee_datos(LEFF,nombre,signal);

	m_w0[0]=1e30;
	m_w0[1]=-1e30;
	if(D_space==1)
	    calcula_multifractal_1D(LEFF,0,signal[0],expon_w[0],&m_w0[0]);
	else
	    calcula_multifractal_2D(LEFF,LEFF,0,signal,expon_w,&m_w0[0]); 

	m_w[0]=fMin(m_w[0],m_w0[0]);
	m_w[1]=fMax(m_w[1],m_w0[1]);

    }


/* Second passage: the histograms are calculated and the shift, corrected    */

     if(D_space==1) sc=escala_wavelet_1D(LEFF);
     else sc=escala_wavelet_2D(LEFF,LEFF);

     for(in=0;in<NSERIES;in++)
     {
	sprintf(nombre,"%s%s-N%05d",path,base,in);
	lee_datos(LEFF,nombre,signal);

	m_w0[0]=1e30;
	m_w0[1]=-1e30;
	if(D_space==1)
	    calcula_multifractal_1D(LEFF,VERBOSE,signal[0],expon_w[0],
				    &m_w0[0]);
	else calcula_multifractal_2D(LEFF,LEFF,VERBOSE,signal,expon_w,
				     &m_w0[0]);
	acumula_histograma(LEFF,&m_w[0],expon_w,histo_w);
     }

     *Nh_w=calcula_Dh_histo(sc,&m_w[0],histo_w,h_w,Dh_w,errDh_w);
     shiftw=moda_por_histo(NBOX,m_w,histo_w);

/*	FREEING MEMORY BEFORE FINISHING */

   
    free(histo_w);
    liberar_matriz(expon_w,dimy);
    liberar_matriz(signal,dimy);

}

void estima_Dh_mom( char *path, char *base, int *Nh_m, double **h_m,
		    double **Dh_m, double **errDh_m)
{
    char nombre[90];

    double **signal;
    double **moments;
    double *taup;
    double sc;

    int dimy;
    int in,ip,it;

    dimy=(D_space==1)?1:LEFF;

    signal=reservar_matriz(dimy,LEFF);
    moments=reservar_matriz(Nmom,Ndist);
    taup=(double *) calloc(Nmom,sizeof(double));
 	
 
/*     First step is to compute the moments      */

    for(in=0;in<NSERIES;in++)
    {
	sprintf(nombre,"%s%s-N%05d",path,base,in);
	lee_datos(LEFF,nombre,signal);
	acumula_momentos(signal,moments);
    }

/*        It is necessary to normalize by the number of series    */

    for(ip=0;ip<Nmom;ip++)
    {
    for(it=0;it<Ndist;it++)
    {
	moments[ip][it]/=(double)NSERIES;
    }
    }

/*      Second step: taup is evaluated    */

    calcula_taup(moments,taup);

/*      Thirs step: from taup, Dh is evaluated    */


    *Nh_m=simple_Legendre_transform(taup,h_m,Dh_m,errDh_m);

/*	FREEING MEMORY BEFORE FINISHING */

    free(taup);
    liberar_matriz(moments,Nmom);
    liberar_matriz(signal,dimy);

}

void estima_Dh_wtmm( char *path, char *base, int *Nh_wtmm, double **h_wtmm,
		     double **Dh_wtmm, double **errDh_wtmm)
{
    char nombre[90];

    double **ExtWTlis, **Z;
    double *scArray;
    double **sTq, **sTqLogT;
    double **signal;
    double **moments;
    double *taup;
    double sc;
    float MinScale,MaxScale;

    int *n_ext;
    int nsc=0; 
    int dimy;
     int ind_q, ind_sc;
    int in,ix,iy,ip;

/*      Initialization of the dimensions    */

    dimy=(D_space==1)?1:LEFF;
    MinScale=MINSCALE;
    MaxScale = ScRatio*(((float)LEFF));
    for (sc=MinScale; sc<=MaxScale; sc*=ScStep) nsc++;

/*     Memory initialization     */

    signal=reservar_matriz(dimy,LEFF+1);
    ExtWTlis=reservar_matriz(nsc,LEFF+1);
    Z=reservar_matriz(nsc,Nmom+1);
    moments=reservar_matriz(Nmom,Ndist);
    scArray=(double*)calloc(nsc, sizeof(double));
    n_ext = (int*)calloc(nsc+1,sizeof(int));
    taup=(double *) calloc(Nmom,sizeof(double));
    *h_wtmm=(double *) calloc(Nmom,sizeof(double));
    *Dh_wtmm=(double *) calloc(Nmom,sizeof(double));
    *errDh_wtmm=(double *) calloc(Nmom,sizeof(double));
    *Nh_wtmm=Nmom;
 	
    if(CANONICAL)
    {
	sTq=reservar_matriz(nsc,Nmom);
	sTqLogT=reservar_matriz(nsc,Nmom);
    }


/*     Value initialization    */

    for (ind_sc=1, scArray[0]=MinScale; ind_sc<nsc; ind_sc++) 
	scArray[ind_sc] = scArray[ind_sc-1] * ScStep;
 
/*    Processing     */
    
    if(VERBOSE) 
    {
	printf("Compute Wavelet Transform, Wavelet Extrema and Partition Function\n" "for all signals and over all scales.\n");    
	if(NSERIES>1) 
	    printf("The statistics of the different signals are accumulated.\n");
    }

    for(in=0;in<NSERIES;in++)
    {
	sprintf(nombre,"%s%s-N%05d",path,base,in);
	lee_datos(LEFF,nombre,signal);

   /* Generate and track maxima lines for extrema extraction */
	for(iy=0;iy<dimy;iy++)
	{
/* clean here the WTMM first : ExtWTlis = 0 everywhere, even if we
   rewrite on it  */
	    limpia(LEFF+1,nsc,ExtWTlis);
	    WTMM_PF_compute(signal[iy],LEFF,(double)wav_wtmm,
			    ord_der_wtmm, Nmom,moms,nsc,scArray,
			    ExtWTlis,Z,n_ext);

	    if(CANONICAL) 
		canonMeanCompute(ExtWTlis,LEFF,n_ext,Nmom,moms,nsc,scArray,
				 sTq, sTqLogT);
	}
    }


  /** Finally compute the exponents tau(q) and the spectrum (h,D(h)) from the
   ** values of the WT over the detected extrema **/
    if(DIRECT ) 
    {
	if(VERBOSE) 
	    printf("Approximate Legendre Transform for estimating multifractal exponents\n"
		   "with DIRECT method\n");
	directSpecCompute(Z,LEFF,Nmom,moms,nsc,scArray,taup,*h_wtmm,*Dh_wtmm);
    
    } 
    else 
    {
	if(VERBOSE)
	    printf("Approximate Legendre Transform for estimating multifractal exponents\n"
		   "with CANONICAL method\n");
	canonSpecCompute(sTq,sTqLogT,LEFF,Nmom,moms,nsc,scArray,
			 taup,*h_wtmm,*Dh_wtmm);
    }

/*       Padding proxy errorbars and shifting spectrum       */

    for(ip=0;ip<Nmom;ip++)
	errDh_wtmm[0][ip]=0.2; // conventionally


/*	FREEING MEMORY BEFORE FINISHING */

    free(scArray);
    free(n_ext);
    free(taup);
    if(CANONICAL)
    {
	liberar_matriz(sTqLogT,nsc);
	liberar_matriz(sTq,nsc);
    }
    liberar_matriz(ExtWTlis,nsc);
    liberar_matriz(Z,nsc);
    liberar_matriz(moments,Nmom);
    liberar_matriz(signal,dimy);

}

void estima_Dh_lw( char *path, char *base, int *Nh_l, double **h_l,
		   double **Dh_l, double **errDh_l)
{
    FILE *chan;
    char nombre[90];
    float aux1,aux2;
    int ip,in;


/*   First step: writing the script which will be linked to LastWave */
 
    chan=fopen("LW.scr","wt");

/*     Lines creating the objects to be employed     */

    fprintf(chan,"pf_lp=[new &PF]\npf_tot=[new &PF]\npf_tot2=[new &PF]\nw_lp=[new &wtrans]\nw=[new &wtrans]\n");


/*    Global parameters defining the wavelet transform   */

    fprintf(chan,"amin = %d\nomax = %d\n",LW_AMIN,LW_OMAX);

/*    Defining the list of moments to be used           */

    fprintf(chan,"qlist = {");
    for(ip=0;ip<Nmom;ip++) fprintf(chan,"%0.2f ",moms[ip]);
    fprintf(chan,"}\n");

/*      Reading the first series and computing the WT  */

    fprintf(chan,"read 0w '%s%s-N%05d.txt'\n",path,base,0);
    fprintf(chan,"cwtd w amin omax 10 -e 0\nextrema w w.extrep\n");

/*  The first time we generate over pf_tot directly, then we sum up  */

    fprintf(chan,"pf wtmm pf_tot w.extrep qlist\n");

    for(in=1;in<NSERIES;in++)
    {
	fprintf(chan,"w=[new &wtrans]\n");
	fprintf(chan,"read 0w '%s%s-N%05d.txt'\n",path,base,in);
	fprintf(chan,"cwtd w amin omax 10 -e 0\nextrema w w.extrep\n");
	fprintf(chan,"pf wtmm pf_lp w.extrep qlist\npf_tot=[pf add pf_lp pf_tot]\n");
    }

/* 
   With the accumulated partition function, we obtain the singularity spectrum.
*/

    fprintf(chan,"log2aMin=0\nlog2aMax=3\n");
    fprintf(chan,"Dh = XY(Zero(pf_tot.qnumber),Zero)\n");
    fprintf(chan,"i = 0\n");
    fprintf(chan,"foreach q pf_tot.qlist {\n");
    fprintf(chan,"Dh.X[i] = [stats fit [pf get h pf_tot q] -x log2aMin log2aMax][0]\n");
    fprintf(chan,"Dh.Y[i] = [stats fit [pf get d pf_tot q] -x log2aMin log2aMax][0]\n");
    fprintf(chan,"i+=1\n}\n");

/*      Sorting the results     */

    fprintf(chan,"sort Dh\n");

/*   Writing the results to a dummy file from which we will re-obtain them  */

    fprintf(chan,"write Dh 'Dh_lw_aux' -h\n");


/*   Forcing exit the script   */

    fprintf(chan,"return 0\nexit\n");

    fclose(chan); // closing the automatically generated script

/*      Creating the executing script    */

    chan=fopen("LW.exe","wt");
    fprintf(chan,"#!/bin/sh\n");
    fprintf(chan,"%s < LW.scr\n",LWCALL);
    fprintf(chan,"rm LW.scr LW.exe\n");
    fclose(chan);
    

/*      Now executing the script      */

    system("chmod a+x LW.exe; LW.exe");


/*   Reading the auxiliary spectrum and transferring it to the data  */

    *Nh_l=Nmom;
    *h_l=(double *) calloc(*Nh_l,sizeof(double));
    *Dh_l=(double *) calloc(*Nh_l,sizeof(double));
    *errDh_l=(double *) calloc(*Nh_l,sizeof(double));

    chan=fopen("Dh_lw_aux","rt");
    for(ip=0;ip<*Nh_l;ip++)
    {
	fscanf(chan,"%lf  %lf",&(h_l[0][ip]),&(Dh_l[0][ip]));
	h_l[0][ip]-=1.; 
              // We shift the singularity range to fit that of derivatives
	errDh_l[0][ip]=0.2; // conventional
    }
    fclose(chan);

/*   Removing the auxiliary file    */

    system("rm Dh_lw_aux");

}

void WTMM_PF_compute( double *signal, int dimx, double wav, int order,
		     int nq, double *qArray, int nsc, double *scArray, 
		     double **ExtWTlis, double **Z, int *n_ext ) {
  /* =====================================================================
   *
   * Find the extrema of WT for each scale of analysis, by performing the  
   * following steps:
   * 1) Wavelet convolution of the signal for increasing wavelet scale.
   * 2) Locate the local maxima of the absolute value of wavelet
   *    coefficient as a function of time for each wavelet scale.
   * 3) Check whether a local maximum at a given wavelet scale is located
   *    close to a maximum at a smaller scale - if yes connect both maxima,
   *    otherwise cancel it. Generate maxima lines.
   * 4) Check that the number of maxima at larger scales is less or equal
   *    to that at a smaller scale. 
   * 5) Track maxima lines for increasing wavelet scale by choosing at each
   *    scale the supremum between all previous values at smaller scales.
   *
   * Parameters:
   *    - signal : original signal,
   *    - order : variable used for the choice of the wavelet,
   *    - maxfilt : maximal size of the wavelet filter (pre-computed),
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nsc, scArray : resp. no of sceles and list of scales.
   * Outputs:
   *    - Z : 2d tabular (scale by moment) storing the factors of the 
   *      partition function,
   *    - ExtWTlis : 2d tabular (scale by list of extrema) storing the WTMM,
   *      i.e. the values of the WT over selected extrema,
   *    - n_ext : tabular storing the number of extrema at each scale,
   * Returns the number icolor of data in color representation, if computed.
   *
   * ===================================================================== */
  
  int i, j, jj;
  double *maxsig, *wavesig, *pt, sc;
  double *wtrans;
  int *wtrans_ind, *maxsig_ind, *ptl;
  int ind_sc, ind_q, n_max;

   /* Allocations of memory
    * Note: making the allocations here (once for all, instead of making
    * it for each scale) enables to reduce time expense and allows
    * comparisons between scales */
  if( /* Arrays of the wavelet transform */
     (wtrans = (double*)calloc(dimx+1,sizeof(double))) == NULL ||
     (wtrans_ind = (int*)calloc(dimx+1,sizeof(int))) == NULL || 
     /* Arrays to manipulate the maxima through the scales */
     (maxsig = (double*)calloc(dimx+1,sizeof(double))) == NULL ||
     (maxsig_ind = (int*)calloc(dimx+1,sizeof(int))) == NULL ){
    fprintf(stderr,"\n Error allocation in WTMM_PF_compute");
    exit(-1);
  }
  
  jj=0; /* used to test wether the number of maxima at larger scales is less 
	 * or equal to that at a smaller scale.*/
  
  for( ind_sc=0; ind_sc<nsc; ind_sc++ ) 
  {
      sc=scArray[ind_sc]; /* current scale of analysis */
      
      /* Wavelet transform of the signal for increasing wavelet scale. */
      WT_transform( signal, dimx, wav, order, sc,  wtrans );
     
    /** Find the local maxima of the wavelet coefficient for each scale **/
    n_max = WTMM_find_extrema( wtrans, wtrans_ind, dimx, sc); 
    /* n_max : number of extrema at current scale */
    
    
    /** Tracking the maxima lines: test for supremum 
     ** Generate maxima lines and compute the partition function. */
    if( ind_sc > 0)    /* i.e.: if sc > min_sc */
      n_ext[ind_sc] =
	PF_extrema_track( wtrans, wtrans_ind, maxsig, maxsig_ind, 
			  jj, n_max, nq, qArray, nsc, scArray, 
			  ExtWTlis[ind_sc], Z[ind_sc] );
    
    /* Update for the next scale */
    jj = n_max;

    pt = maxsig;
    maxsig = wtrans; /* store in maxsig the WT of the current scale in order
		      * to compare it with the next scale */
    wtrans = pt; /* change the variable pointed by wtrans, so that
		  * further modification of wtrans won't modify maxsig */
    ptl = maxsig_ind;
    maxsig_ind = wtrans_ind; /* ibid with the indexes */ 
    wtrans_ind = ptl;  /* ibid with the indexes */
    
  } /* end loop over the scales: "for (ind_sc=0;..." */

  /* Free memories */
  if(wtrans) free(wtrans);
  if(wtrans_ind) free(wtrans_ind);
  if(maxsig) free(maxsig);
  if(maxsig_ind) free(maxsig_ind);
  
}


/* ===================================================================== */
int WTMM_find_extrema( double *wtrans, int *wtrans_ind, int dimx, 
			double sc ) {
  /* =====================================================================
   *
   * Find the local maxima of the wavelet coefficient for each scale.
   *
   * Parameters:
   *    - sc : current scale of analysis,
   *    - wtrans : 1d array storing  the values of the WT at scale sc,
   *    - dimx : common size of both  wtrans and wtrans_ind.
   * Outputs:
   *    - wtrans : it is modified to finally store the values of the extrema
   *      only,
   *    - wtrans_ind : 1d array storing the indexes of the selected 
   *      extrema.
   * Returns the number n_max of selected extrema at scale sc.
   *
   * ===================================================================== 
   * Called by : WTMM_PF_compute  
   * ===================================================================== */
  
  int n_max=0;
  double temp, temp1;
  int i,j,sign;
  
  /* useless: for (i=0; i<dimx+1; i++) wtrans_ind[i] = 0; */
  
  sign = SIGN(wtrans[1] - wtrans[0]);
  temp1 = 0.; 
  for( j=0, i=2; i<dimx; i++ ) {  
    if ((fabs(wtrans[i] - wtrans[i-1]) > 0.0) && 
	((sign == 1) && ((SIGN(wtrans[i] - wtrans[i-1])) == -1))) {
      temp = wtrans[i-1];
            
      n_max++;
      /*  maximum condition */
      wtrans[j]=temp;
      wtrans_ind[j]=i-1;
      temp1=temp;
      j++;
    }
    sign=SIGN(wtrans[i]-wtrans[i-1]);
  } /* end loop on the filter "for( j=0,..." */
  
  return n_max;
}



/* ===================================================================== */
int PF_extrema_track( double *wtrans, int *wtrans_ind, 
		      double *maxsig, int *maxsig_ind,
		      int jj, int  n_max,
		      int nq, double *qArray, int nsc, double *scArray, 
		      double *ExtWTlis, double *Z ) {
    /* ===================================================================== 
   *
   * Check whether the local maximum is located close to a maximum at a 
   * smaller scale. if yes connect both maxima, otherwise cancel it. 
   *
   * Parameters:
   *    - wtrans, wtrans_ind : 1d arrays storing resp. the values and the
   *      indexes of the WT at the current scale,
   *    - maxsig, maxsig_ind : 1d arrays storing resp. the values and the
   *      indexes of the WT at the previous scale,
   *    - n_max : number of local extrema,
   *    - jj : index of the last detected extrema in the list of extremas
   *      already selected.
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nsc, scArray : resp. no of sceles and list of scales,
   * Outputs:
   *    - Z : 1d tabular (indexed by moment) storing the factors of the 
   *      partition function at current scale,
   *    - ExtWTlis : 1d array storing the WTMM over selected extrema
   *      at current scale,
   * Returns the final number of tracked global extrema at current scale. 
   *
   * ===================================================================== 
   * Called by : WTMM_PF_compute  
   * ===================================================================== */
  
  int i1=0, i2=0;
  int n_ext=0,ind_q;
  
  /* With jj, check that the number of maxima at larger scales is less or 
   * equal to that at a smaller scale. */
  while (((i1-1) < n_max) && ((i2-1) < jj)) 
  {

      /* Choose the supremum between all previous values at smaller scales. */
    /* Activated only if the SUP_WTMM option is given                     */
      
      if(SUP_WTMM)
      {
	  if ((wtrans_ind[i1] - maxsig_ind[i2]) <=
	      (maxsig_ind[i2+1] - wtrans_ind[i1]))
	      wtrans[i1] = fMax( wtrans[i1], maxsig[i2] );
	  else
	      wtrans[i1] = fMax( wtrans[i1], maxsig[i2+1] );
      }

    /* Store in ExtWTlis the list of the values of the WT over the detected 
     * maxima, i.e. the final WTMM */
    
      ExtWTlis[n_ext++] = wtrans[i1];
    
    /* ...and also compute directly the factors of the partition funtion
     * from the WTMM */
      for ( ind_q=0; ind_q<nq; ind_q++ )
	  Z[ind_q] += pow( wtrans[i1], qArray[ind_q]);
    /* Note that in the case of several signals (NSERIES>1), we simply
     * add the factors of the partition function */

      i1++;
      i2++;
      while ((i2 < jj) && (wtrans_ind[i1] >= maxsig_ind[i2]))
      i2++;
      i2--;
  }	  
  
  return n_ext;
}


/* ===================================================================== */
void WT_transform( double *signal, int dimx, double wav, int order,
		   double sc, double *wtrans ){
  /* ===================================================================== 
   *
   * Compute the wavelet transform of the signal at scale sc
   *
   * Parameters:
   *    - signal : original signal,
   *    - dimx : size of the signal,
   *    - order : variable used for the choice of the wavelet,
   *    - sc : current scale of analysis,
   *    - maxfilt : maximal size of the wavelet filter.
   * Output:
   *    - wtrans : 1d tabular storing the values of the WT at scale sc,
   *
   * ===================================================================== 
   * Note : signal is zero padded to remove aliasing. FFT is used
   * Called by : WTMM_PF_compute  
   * ===================================================================== */
 
  int i, j;
  int tempi;
  double *wave,*auxsig,*auxwtrans;
  

 /* allocate memory for the wavelet and extended signal and wavele transform */
  if( 
      ((wave = (double*)calloc(2*dimx,sizeof(double))) == NULL )|| 
      ((auxsig = (double*)calloc(2*dimx,sizeof(double))) == NULL )|| 
      ((auxwtrans = (double*)calloc(2*dimx,sizeof(double))) == NULL ) 
      )
  {
    fprintf(stderr,"\n Error allocation in WT_transform");
    exit(-1);
  }
  
  asigna_lista(dimx,signal,auxsig); // auxsig is zero-padded

  /* Create the wavelet filter */
  define_wavelet( 2*dimx, sc, wav, order, wave ); 
  
   /* Obtaining the wavelet transforms   */
  W_projection( auxsig, wave, 2*dimx, sc, auxwtrans );
   
   /*  We just keep the first half */

  asigna_lista(dimx,auxwtrans,wtrans);
      
  if(wave)  free(wave);
  if(auxsig)  free(auxsig);
  if(auxwtrans)  free(auxwtrans);
}


/* ===================================================================== */
void W_projection( double *signal,  double *filt, int dimx, double sc, 
		   double *wtrans ) 
{
  /* ===================================================================== 
   *
   * Naive convolution + normalization to compute the wavelet coefficient
   *
   * Parameters:
   *    - signal : original signal,
   *    - dimx : size of the signal, the filter and the wavelet transform
   *    - sc : current scale of analysis,
   *    - filt : wavelet filter used to convolve the signal at scale sc,
   * Output:
   *    - wtrans : values of the WT at scale sc,
   *
   * ===================================================================== 
   * Called by : WT_transform
   * ===================================================================== */
  
  int i;
  
  convuelve_1D(dimx,signal,filt,wtrans);
  /* absolute value of wavelet coefficient */
  for(i=0;i<dimx;i++) wtrans[i] = fabs(wtrans[i]);  
}


/* ===================================================================== */
void define_wavelet( int dimx, double sc, double wav, int order, 
			double *wavesig ) {
  /* =====================================================================
   *
   * Function computing the wavelet function, which can be a lorentzian
   * wavelet (wav>0) or a gaussian wavelet (wav<0) - version II
   *
   * Parameters :
   *    - n : dimension of the output wavelet function (this will be
   *      typically TIMESSC*sc),
   *    - sc : scale of analysis, 
   *    - wav : order of the lorentzian wavelet when >0, otherwise 
   *      it is a gaussian wavelet,
   *    - order : order of derivative of the gaussian.
   * Returns the wavelet function in wavesig.
   *
   * ===================================================================== 
   * Called by : WT_transform
   ===================================================================== */

  int ix;
  double t;
  
  for( ix=0; ix<dimx; ix++ )
  {
    t = (double)ix;
    if(ix>=dimx/2) t -= (double)dimx;
    t = t/sc;
    if(wav<0.) wavesig[ix] = pow(1.+t*t,-1); /* lorentzian */
    else if (wav == 0.) wavesig[ix] = exp(-t*t/2.) * cos(5.*t); /* morlet */
    else /*if(wav>0.)*/  wavesig[ix] = comp_coefficient(order,t); /* gaussian */
    wavesig[ix]/=sc;
  }

}



/* ===================================================================== */
double comp_coefficient(int order, double t) {
  /* =====================================================================
   *
   * Computes continuous Gaussian wavelet functions (0 to 5th derivative).
   *
   * ===================================================================== 
   * Called by : WT_define_wavelet 
   ===================================================================== */
  
    switch (order) {
	case 0:  return exp(-0.5*t*t);
	case 1:  return -t * exp(-0.5*t*t);
	case 3:  return t * exp(-0.5*t*t) * (3-t*t);
	case 4:  return exp(-0.5*t*t) * (pow(t,4.0) - 6*t*t + 3);
	case 5:  return -t * exp(-0.5*t*t)* (pow(t,4.0) - 10*t*t + 15);
	default:
	case 2:  return (1-t*t) * exp(-0.5*t*t);
    }
}


/**
 ** ===================================================================== **
 ** METHOD I : direct computation of the density spectrum through the 
 **            exponent tau(q) and the partition function
 ** ===================================================================== **
 **/

/* ===================================================================== */
void directSpecCompute( double **Z, int dimx, 
			int nq, double *qArray, int nws, double *wsArray,
			double *Tauq, double *H, double *Dh ) {
  /* ===================================================================== 
   *
   * Direct method to compute the spectrum thanks to the Legendre transform
   *
   * Parameters:
   *    - Z : 2d tabular (scale by moment) storing the factors of the
   *      partition function,
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nws, wsArray : resp. no of scales and list of scales,
   * Outputs:
   *    - Tauq : multifractal exponents,
   *    - H : singularity exponents,
   *    - Dh : spectrum.
   *
   * =====================================================================
   * The Legendre transform is approximated through the relations:
   *            h = \frac{d\tau}{dq} 
   *            D(h) = q h - tau(q) 
   * ===================================================================== 
   * Note : not implemented in the LastWave toolbox.
   * ===================================================================== */
  
  double *sumpf, *sumscpf;
  int i, count=0, ind=0;
  double sumsc=0., sumsqsc=0.;
  int ind_q, ind_ws;
  double min_q, max_q;
  double Lws, LZ;
  
  int ind_ext; /* index of extrema */
  double wtext; /* local variable to store the value of the wavelet transform
		 * over the extrema */
  
  /* minimum considered moment */
  min_q=qArray[0];
   
  /* Local allocations */
  if( (sumpf=(double*)calloc(nq,sizeof(double))) == NULL ||
      (sumscpf=(double*)calloc(nq,sizeof(double))) == NULL ) {
    fprintf(stderr,"\n Error allocation in function directSpecCompute");
    exit(-1);
  }

  /* Initialization of both local tabulars */
  for(ind_q=0; ind_q<nq; ind_q++) sumpf[ind_q] = sumscpf[ind_q] = 0.;
  
  for (ind_ws=1; ind_ws<nws; ind_ws++) {
    /* log of the current scale */
    Lws = log(wsArray[ind_ws])/log(10.);

      /* number of scales considered till now */
    count ++;
    /* sum of scales */
    sumsc += Lws;
    /* sum of squared scales */
    sumsqsc += Lws*Lws;
    for( ind_q=0; ind_q<nq; ind_q++ ) {
	/* log of the partition function */
	LZ = log( Z[ind_ws][ind_q])/log(10.);
	/* sum of the partition function */
	sumpf[ind_q] += LZ;
	/* sum of partition function weightened by the scale  */
	sumscpf[ind_q] += Lws * LZ;	
    }
  } /* end loop over the scales: "for (ind_ws=0;..." */
  
   /** Compute the tauq */
  for(ind_q=0; ind_q<nq; ind_q++) {
    /** Compute the tauq exponent corresponding to qArray[ind_q]
     * */
    Tauq[ind_q] = (sumscpf[ind_q]*count - sumsc*sumpf[ind_q]) 
      / (count*sumsqsc - sumsc*sumsc);
  } /* end of the 1st loop over the moments: "for (ind_q=0;..." */
    

 /* Compute the spectrum (h,D(h)) */
  for(ind_q=1; ind_q<nq-1; ind_q++) {

    /** First compute the singularity exponent h:
     *        h = \frac{d\tau}{dq}  */
    H[ind_q] = (Tauq[ind_q+1]-Tauq[ind_q-1])/(qArray[ind_q+1]-qArray[ind_q-1]);
    /* (tau[iq+1]-tau[iq-1])/(qArray[iq+1]-qArray[iq-1]) is an approximation
     * of the derivative \frac{d\tau}{dq} in iq */

    /** Then compute the density spectrum:
     *        D(h) = q h -tau(q) */
    Dh[ind_q] = qArray[ind_q]*H[ind_q] - Tauq[ind_q]; 
    /* variant for numerical percision : */
    /* Dh[ind_q] = qArray[ind_q] / (qArray[ind_q+1]-qArray[ind_q-1])
     * (Tauq[ind_q+1]-Tauq[ind_q-1]) - Tauq[ind_q]; */
    
  } /* end of the 2nd loop over the moments: "for (ind_q=0;..." */
  
}



/**
 ** ===================================================================== **
 ** METHOD II : canonical computation through the canonical formula for
 **             the partition function and averaged exponents 
 ** ===================================================================== **
 **/


/* ===================================================================== */
void canonMeanCompute( double **ExtWTlis, int dimx, int *n_ext,
		       int nq, double *qArray, int nws, double *wsArray,
		       double **sTq, double **sTqLogT 
		       /* double **logSTq */ ) {
  /* ===================================================================== 
   *
   * Function used in the canonical method proposed to compute the spectrum 
   * through to the Legendre transform - part I
   * Averages of multifractal are computed for the different scales.
   *
   * Parameters:
   *    - ExtWTlis : 2d tabular (scale by moment) storing the WTMM, i.e.  
   *      the values of the WT computed over selected extrema,
   *    - n_ext : number of extrema at each scale, 
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nws, wsArray : resp. no of sceles and list of scales.
   * Outputs:
   *    - sTq, sTqLogT : 2d (scale by moment) intermedary tabulars used
   *      further to compute the different multifractal exponents and
   *      parameterized similarly to the factors of the partition function.
   *      They look like:
   *         sTq[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q
   *         sTqLogT[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q * log(|WT(x,ws)|)
   *         logSTq[ws,q] =  log(\sum{x \in L(ws)} |WT(x,ws)|^q)
   *
   * =====================================================================     
   *
   * The Legendre transform is approximated thanks to the method inspired
   * by Chhabra and Jensen: 
   *       \tilde{T_\psi} [s](q,x,a)=
   *                         \frac{T_\psi[s](x,a)}{Z(q,a)}
   * then the averages below are computed:
   *       < h >(q,a) = \sum_{(x,a)} 
   *                    \tilde{T_\psi} [s](q,x,a) \ln |T_\psi[s](x,a)|
   *       Dh(q,a) = \sum_{(x,a)} 
   *                 \tilde{T_\psi} [s](q,x,a) \ln \tilde{T_\psi}[s](x,a)
   * and, finally, the slopes of these quantities will provide the 
   * multifractal exponents (see canonSpecCompute below).
   *
   * ===================================================================== 
   * Note: in LastWave toolbox, this is the goal of the function
   *                PFComputeOneScaleF (one scale only)
   * in file pf_lib.c (package package_wtmm) called by: 
   *                ComputePartFuncOnExtrep (all scales)
   * in file pf_functions (package package_wtmm).
   * ===================================================================== */
  
  double *tempTq,*tempLogT;
  double *tempWT; 
  int ind_ws, ind_q, i, imin;
  int maxsize=-1, size;
  double tm, q;
  
  /* Find the maximum number of extrema over scales for further useful
   * allocations */
  for( i=0; i<nws; i++ )   maxsize = Max( maxsize, n_ext[i] );
  
  /* Local allocation of space for tempTq and tempLogT */
  if( (tempTq = (double *) malloc(2*maxsize*sizeof(double))) == NULL) {
    fprintf(stderr,"\n Error allocation in canonMeanCompute");
    exit(-1) ;
  }
  tempLogT = tempTq + maxsize;
  
  /* DEBUG   for( ind_ws=0; ind_ws<nws; ind_ws++ )... */
  for( ind_ws=1; ind_ws<nws; ind_ws++ ) {
    
    /* pointer on the list of selected extrema at scale ws */
    tempWT = ExtWTlis[ind_ws];
    /* number of extrema at this scale */
    size = n_ext[ind_ws];

    /* Rearrange in incresing order the values of the WTMM over the extrema */
    qsort((void*)tempWT,size,sizeof(double),
	  (int (*)(const void*,const void*))compare);
    /* check that it is the same qsort for your compilator: it may depend
     * on the libraries used... */

    /* Find the first occurence of non null WTMM */
    imin = 0;
    while(imin<size && tempWT[imin] == 0.)  imin++;
    
    /* Approximation of the Legendre transform */
    for(ind_q=0; ind_q<nq; ind_q++) {
      
      /* First determine the current q */
      q = qArray[ind_q];
      
      /* Do we want T/tm to be >= 1 or <= 1 ? */
      if(q >= 0.)	tm = tempWT[imin];
      else	tm = tempWT[size-1];
      /* DEBUG tm=1.; */
      
      /* We compute Tq and Log(T/tm) */
      for( i=imin; i<size; i++ )      {
	tempTq[i] = pow( tempWT[i], q );
	tempLogT[i] =log(tempWT[i]/tm)/log(10.);
      }
      
      /** Note that the sum below allow to compute statistical variables 
       ** over several signals (NSERIES>1) by simply adding their values **/

      /* We compute sTq and sTqLogT */
      if(q >= 0.)      
	for( i=imin; i<size; i++ )	{
	  sTq[ind_ws][ind_q] += tempTq[i]; /* pow( tempWT[i], q ); */
	  sTqLogT[ind_ws][ind_q] += tempTq[i] * tempLogT[i];
	  /* pow( tempWT[i], q ) *log( tempWT[i]/tm )/log(10.); */
	}
      
      else 
	for( i=size-1; i>=imin; i-- )	{
	  sTq[ind_ws][ind_q] += tempTq[i];/* pow( tempWT[i], q ); */
	  sTqLogT[ind_ws][ind_q] += tempTq[i] * tempLogT[i]; 
	  /* pow( tempWT[i], q ) * log( tempWT[i]/tm )/log(10.); */
	}
      
      /* sTqLogT = sTqLogT + log(tm)*sTq */
      sTqLogT[ind_ws][ind_q] += log(tm)/log(10.)*sTq[ind_ws][ind_q];
      
      /* We compute LogSTq
	 if(sTq[ind_ws][ind_q] != 0.0)	{
	 logSTq[ind_ws][ind_q] = 
	 log(sTq[ind_ws][ind_q]/((double) (size-imin)))/log(10.);
	 }    else 	{
	 logSTq[ind_ws][ind_q] = 0.;
	 } 
      */ 
      
    } /* end loop over the moments: "for (ind_q=0;..." */
    
    
  } /* end loop over the scales: "for (ind_ws=0;..." */

  
  if(tempTq) free(tempTq); /* then tempLogT is automatically free */
  
}


/* ===================================================================== */
void canonSpecCompute( double **sTq, double **sTqLogT, /* double **logSTq */
		       int dimx, 
		       int nq, double *qArray, int nws, double *wsArray,
		       double *Tauq, double *H, double *Dh ) {
  /* ===================================================================== 
   *
   * Function used in the canonical method proposed to compute the spectrum 
   * through to the Legendre transform - part II
   * Estimation of multifractal exponents are realized through linear
   * regression over the scales.
   *
   * Parameters:
   *    - sTq, sTqLogT : 2d (scale by moment) tabulars parameterized
   *      similarly to the factors of the partition function,
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nws, wsArray : resp. no of scales and list of scales,
   * Outputs:
   *    - Tauq : multifractal exponents,
   *    - H : singularity exponents,
   *    - Dh : spectrum of singularity.
   *
   * ===================================================================== 
   * Note: in LastWave toolbox, this operation is realized through the 
   * functions 
   *        tauqSpectrum and singSpectrum
   * of the script file wtmm1d.pkg (package scripts). 
   * These scripts make an implicit use of the functions:
   *        PFAccessTQFloat    PFAccessHQFloat    PFAccessDQFloat
   * of file pf_lib.c (package package_wtmm), and of the function :
   *        LineFitSig
   * of file signal_function.c (package package_signal), rewritten below.
   * ===================================================================== */
  
  int ind_ws, ind_q;
  double q, stq;
  double *tq, *hq, *dq;
  double *LwsArray;
  int ss;
  
  /* Local allocations */
  if( (tq=(double*)calloc(3*nws,sizeof(double))) == NULL ||
      (LwsArray=(double*)calloc(nws,sizeof(double))) == NULL ) {
    fprintf(stderr,"\n Error allocation in function canonSpecCompute");
    exit(-1);
  } 
  hq = tq + nws;
  dq = tq + 2*nws;
  
  /* Array storing the log of the scales */
  for( ind_ws=0; ind_ws<nws; ind_ws++ )
    LwsArray[ind_ws] = log(wsArray[ind_ws])/log(10.);
  
  for(ind_q=0; ind_q<nq; ind_q++) { 
    
    /* consider the current moment for computing the variables */
    q=qArray[ind_q];
    
    for( ind_ws=0; ind_ws<nws; ind_ws++ ) { 
      stq = sTq[ind_ws][ind_q];
      
      /** Compute the canonical multifractal exponents tq */ 
      if(stq == 0.)	tq[ind_ws] = 0.;
      else   tq[ind_ws] = (double) log(stq)/log(10.);
      /* equivalent to the function PFAccessTQFloat */
      
      /** Compute the canonical singularity exponents hq:
       *      <h>(q,ws) = \sum_{(x,ws)} 
       *         \tilde{WT}(q,x,ws) \ln |WT(x,ws)|
       * where:
       *   \tilde{WT}(q,x,ws)= WT(x,ws) / Z(q,ws)
       * by using: 
       *     sTq[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q
       *     sTqLogT[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q * log(|WT(x,ws)|)
       */
      if(stq == 0.)	hq[ind_ws] = 0.;
      else   hq[ind_ws] = (double)(sTqLogT[ind_ws][ind_q] / stq);
      /* equivalent to the function PFAccessHQFloat in file pf_lib.c */
      
      /* Compute the canonical spectrum dq:
       *      dq(q,ws) = \sum_{(x,ws)} 
       *        \tilde{WT}(q,x,ws) \ln \tilde{WT}(x,ws)
       */
      if(stq == 0.)	dq[ind_ws] = 0.;
      else   dq[ind_ws] = (double)(q * sTqLogT[ind_ws][ind_q] / stq 
				   - log(stq)/log(10.));
      /* equivalent to the function PFAccessDQFloat */
      
    } /* end loop over the scales: "for (ind_ws=0;..." */
    
    
    /* Regression over the scales to get the different exponents.
     * For each moment, the tq, hq and dq are tabulars of nws values
     * whose slopes give the corresponding exponent Tauq, H and Dh. */
    ss = linefit( tq, LwsArray, nws,  &(Tauq[ind_q]) );
    
    linefit( hq, LwsArray, nws, &(H[ind_q]) );
    
    linefit( dq, LwsArray, nws, &(Dh[ind_q]) );
    
  } /* end loop over the moments: "for (ind_q=0;..." */
  
  if(VERBOSE) 
    printf("%d scales used for the estimation of multifractal exponents\n", ss);
  
  /* Free memory */ 
  if(tq) free(tq);
  if(LwsArray) free(LwsArray);
}



/* ===================================================================== */
int linefit( double *yValues, double *xValues, int dimx, double *slope )  {
  /* ===================================================================== 
   *
   * Fit a signal with a straight line by regression, i.e. it computes the
   * slope a of the line y=ax+b that best fits the data [xValues,yValues].
   *
   * Parameters: 
   *      - xValues : the list of indexes of the signal,
   *      - yValues : the values of the signal to fit, 
   *      - dimx : lenght of the signal,
   *      - slope : slope a of the approximation line, what we want
   * Returns the number of scales used in the approximation
   *
   * ===================================================================== 
   * Note: in LastWave toolbox, this operation is mainly realized by the 
   * function:
   *      LineFitSig
   * in file signal_function.c (package package_signal), except that we
   * consider only signals with regular time intervals.
   * ===================================================================== */

  int i;
  double t,sxoss,sx=0.,sy=0.,st2=0.;
  double a = 0.;
  /* a : the equation line is y = a*x+b */
  double x, y;
  int ss=0;

  for( i=0; i<dimx; i++ ) 
  {
      sx += xValues[i];
      sy += yValues[i];
      ss++;
  }

  sxoss = sx/(double)ss;

  for( i=0; i<dimx; i++ ) 
  {
      t = xValues[i] - sxoss;
      st2 += t*t;
      a += t*yValues[i];
  }
  a /= st2;

  *slope = a;

  return ss;
}


/* ===================================================================== */
int compare(const double *d1,const double *d2) {
  /* ===================================================================== 
   *
   * Stupid function to compare double values and used by qsort in function
   * canonMeanCompute below.
   *
   * ===================================================================== */
  
  if(*d1<*d2)    return -1;
  else if(*d1 == *d2)    return 0;
  else    return +1;
}

void acumula_momentos( double **signal, double **moments)
{
    double partial_mom[Nmom];
    double *grad;
    double eps,eps0,poweps;

    int xeff,dimy;
    int dmax;
    int ix,j,iy,ip;
    int it;


/*     Initialization    */

    dimy=(D_space==1)?1:LEFF;
    grad=(double *) calloc(LEFF,sizeof(double));

    for(it=0;it<Ndist;it++)
    {

/*        Defining the number of intervals of size t that can be included    */

	dmax=LEFF-dist[it]-1;

/*  
     Cleaning the partial moments in which the results at this distance
     are accumulated. We proceed in this way to avoid losing contributions
     which are relevant for the single series under process but what could 
     be lost when compared to all the accumulated course.
*/

	for(ip=0;ip<Nmom;ip++) partial_mom[ip]=0.;

/*   
     We will now run all the (horizontal) ranges of point of lenght t, and for 
     them we evaluate their contribution to the different moments at this 
     size.
*/

 	for(iy=0;iy<dimy;iy++)
	{

/*    Initializating the weight list   */

	    asigna_lista(LEFF,signal[iy],grad);
	    gradiente_1D(LEFF,grad);
	    for(ix=0;ix<LEFF;ix++) grad[ix]=fabs(grad[ix]);

/*    Computing the contributions of all the ranges of size dist[it]  */
/*
      To increase computational speed, to compute the contribution of
      a new interval we update it with respect to the previous one; so,
      we just need to sum the contribution of the new point entering and
      to discount that of the first point in the previous interval
*/

/*      Contribution for ix=0    */

	    eps=0.;
	    for (j=0;j<dist[it];j++) eps+=grad[j];
	    eps/=(double)dist[it];
	    for(ip=0;ip<Nmom;ip++)
	    {
		if(eps>1e-10) partial_mom[ip]+=pow(eps,moms[ip]);
	    }
	    
	    for (ix=0;ix<dmax-1;ix++)
	    {
		eps+=(grad[ix+dist[it]]-grad[ix])/((double)dist[it]);
		for(ip=0;ip<Nmom;ip++)
		{
		    if(eps>1e-10) partial_mom[ip]+=pow(eps,moms[ip]);
		}
	    }
	}
	
/*   Transferring the partial moments to the accumulated moments  */

	for(ip=0;ip<Nmom;ip++) 
	    moments[ip][it]+=partial_mom[ip]/((double)dimy*dmax);

    }

/*         Releasing memory before ending   */

    free(grad);
}

void calcula_taup( double **moments, double *taup)
{
    double lgd[Ndist],rows[Ndist];
    double a,b,corr;
    int ip,it;

    for(ip=0;ip<Nmom;ip++)
    {
	for (it=0;it<Ndist;it++)
	{
	    lgd[it]=log(dist[it]);
	    rows[it]=log(moments[ip][it]);
	}
	fit(&lgd[0],&rows[0],Ndist,&a,&b,&corr);
	taup[ip]=a;
    }

}

int simple_Legendre_transform( double *taup, double **h_m, double **Dh_m,
				double **errDh_m)
{
    int N_m;
    int ip;


/*      Defining sizes and reserving space           */

    N_m=Nmom-1;
    *h_m=(double *) calloc(N_m,sizeof(double));
    *Dh_m=(double *) calloc(N_m,sizeof(double));
    *errDh_m=(double *) calloc(N_m,sizeof(double));


/*        Obtaining the spectrum and its singularities   */

    h_m[0][0]=(taup[1]-taup[0])/(moms[1]-moms[0]);
    Dh_m[0][0]=moms[0]*h_m[0][0]+1-taup[0];
    errDh_m[0][0]=0.2; // conventionally
    for(ip=1;ip<N_m;ip++)
    {
	h_m[0][N_m-1-ip]=(taup[ip+1]-taup[ip-1])/(moms[ip+1]-moms[ip-1]);
	Dh_m[0][N_m-1-ip]=moms[ip]*h_m[0][N_m-1-ip]+1-taup[ip];
	errDh_m[0][N_m-1-ip]=0.2; // conventionally
    }

    return(N_m);
}


void genera_expon( int leff, double **signal, double **expon)
{
    if(D_space==1) genera_expon_1D(leff,signal[0],expon[0]);
    else genera_expon_2D(leff,signal,expon);

}

void genera_expon_1D( int leff, double *serie, double *expon)
{
    double meang;
    int ix;

    asigna_lista(leff,serie,expon);
    gradiente_1D(leff,expon);
    for(ix=0;ix<leff;ix++) expon[ix]=fabs(expon[ix]);
    meang=1.;


    for(ix=0;ix<leff;ix++)
    {
	if(expon[ix]/meang>1e-30) 
	    expon[ix]=-log(expon[ix]/meang)/log((double)leff);
	else expon[ix]=30.*log(10.)/log((double)leff);
    }

}


void genera_expon_2D( int leff, double **image, double **expon)
{
       double **gy;
       double meang;
       int ix,iy;

       gy=reservar_matriz(leff,leff);

       asigna(leff,leff,image,expon);
       gradiente_2D(leff,leff,expon,gy);
       for(iy=0;iy<leff;iy++)
       { 
       for(ix=0;ix<leff;ix++)
       { 
	 expon[iy][ix]=sqrt(expon[iy][ix]*expon[iy][ix]
			      +gy[iy][ix]*gy[iy][ix]);
       }
       }	
       meang=((double)leff);
 
       for(iy=0;iy<leff;iy++)
       {
       for(ix=0;ix<leff;ix++)
       {
	 if(expon[iy][ix]/meang>1e-30) 
	   expon[iy][ix]=-log(expon[iy][ix]/meang)/log((double)leff);
	 else expon[iy][ix]=30.*log(10.)/log((double)leff);
       }
       }
	

       liberar_matriz(gy,leff);
}


void acumula_histograma( int leff, double *mh, double **expon, double *histo)
{
    int dimy;
    int ix,iy,ip;

    dimy=(D_space==1)?1:leff;

    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<leff;ix++)
    {
	ip=(int)((expon[iy][ix]-mh[0])/(mh[1]-mh[0])*NBOX);
	if(ip<0) ip=0;
	if(ip>=NBOX) ip=NBOX-1;
	histo[ip]+=1.;		
    }
    }
}

int calcula_Dh_histo( double sc0, double *mh, double *histo, double **h, 
		      double **Dh, double **errDh)
{
    double *histo_r,*h_r;
    double prob,dp,Nev;
    double maxprob;
    double cump;
    double h0;
    const double min_ev=60;
    const double Ks=3; // 3 sigmas corresponds to a 99.7% confidence level

    int Nh;
    int ip,ir;

	
/*          Inicialization                    */

    histo_r=(double *)calloc(NBOX,sizeof(double));
    h_r=(double *)calloc(NBOX,sizeof(double));

/*        Filtering histogram to avoid low-probability distorsions    */

    cump=0.;
    ir=0;
    for(ip=0;ip<NBOX;ip++)
    {
	cump+=histo[ip];
	histo_r[ir]+=histo[ip]*histo[ip];
	h_r[ir]+=(mh[0]+(0.5+(double)ip)/((double)NBOX)*(mh[1]-mh[0]))
	    *histo[ip];
	if(cump>min_ev)
	{
	    histo_r[ir]/=cump;
	    h_r[ir]/=cump;
	    cump=0.;
	    ir++;
	}
    }

    if(cump>0.)
    {
	histo_r[ir]/=cump;
	h_r[ir]/=cump;
	ir++;	  
    }
    Nh=ir;

/*  The number of exponents is now known; we reservate memory accordingly */

    *h=(double *) calloc(Nh,sizeof(double));
    *Dh=(double *) calloc(Nh,sizeof(double));
    *errDh=(double *) calloc(Nh,sizeof(double));
    

/*                Producing the error bars                         */ 

//  First, we compute the total number of events

    Nev=0.;
    for(ip=0;ip<Nh;ip++) Nev+=histo_r[ip];

// The confidence range is taken as +- Ks sigmas in the distribution 
// of probability boxes, which is a renormalized binomial. In this routine
// we directly propagate the confidence range to the D(h) calculation 

    for(ip=0;ip<Nh;ip++)
    {
	h[0][ip]=h_r[ip]; // associated singularity value

	prob=histo_r[ip]/Nev; // probability of that interval
	if(prob>1e-30) dp=Ks*sqrt((1.-prob)/(prob*Nev)); 
                   // Ks sigmas in the prob. distribution
	else dp=0.;
	

// The error bar is constructed by propagation; as the propagated interval is
// asymetric, we take the bar as the maximum of the two distances to the 
// central value. This is given by the lower bound in the logarithm

	if(dp>=1.) errDh[0][ip]=(double)D_space;
	else errDh[0][ip]=log(1.-dp)/log(sc0);
    } 
   	

/*           Finding and normalizing by the mode                   */

    maxprob=histo_r[0];
    for(ip=1;ip<Nh;ip++) maxprob=fMax(maxprob,histo_r[ip]);
    if(maxprob>1e-30) for(ip=0;ip<Nh;ip++) histo_r[ip]/=maxprob;


/*     Evaluating experimental D(h)         */
    
    for(ip=0;ip<Nh;ip++)
    {
	if(histo_r[ip]>1e-30) 
	    Dh[0][ip]=(double)D_space-log(histo_r[ip])/log(sc0);
	else Dh[0][ip]=0.;
    }

    free(histo_r);
    free(h_r);


    return(Nh);
}

int calcula_Dh_histo_no_ponderado( double sc0, double *mh, double *histo, 
				   double **h, double **Dh, double **errDh)
{
    double *histo_r,*h_r;
    double prob,dp,Nev;
    double maxprob;
    double cump;
    double h0;
    const double min_ev=60;
    const double Ks=3; // 3 sigmas corresponds to a 99.7% confidence level

    int *width;
    int Nh;
    int ip,ip0,ir;

	
/*          Inicialization                    */

    histo_r=(double *)calloc(NBOX,sizeof(double));
    h_r=(double *)calloc(NBOX,sizeof(double));
    width=(int *)calloc(NBOX,sizeof(int));

/*        Filtering histogram to avoid low-probability distorsions    */

    cump=0.;
    ip0=-1;
    ir=0;
    for(ip=0;ip<NBOX;ip++)
    {
	cump+=histo[ip];
	h_r[ir]+=(mh[0]+(0.5+(double)ip)/((double)NBOX)*(mh[1]-mh[0]))
	    *histo[ip];
	if(cump>min_ev)
	{
	    width[ir]=ip-ip0;
	    histo_r[ir]=cump/((double)width[ir]);
	    h_r[ir]/=cump;
	    cump=0.;
	    ip0=ip;
	    ir++;
	}
    }

    if(cump>0.)
    {
	width[ir]=NBOX-1-ip0;
	histo_r[ir]=cump/((double)width[ir]);
	h_r[ir]/=cump;
	ir++;	  
    }
    Nh=ir;

/*  The number of exponents is now known; we reservate memory accordingly */

    *h=(double *) calloc(Nh,sizeof(double));
    *Dh=(double *) calloc(Nh,sizeof(double));
    *errDh=(double *) calloc(Nh,sizeof(double));
    

/*                Producing the error bars                         */ 

//  First, we compute the total number of events

    Nev=0.;
    for(ip=0;ip<Nh;ip++) Nev+=histo_r[ip]*(double)width[ip];

// The confidence range is taken as +- Ks sigmas in the distribution 
// of probability boxes, which is a renormalized binomial. In this routine
// we directly propagate the confidence range to the D(h) calculation 

    for(ip=0;ip<Nh;ip++)
    {
	h[0][ip]=h_r[ip]; // associated singularity value

	prob=histo_r[ip]*(double)width[ip]/Nev; // probability of that interval
	if(prob>1e-30) dp=Ks*sqrt((1.-prob)/(prob*Nev)); 
                   // Ks sigmas in the prob. distribution
	else dp=0.;
	

// The error bar is constructed by propagation; as the propagated interval is
// asymetric, we take the bar as the maximum of the two distances to the 
// central value. This is given by the lower bound in the logarithm

	if(dp>=1.) errDh[0][ip]=(double)D_space;
	else errDh[0][ip]=log(1.-dp)/log(sc0);
    } 
   	

/*           Finding and normalizing by the mode                   */

    maxprob=histo_r[0];
    for(ip=1;ip<Nh;ip++) maxprob=fMax(maxprob,histo_r[ip]);
    if(maxprob>1e-30) for(ip=0;ip<Nh;ip++) histo_r[ip]/=maxprob;


/*     Evaluating experimental D(h)         */

    for(ip=0;ip<Nh;ip++)
    {
	if(histo_r[ip]>1e-30) 
	    Dh[0][ip]=(double)D_space-log(histo_r[ip])/log(sc0);
	else Dh[0][ip]=0.;
    }

    free(histo_r);
    free(h_r);
    free(width);


    return(Nh);
}

int carga_Dh( char *nombre, int *Nh, double **h, double **Dh, double **errDh)
{
    double **series;
    int error;
    int dims,dimt;
    int it;

    lee_serie_temp(nombre,&dims,&dimt,&series);

    if((dims>3)&&(dimt>0))
    {
	error=0;
	*Nh=dimt;
	*h=(double *) calloc(*Nh,sizeof(double));
	*Dh=(double *) calloc(*Nh,sizeof(double));
	*errDh=(double *) calloc(*Nh,sizeof(double));
	for(it=0;it<*Nh;it++)
	{
	    h[0][it]=series[0][it];
	    Dh[0][it]=series[1][it];
	    errDh[0][it]=series[2][it];
	}
    }
    else
    {
	error=1;
	*Nh=0;
    }

/*      Memory release and end     */

    liberar_matriz(series,dims);
    return(error);
}

void lee_serie_temp( char *name_in, int *dims, int *dimt, double ***series)
{
    FILE *canal;
    char campo[90];
    float dato;
    int it,is,pass;

    canal=fopen(name_in,"rt");
    *dims=columnas_serie_temp(canal);
    *dimt=lineas_serie_temp(canal);

    *series=reservar_matriz(*dims,*dimt);

    pass=1;
    for(it=0;(it<*dimt)&&pass;it++)
    {
    for(is=0;(is<*dims)&&pass;is++)
    {
	if(fscanf(canal,"%s",campo)==EOF) pass=0;
	if(pass>0)
	{
	    sscanf(campo,"%f",&dato);
	    series[0][is][it]=(double)dato;
	}
    }
    }
    fclose(canal);

}

int columnas_serie_temp( FILE *canal)
{
    char campo[90];
    int status;
    int dims;
    int pass;

    for(dims=0,pass=1;pass;dims++)
    {
	if(fscanf(canal,"%s",campo)==EOF) pass=0;
	status=fgetc(canal);
	if((status==EOF)||(status==0x0a)) pass=0;       
    }

    return(dims);
}

int lineas_serie_temp( FILE *canal)
{
    char campo[90];
    int status;
    int dimt;
    int pass;

    fseek(canal,0,SEEK_SET);
    for(dimt=0,pass=1;pass;dimt++)
    {
	columnas_serie_temp(canal);
	status=fgetc(canal);
	if(status==EOF) pass=0;
    }
    fseek(canal,0,SEEK_SET);

    return(dimt);
}

double registra_Dh( int graba, int Nr, char *nombre, double sc0, double *h, 
		    double *Dh, double *errDh)
{
    FILE *canal;

    double *errDh_unif;
    double Dh_th;
    double hmin,hmax,hrad;
    double mean_err,av_err,quad_err,sigma,norma;
    double weight;
    double y;
    double h_unif;
    double sc;

    const int Nh=20; // Number of points to be solved in the spectrum
    int Nerr;
    int ix,ip,ip1,ih,ih1,ih2;
	


/*                    Generating outputs               */

    if(graba)
    {
	canal=fopen(nombre,"wt");
	for(ip=0;ip<Nr;ip++)
	{
	    Dh_th=theoretical_Dh(h[ip]);
	    fprintf(canal,"%f  %f  %f  %f\n",h[ip],Dh[ip],errDh[ip],Dh_th);
	}
	fclose(canal);
    }
 
/*      Generating close to uniformly sampled D(h) from data  */
    
    errDh_unif=(double *) calloc(Nh,sizeof(double));
    theoretical_Deltah(&hmin,&hmax);
    hrad=(hmax-hmin)/((double)Nh); // Step and uncertainty radius
    ip=0;
    for(ih=0,h_unif=hmin;ih<Nh;ih++,h_unif+=hrad)
    {
	Nerr=0;
	errDh_unif[ih]=0.;
	if(ip<Nr)
	{
	    for(;(h[ip]<h_unif-hrad/2)&&(ip<Nr);ip++) ;
	    for(;(h[ip]<h_unif+hrad/2)&&(ip<Nr);ip++)
	    {
		Dh_th=theoretical_Dh(h[ip]);
		if((Dh[ip]<=D_space)&&(Dh_th>0.)) // Reasonable points
		{
		    Nerr++;
		    errDh_unif[ih]+=Dh_th-Dh[ip];
		}
	    }
	}
	if(Nerr>0) errDh_unif[ih]/=(double)Nerr;
	else errDh_unif[ih]=theoretical_Dh(h[ip]); 
    }


/*       Calculating errors                */
	
    norma=mean_err=av_err=quad_err=0.;
    for(ih=0,h_unif=hmin+hrad/2;ih<Nh;ih++,h_unif+=hrad)
    {
	weight=1.;
	mean_err+=errDh_unif[ih]*weight;
	av_err+=fabs(errDh_unif[ih])*weight;
	quad_err+=errDh_unif[ih]*errDh_unif[ih]*weight;
	norma+=weight;
    }
    
    mean_err/=norma;
    av_err/=norma;
    sigma=sqrt(quad_err/norma-mean_err*mean_err);
    quad_err=sqrt(quad_err/norma);
  
    if(VERBOSE)
    {
	printf("Mean err: %0.2f; std. deviation: %0.2f\n",mean_err,sigma);
	printf("Typical deviation: %0.2f; quad. error: %0.2f\n",av_err,quad_err);
	printf("\n");
    }

/*         Freeing memory before finishing   */ 

    free(errDh_unif);
    return(av_err/(hmax-hmin));
}

double theoretical_Dh( double h)
{
    double Dh_th;
    double beta,omega,y;

    switch(TYPE)
    {
	case 0:
	    beta=1+HINF/CODINF;
	    omega=-(h-HINF)/(CODINF*log(beta));
	    if(omega>1e-30) 
		Dh_th=(double)D_space+CODINF*(omega-omega*log(omega)-1.);
	    else Dh_th=0;
	    break;
	case 1:
	    Dh_th=(double)D_space-(h-MU)*(h-MU)/(2.*SIGMA*SIGMA);
	    break;
	case 2:
	    Dh_th=D_space-pow(fabs(h-MU)/SIGMA,ALPHA)/ALPHA;
	    break;
	case 3:
	    y=(h-HINF)/(H1-HINF);
	    if((y<0.01)||(y>0.99))  Dh_th=0.;
	    else Dh_th=(double)D_space-1.-(y*log(y)+(1-y)*log(1-y))/log(2.);
	    break;
	default:
	    Dh_th=0.;
	    break;
    }
    return(Dh_th);
}

void theoretical_Deltah( double *hmin, double *hmax)
{
    const double steph=0.01;

    switch(TYPE)
    {
	case 0:
	    *hmin=HINF;
/*
      As the equation is trancendent, we look for the second zero-crossing 
         in an iterative way 
*/
	    for(*hmax=0;theoretical_Dh(*hmax)>0.;*hmax+=steph) ;
	    *hmax-=steph;
	    break;
	case 1:
	    *hmin=MU-SIGMA*sqrt(2.*(double)D_space);
	    *hmax=MU+SIGMA*sqrt(2.*(double)D_space);
	    break;
	case 2:
	    *hmin=MU-SIGMA*pow(ALPHA*(double)D_space,1./ALPHA);
	    *hmax=MU+SIGMA*pow(ALPHA*(double)D_space,1./ALPHA);
	    break;
	case 3:
	    *hmin=HINF;
	    *hmax=H1;
	    break;
	default:
	    *hmin=-1.;
	    *hmax=1.;
	    break;
    }
}

void lee_datos( int leff, char *nombre_in, double **signal)
{
    FILE *canal;
    char nombre[90];
    int iy,dimy;


    dimy=(D_space==1)?1:leff;

    sprintf(nombre,"%s.dat",nombre_in);
    canal=fopen(nombre,"rb");
    for(iy=0;iy<dimy;iy++) fread(signal[iy],sizeof(double),leff,canal);
    fclose(canal);


}

void graba_datos( int leff, char *nombre_in, double **signal)
{
    FILE *canal;
    char nombre[90];
    int iy,dimy;


    dimy=(D_space==1)?1:leff;

    sprintf(nombre,"%s.dat",nombre_in);
    canal=fopen(nombre,"wb");
    for(iy=0;iy<dimy;iy++) fwrite(signal[iy],sizeof(double),leff,canal);
    fclose(canal);

    sprintf(nombre,"%s.hdr",nombre_in);
    canal=fopen(nombre,"wt");
    fprintf(canal,"%s.dat\n%d",nombre_in,leff);
    fclose(canal);


}

void lee_datos_float( int leff, char *nombre_in, double **signal)
{
    FILE *canal;
    char nombre[90];
    float val;
    int ix,iy,dimy;


    dimy=(D_space==1)?1:leff;
 
    sprintf(nombre,"%s.fdat",nombre_in);
    canal=fopen(nombre,"rb");
    for(iy=0;iy<dimy;iy++)
    { 
    for(ix=0;ix<leff;ix++)
    { 
	fread(&val,sizeof(float),1,canal);
	signal[iy][ix]=(double)val;
    }
    }
    fclose(canal);


}

void graba_datos_float( int leff, char *nombre_in, double **signal)
{
    FILE *canal;
    char nombre[90];
    float val;
    int ix,iy,dimy;


    dimy=(D_space==1)?1:leff;
 
    sprintf(nombre,"%s.fdat",nombre_in);
    canal=fopen(nombre,"wb");
    for(iy=0;iy<dimy;iy++)
    { 
    for(ix=0;ix<leff;ix++)
    { 
	val=(float)signal[iy][ix];
	fwrite(&val,sizeof(float),1,canal);
    }
    }
    fclose(canal);

    sprintf(nombre,"%s.hdr",nombre_in);
    canal=fopen(nombre,"wt");
    fprintf(canal,"%s.dat\n%d",nombre_in,leff);
    fclose(canal);


}

void graba_serie( int leff, char *nombre_in, double *datos)
{
    FILE *canal;
    char nombre[90];
    int ix;

    sprintf(nombre,"%s.txt",nombre_in);
    canal=fopen(nombre,"wt");
    for(ix=0;ix<leff;ix++)
	fprintf(canal,"%f\n",datos[ix]);
    fclose(canal);

}


void graba_imagen( int leff, char *nombre_in, double **datos)
{
    char nombre[90];

    sprintf(nombre,"%s.gif",nombre_in);
    graba_foto(leff,leff,nombre,datos);

}

void graba_geomapa( char *nombre, double **Q)
{
    FILE *canal;
    int ileffs,inums;

    canal=fopen(nombre,"wt");
    for(inums=0;inums<Nnums;inums++)
    {
	for(ileffs=0;ileffs<Nleffs;ileffs++)
	    fprintf(canal,"%f  ",Q[inums][ileffs]);
	fprintf(canal,"\n");
    }
    fclose(canal);
}

void graba_typemap( char *nombre, double **Q)
{
    FILE *canal;
    int dim1,dim2;
    int i1,i2;

    canal=fopen(nombre,"wt");
    if(TYPE==0)
    {
	i1=0;
	dim1=1;
	dim2=NLPs;
	fprintf(canal,"%f  ",Q[i1][0]);
 	for(i2=1;i2<dim2;i2++)
	{
	    if(codinfs[i2-1]>codinfs[i2]) fprintf(canal,"\n");
	    fprintf(canal,"%f  ",Q[i1][i2]);
	}
    }
    else
    {
	dim1=Nmeans;
	dim2=Nsigmas;
	for(i1=0;i1<dim1;i1++)
	{
	    for(i2=0;i2<dim2;i2++)
		fprintf(canal,"%f  ",Q[i1][i2]);
	    fprintf(canal,"\n");
	}
    }

    fclose(canal);
}

double moda(int leff, double **signal)
{
    if(D_space==1) return(moda_1D(leff,signal[0]));
    else return(moda_2D(leff,leff,signal));
}

double moda_1D(int dimx, double *datos)
{
    double *histo;
    double mm[2];
    double out;
    int Nhisto;
    int ix,ip;

    Nhisto=dimx/30; // 30 events by bin in average
    histo=(double *) calloc(Nhisto,sizeof(double));


    extrema_lista(dimx,datos,&mm[0]);
    for(ix=0;ix<dimx;ix++)
    {
	ip=(int)(((double)Nhisto)*(datos[ix]-mm[0])/(mm[1]-mm[0]));
	if(ip>Nhisto-1) ip=Nhisto-1;
	histo[ip]+=1.;
    }

    out=moda_por_histo(Nhisto,mm,histo);
 
    free(histo);
    return(out);
}

double moda_2D(int dimx, int dimy, double **datos)
{
    double *histo;
    double mm[2];
    double out;
    int Nhisto;
    int ix,iy,ip;

    Nhisto=(dimx*dimy)/30; // 30 events by bin in average
    histo=(double *) calloc(Nhisto,sizeof(double));


    extrema(dimx,dimy,datos,&mm[0]);
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	ip=(int)(((double)Nhisto)*(datos[iy][ix]-mm[0])/(mm[1]-mm[0]));
	if(ip>Nhisto-1) ip=Nhisto-1;
	histo[ip]+=1.;
    }
    }

    out=moda_por_histo(Nhisto,mm,histo);


    free(histo);
    return(out);
}

double moda_por_histo( int Nhisto, double *mm, double *histo)
{
    double out;
    int ip,ip0;

    for(ip0=0,ip=0;ip<Nhisto;ip++) if(histo[ip]>histo[ip0]) ip0=ip;
    out=mm[0]+(mm[1]-mm[0])*(0.5+(double)ip0)/((double)Nhisto);

    return(out);
}
