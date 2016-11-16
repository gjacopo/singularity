
/*      Program multifractal_generator.c. Version: 22/09/2005   */
 
 
#include <stdio.h>
#include <math.h>
#include <string.h>
 
                               

//      My libraries

#include <multifractal.c>  // Version del  24 de Agosto, 2004


/*      Global variables               */

int VERBOSE=0;     // Flag controlling the verbose information
int NSERIES=1;     // Number of output series
int D_space=1;     // Dimension of the space
int LMAX=512;     // Size of the generated series
int OUTRES=0;      // Number of resolutions to be elliminated to smoothen. 
int TYPE=0;        // Type of multifractal to be calculated.
                   // 0: Log-Poisson
                   // 1: Log-Normal
                   // 2: Log-Levi
                   // 3: Binomial
int INV_TRANS=0;   // Enables implement explicit translational invariance
float HINF=-0.5;   // Most singular exponent (Types 0 and 3)
float CODINF=1.;   // Most singular co-dimension (Type 0)
float H1=0.5;      // Upper singularity bound in binomial model
float MU=0.5;      // Average singularity (for log-Normal and log-Levi)
float SIGMA=1.;    // Singularity dispersion in log-Normal and log-Levi 
float ALPHA=1.5;   // Exponent defining log-Levi

float TCH=3.;      // Maximum effective dispersion (in sigmas)



int WAV_BASE=0;    // Wavelet of choice.
                   // 0: Gaussian
                   // 1: Lorentian
                   // 2: Haar
int DER_WAV_BASE=2.;// Derivative order of the wavelet of choice


int NBOX=128;       // Number of bins in the histogram


/*         Non-externally accessed variables   */
 

const float DISP_ETA=0.5;     // Log-Normal processing variables

const float D0[3]= {1., .5,  1.}; // Minimum scale for each wavelet.
float DH=0.1;      // Accepted conventional error

double *prob_levi; // To keep a copy of a Levi distribution of parameter alpha

/*      Parametros de compatibilidad interplataforma        */

#define RAND_PER 2147483647 // Substituye a RAND_MAX para la rutina random();

/*      Prototipos de Funciones         */
 
int main(int argc, char *argv[]);
void parse_arguments(int argc, char *argv[]);

void genera_levi(void);
void genera_multifractal( int leff, double **signal);
void genera_multifractal_1D( int leff, double *serie);
void genera_multifractal_2D( int leff, double **image);
void genera_alpha( int dimalpha, int Nwav, double *alpha);
void genera_alpha_1D( int dimalpha, int Nwav, double *alpha0);
void genera_alpha_2D( int dimalpha, int Nwav, double *alpha0);
double genera_h( void);
double poisson( double lambda);
double normal_standard(void);
double levi_standard(void);
double binomial( double h0, double h1);


void wavelet_1D( int leff, double sc, double *wave);
void wavelet_2D( int leff, double sc, double **wave);

void lee_datos( int leff, char *nombre_in, double **signal);
void graba_datos( int leff, char *nombre_in, double **signal);
void lee_datos_float( int leff, char *nombre_in, double **signal);
void graba_datos_float( int leff, char *nombre_in, double **signal);
void graba_serie( int leff, char *nombre_in, double *datos);
void graba_imagen( int leff, char *nombre_in, double **datos);

double moda(int leff, double **signal);
double moda_1D(int dimx, double *datos);
double moda_2D(int dimx, int dimy, double **datos);

int main(int argc, char *argv[])
{

/*	PARAMETERS		*/


/*	DATA			*/

    char nombre[90],base[90];
    
/*	QUANTITIES TO BE COMPUTED	*/

    double **signal;
	
/*	AUXILIAR VARIABLES	*/

    FILE *canal;
    double val,beta,hmax;
    double maxprob;
    double sc;
	
    int block,leff,dimy;
    int in,ip,ir;



/*		Program		*/


    parse_arguments(argc,argv);


/*      Initializating variables and memory           */


    if(TYPE==2)
    {
	genera_levi();
	if(ALPHA<0.01)
	{
	    if(VERBOSE) printf("Too small alpha value. Quantizing...\n");
	    ALPHA=0.01;
	}
    }
    if(TYPE==3)
    {
	if(H1<=HINF)
	{
	    H1=HINF+DH;
	    printf("Unacceptable greater exponent; changing it to %f\n",
		       H1);
	}
    }
    
/*        Constraint of translational invariance for log-Normal    */
    
    if(INV_TRANS)
    {
	switch(TYPE)
	{
	    case 0:
		if(VERBOSE)
			printf("Log-Poisson MF are scale invariant, flag ignored at this point\n");
		break;
	    case 1:
		MU=SIGMA*SIGMA/4.;
		if(VERBOSE)
		    printf("Mean singularity changed to %f\n",MU);
		break;
	    case 2:
		if(fabs(ALPHA-1.)>1e-30)
		{
		    MU=pow(SIGMA/ALPHA,1.+1./(ALPHA-1.))*(ALPHA-1.);
		    if(VERBOSE)
			printf("Mean singularity changed to %f\n",MU);
		}
		else if(VERBOSE) printf("Scale invariance impossible to implement for this MF\n");
		break;
	    case 3:
		if(HINF>1)
		{
		    H1=-log(1-pow(0.5,HINF));
		    if(VERBOSE) 
			printf("Maximum singularity changed to %0.2f\n",H1);
				       
		}
		else
		{
		    if((H1>0.)&&(H1<1.))
		    {
			H1=0.75;
			printf("Maximum singularity changed to %0.2f\n",H1);
		    }
		    HINF=-log(1-pow(0.5,H1));
		    printf("Minimum singularity changed to %0.2f\n",HINF);
		}
		break;
	}
    }
    
/*		Data reading		*/


    if(VERBOSE)
	printf("Processing for %d output %dD signals\n",NSERIES,D_space);
    
    block=(int)pow(2.,(double)OUTRES);
    leff=dimensiona(LMAX)/block;
    dimy=(D_space==1)?1:leff;

    if(leff<2)
    {
	printf("Too reduced resolution!!!\n");
	exit(-1);
    }

    signal=reservar_matriz(dimy,leff);

    switch(TYPE)
    { 
	case 0:
	    sprintf(base,"Log-Poisson.h%0.2f-coD%0.2f-size%d",
		    HINF,CODINF,leff);
	    if(VERBOSE)
		printf("of the type Log-Poisson; hinf= %0.2f, Dinf= %0.2f\n",
		       HINF,(double)D_space-CODINF);
	    break;
	case 1:
	    sprintf(base,"Log-Normal.mean%0.2f-sigma%0.2f-size%d",
		    MU,SIGMA,leff);
	    if(VERBOSE)
		printf("of the type Log-normal; mean= %0.2f, sigma= %0.2f\n",
		       MU,SIGMA);
	    break;
	case 2:
	    sprintf(base,"Log-Levi_mean%0.2f-sigma%0.2f-alpha%0.2f-size%d",
		    MU,SIGMA,ALPHA,leff);
	    if(VERBOSE)
		printf("of the type Log-Levi; mean= %0.2f, sigma= %0.2f, alpha= %0.2f\n",
		       MU,SIGMA,ALPHA);
	    break;
	case 3:
	    sprintf(base,"Binomial_h0%0.2f-h1%0.2f-size%d",HINF,H1,leff);
	    if(VERBOSE)
		printf("of the type binomial; h0= %0.2f, h1= %0.2f\n",
		       HINF,H1);
	    break;
	default:
	    sprintf(base,"bugged");
	    break;
    }
    if(D_space==1) strcat(base,"_1D");
    else strcat(base,"_2D");
	
    for(in=0;in<NSERIES;in++)
    {
	genera_multifractal(block*leff,signal);
	sprintf(nombre,"%s-N%05d",base,in);
	graba_datos(leff,nombre,signal);
	if(VERBOSE)
	{
	    if(D_space==1) graba_serie(leff,nombre,signal[0]);
	    else graba_imagen(leff,nombre,signal);
	}
    }
    liberar_matriz(signal,dimy);
    
    
/*	FREEING MEMORY BEFORE FINISHING	*/


    if(TYPE==2) free(prob_levi);
    
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
	int Narg0=30;  // Initialization value; it should be greater than (but not necessarily equal to) the number or arguments
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

//    Argument LMAX. Type 1: integer

	sprintf(olarg[in],"%s","-dim");
	sprintf(olval[in],"%s","size");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Linear size for generated series. It will be rounded to\nthe least power of 2 greater than this value. Default:",LMAX); 
	type[in]=1;
	ptrvar_i[in]=&LMAX;
	ptrval_i[in][0]=2;
	ptrval_i[in][1]=1000000;
	in++;

//    Argument OUTRES. Type 1: integer

	sprintf(olarg[in],"%s","-outres");
	sprintf(olval[in],"%s","#res_levels");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Number of dyadic resolutions to be smoothened. Default:",
		OUTRES); 
	type[in]=1;
	ptrvar_i[in]=&OUTRES;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=10;
	in++;

//    Argument INV_TRANS. Type 0: flag

	sprintf(olarg[in],"%s","-invtrans");
	sprintf(olexp[in]," %s : %s\n",
		olarg[in],
		"Flag. If enabled, the program produces translational invariant multifractals.\nDefault: DISABLED"); 
	type[in]=0;
	ptrflag[in]=&INV_TRANS;
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
		"Type of multifractal to be generated.\n   0: Log-Poisson\n   1: Log-Normal\n   2: Log-Levi\n   3: Binomial\nDefault:",TYPE); 
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

//    Argument TCH. Type 2: float

	sprintf(olarg[in],"%s","-max_disp");
	sprintf(olval[in],"%s","#sigmas");
	sprintf(olexp[in]," %s : %s %0.2f\n",
		olarg[in],
		"Singularity range, expressed in sigmas. Valid for log-Normal and log-Levi.\nDefault:",TCH); 
	type[in]=2;
	ptrvar_f[in]=&TCH;
	ptrval_f[in][0]=1.;
	ptrval_f[in][1]=20.;
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

//   Derivative parameters (included in <derivacion_1D.c>)

	in=parsing_derivacion(in,0,0,olarg,olval,olexp,ptrvar_f,ptrval_f,
			  ptrvar_i,ptrval_i,ptrflag,type);

//   multifractal parameters (included in <multifractal_1D.c>)

	in=parsing_multifractal_1D(in,0,0,olarg,olval,olexp,ptrvar_f,ptrval_f,
			  ptrvar_i,ptrval_i,ptrflag,type);


//    Argument WAV_BASE. Type 1: integer

	sprintf(olarg[in],"%s","-wv_basis");
	sprintf(olval[in],"%s","choice");
	sprintf(olexp[in],"%s\n %s : %s %d\n",
		"\nPARAMETERS DEFINING THE REPRESENTATION BASIS\n===========================================",
		olarg[in],
		"Wavelet of choice for the basis.\n    0: Gaussian wavelet\n    1: Lorentzian wavelet\n    2: Diagonal Haar\nDefault:",WAV); 
	type[in]=1;
	ptrvar_i[in]=&WAV_BASE;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=2;
	in++;

//    Argument DER_WAV_BASE. Type 1: integer

	sprintf(olarg[in],"%s","-der_wv_basis");
	sprintf(olval[in],"%s","order");
	sprintf(olexp[in]," %s : %s %d\n",
		olarg[in],
		"Order of the derivatives in the wavelet. Default:",DER_WAV_BASE); 
	type[in]=1;
	ptrvar_i[in]=&DER_WAV_BASE;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=2;
	in++;


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

void genera_levi( void)
{
    double *prob_leviI;
    double x,dx,norma;
    int ix;

    prob_levi=(double *) calloc(NBOX,sizeof(double));
    prob_leviI=(double *) calloc(NBOX,sizeof(double));

    for(ix=0;ix<NBOX;ix++)
    {
	x=(double)ix;
	if(ix>NBOX/2) x-=(double)NBOX;
	dx=fabs(x*PI/(2.*TCH));
	if(dx>1e-30) prob_levi[ix]=exp(-ALPHA*pow(dx,ALPHA));
	prob_leviI[ix]=prob_levi[ix]*sin(PI*x);
	prob_levi[ix]=prob_levi[ix]*cos(PI*x);
    }
    FFFT1D(NBOX,prob_levi,prob_leviI,-1);
    for(ix=0,norma=0.;ix<NBOX;ix++) norma=fMin(norma,prob_levi[ix]);
    for(ix=0;ix<NBOX;ix++) prob_levi[ix]-=norma;

    for(ix=0,norma=0.;ix<NBOX;ix++) norma+=prob_levi[ix];
    for(ix=0;ix<NBOX;ix++) prob_levi[ix]/=norma;


    free(prob_leviI);
}


void genera_multifractal( int leff, double **signal)
{
    if(D_space==1) genera_multifractal_1D(leff,signal[0]);
    else genera_multifractal_2D(leff,signal);
}

void genera_multifractal_1D( int leff, double *serie)
{
       double *wave,*aux,*serie2;
       double *alpha;
       double block;
       int dimalpha,linfalpha,Nwav;
       int dima,bl,db;
       int ix,ip,iwav;

       
       Nwav=adimensiona(leff)-1;
       dimalpha=leff-1;
       block=pow(2.,(double)OUTRES);

       alpha=(double *) calloc(dimalpha,sizeof(double));
       wave=(double *) calloc(leff,sizeof(double));
       aux=(double *) calloc(leff,sizeof(double));
       serie2=(double *) calloc(leff,sizeof(double));


       genera_alpha(dimalpha,Nwav,alpha);
       for(iwav=0,dima=1,linfalpha=0;iwav<=Nwav;iwav++,dima*=2)
       {
	   wavelet_1D(leff,(double)dima,wave);
	   limpia_lista(leff,aux);
	   bl=leff/dima;
	   db=bl/2;
	   for(ix=0;ix<dima;ix++) aux[bl*ix+db]=alpha[linfalpha+ix];
	   convuelto_1D(leff,wave,aux);
	   for(ix=0;ix<leff;ix++) serie2[ix]+=aux[ix];
	   linfalpha+=dima;
       }
       coarse_resolution_lista(leff,block,serie2,serie);

       free(alpha);
       free(wave);
       free(serie2);
       free(aux);
}

void genera_multifractal_2D( int leff, double **image)
{
       double **wave,**aux,**image2;
       double *alpha;
       double block;
       int dimalpha,linfalpha,Nwav;
       int dima,bl,db;
       int ix,iy,iwav;

       
       Nwav=adimensiona(leff)-1;
       dimalpha=(leff*leff-1)/3;
       block=pow(2.,(double)OUTRES);

       alpha=(double *) calloc(dimalpha,sizeof(double));
       wave=reservar_matriz(leff,leff);
       aux=reservar_matriz(leff,leff);
       image2=reservar_matriz(leff,leff);


       genera_alpha(dimalpha,Nwav,alpha);
       for(iwav=0,dima=1,linfalpha=0;iwav<=Nwav;iwav++,dima*=2)
       {
	 wavelet_2D(leff,(double)dima,wave);
	 limpia(leff,leff,aux);
	 bl=leff/dima;
	 db=bl/2;
	 for(iy=0;iy<dima;iy++)
	 { 
	 for(ix=0;ix<dima;ix++)
	 { 
	   aux[bl*iy+db][bl*ix+db]=alpha[linfalpha+ix+iy*dima];
	 }
	 }
	 convuelto_2D(leff,leff,wave,aux);
	 asigna_suma(leff,leff,aux,image2);

	 linfalpha+=dima*dima;
       }
       coarse_resolution(leff,leff,block,image2,image);

       free(alpha);
       liberar_matriz(wave,leff);
       liberar_matriz(image2,leff);
       liberar_matriz(aux,leff);
}

void genera_alpha( int dimeta, int Nwav, double *alpha)
{
    if(D_space==1) genera_alpha_1D(dimeta,Nwav,alpha);
    else genera_alpha_2D(dimeta,Nwav,alpha);
}

void genera_alpha_1D( int dimeta, int Nwav, double *alpha0)
{
    double *alpha,*lalpha;
    double lalphaexp;
    double heff;
    double signo;
    int linfalpha,linfalpha0;
    int dima,dima0;
    int ix,ix0,iwav;
    int ip;

    alpha=(double *) calloc(dimeta,sizeof(double));
    lalpha=(double *) calloc(dimeta,sizeof(double));


    alpha[0]=1.;
    lalpha[0]=0.;
    linfalpha=0;
    dima=1;
    for(iwav=1;iwav<=Nwav;iwav++)
    {
	
/*         Update of parameters              */ 

	linfalpha0=linfalpha;
	linfalpha+=dima;
	dima0=dima;
	dima*=2;
	    
	 
/*       Generation                 */

	for(ix=0;ix<dima;ix++)
	{
	    ix0=ix/2;
		
/*        Producing the accompanying multiplicative factor     */

	    signo=((double)random())/RAND_PER;
	    if(signo<0.5) signo=-1.;
	    else signo=1.;
	    alpha[linfalpha+ix]=signo*sqrt(0.5)*alpha[linfalpha0+ix0];
	    

/*       Producing the logarithm of alphas            */
	    
	    lalpha[linfalpha+ix]=genera_h()+lalpha[linfalpha0+ix0];
		
/*      Accumulating alphas      */

	    alpha0[linfalpha+ix]+=
		alpha[linfalpha+ix]*exp(-log(2.)*lalpha[linfalpha+ix]);
		
	}

/*     End of scale loop        */

    }
 
    free(lalpha);
    free(alpha);


}

void genera_alpha_2D( int dimeta, int Nwav, double *alpha0)
{
    double *alpha,*lalpha;
    double lalphaexp;
    double heff;
    double signo;
    int linfalpha,linfalpha0;
    int dima,dima0;
    int ix,iy,ix0,iy0,iwav;
    int ip;

    alpha=(double *) calloc(dimeta,sizeof(double));
    lalpha=(double *) calloc(dimeta,sizeof(double));


    alpha[0]=1.;
    lalpha[0]=0.;
    linfalpha=0;
    dima=1;
    for(iwav=1;iwav<=Nwav;iwav++)
    {

/*         Update of parameters              */ 
	
	linfalpha0=linfalpha;
	linfalpha+=dima*dima;
	dima0=dima;
	dima*=2;

	 
/*       Generation                 */

	for(iy=0;iy<dima;iy++)
	{
	    iy0=iy/2;
	    for(ix=0;ix<dima;ix++)
	    {
		ix0=ix/2;
		
/*        Producing the accompanying multiplicative factor     */
		
		signo=((double)random())/RAND_PER;
		if(signo<0.5) signo=-1.;
		else signo=1.;
		
		alpha[linfalpha+ix+iy*dima]=
		    signo*alpha[linfalpha0+ix0+iy0*dima0];
		
/*       Producing the logarithm of alphas            */
		
		lalpha[linfalpha+ix+iy*dima]=genera_h()+
		    lalpha[linfalpha0+ix0+iy0*dima0];
		
		
/*      Accumulating alphas      */
		
		alpha0[linfalpha+ix+iy*dima]+=alpha[linfalpha+ix+iy*dima]*
		    exp(-log(2.)*lalpha[linfalpha+ix+iy*dima]);
		
	    }
	    
	}

/*     End of scale loop        */

    }


    free(alpha);
    free(lalpha);


}

double genera_h( void)
{
       double h;
       double beta,s;
       double prob;
       int ip;
       
       switch(TYPE)
       {
	   case 0:
	       beta=1+HINF/CODINF;
	       s=CODINF*log(2.);
	       h=HINF-poisson(s)*log(beta)/log(2.);
	       break;
	   case 1:
	       h=MU+SIGMA*normal_standard()/sqrt(log(2.));
	       break;
	   case 2:
	       h=MU+SIGMA*levi_standard()*pow(log(2.),-1./ALPHA);
	       break;
	   case 3:
	       h=binomial(HINF,H1);
	       break;
	   default:
	       h=0.;
	       break;
       }

 

       return(h);
}

double poisson( double lambda)
{
       double tirada,prob,weight;
       int ip;


       tirada=((double)random())/RAND_PER;
       weight=exp(-lambda);
       ip=0;

       for(prob=weight;((prob<tirada)&&(ip<100));ip++,prob+=weight)
	 weight*=lambda/((double)(ip+1));
  
       return((double)ip);
}

double normal_standard( void)
{
       double norma,dx;
       double tirada,prob;
       double x;
       int ib;
 
       dx=2.*TCH/((double)NBOX);
       norma=0.;
       x=-TCH;
       for(ib=0;ib<NBOX;ib++,x+=dx) norma+=exp(-0.5*x*x);
       norma=1./norma;


       tirada=((double)random())/RAND_PER; 
       x=-TCH;
       prob=0.;
       for(ib=0;(ib<NBOX)&&(prob<=tirada);ib++,x+=dx)
	 prob+=norma*exp(-0.5*x*x);
 
       return(x);

}

double levi_standard(void)
{
       double norma;
       double tirada,prob,x,dx;
       int ib;
 
       tirada=((double)random())/RAND_PER;

       x=-TCH;
       dx=2.*TCH/((double)NBOX);
       prob=0.;
       for(ib=0;(ib<NBOX)&&(prob<tirada);ib++,x+=dx)
	   prob+=prob_levi[ib];

 
       return(x);
}


double binomial( double h0, double h1)
{
    double tirada;
    double h;

    tirada=((double)random())/RAND_PER;
    if(tirada<0.5) h=h0;
    else h=h1;

    return(h);
}


void wavelet_1D( int leff, double sc, double *wave)
{
       double norma;
       double x;
       int ix,ix0,ix1;

       for(ix=0;ix<leff;ix++)
       {
	 x=((double)ix)/((double)leff);
	 if(x>0.5) x-=1.;
	 x=x*sc*D0[WAV_BASE];
	 switch(WAV_BASE)
	 {
	    case 0:
	      switch(DER_WAV_BASE)
	      {
	         case 0:
		   wave[ix]=exp(-0.5*x*x);
		   break;
	         case 1:
		   wave[ix]=-x*exp(-0.5*x*x);
		   break;
	         case 2:
	         default:
		   wave[ix]=(x*x-1.)*exp(-0.5*x*x);
		   break;
	      }
	      break;
	    case 1:
	      switch(DER_WAV_BASE)
	      {
	         case 0:
		   wave[ix]=1./(1.+x*x);
		   break;
	         case 1:
		   wave[ix]=-x/pow(1.+x*x,2.);
		   break;
	         case 2:
	         default:
		   wave[ix]=(1.-x*x)/pow(1.+x*x,3.);
		   break;
	      }
	      break;
	    case 2:
	    default:
	      switch(DER_WAV_BASE)
	      {
	         case 0:
		   if(fabs(x)>0.5) wave[ix]=0.;
		   else if(x>0) wave[ix]=1.;
		   else wave[ix]=-1.;
		   break;
	         case 1:
	         case 2:
	         default:
		   wave[ix]=0.;
		   break;
	      }
	      break;
	 }
       }
       
       if(WAV_BASE==2)
       {
	 ix0=Mod((int)( ((double)leff)/(2*sc*D0[WAV_BASE]))-1,leff);
	 ix1=Mod(leff-1-ix0,leff);
	 if(DER_WAV_BASE==1)
	 {
	   wave[0]=1.;
	   wave[ix0]=wave[ix1]=-0.5;
	 }
	 else if(DER_WAV_BASE==2)
	 {
	   wave[0]=-1.;
	   wave[leff-1]=1.;
	   wave[ix0]=wave[ix1]=0.5;
	   wave[Mod(ix0-1,leff)]=wave[Mod(ix1-1,leff)]=-0.5;
	 }
       }

       norma=0.;
       for(ix=0;ix<leff;ix++) norma+=wave[ix]*wave[ix];
       norma=sqrt(norma/((double)leff));
       for(ix=0;ix<leff;ix++) wave[ix]/=norma;
}

void wavelet_2D( int leff, double sc, double **wave)
{

       double x,y;
       int ix,iy;

       for(iy=0;iy<leff;iy++)
       {
	 y=((double)iy)/((double)leff);
	 if(y>0.5) y-=1.;
	 y=y*sc*D0[WAV_BASE];
	 for(ix=0;ix<leff;ix++)
	 {
	   x=((double)ix)/((double)leff);
	   if(x>0.5) x-=1.;
	   x=x*sc*D0[WAV_BASE];
	   switch(WAV_BASE)
	   {
	      case 0:
		wave[iy][ix]=exp(-0.5*(x*x+y*y));
		break;
	      case 1:
		wave[iy][ix]=1./(1.+x*x+y*y);
		break;
	      case 2:
	      default:
		if((fabs(x)>0.5)||(fabs(y)>0.5)) wave[iy][ix]=0.;
		else if(x*y>0) wave[iy][ix]=1.;
		else wave[iy][ix]=-1.;
		break;
	   }
       }
       }
       if(DER_WAV_BASE) filtro_2D(leff,leff,(double)DER_WAV_BASE,wave);

       anorma2(leff,leff,wave);
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
    int ix,ip,ip0;

    Nhisto=dimx/30; // 30 events by bin in average
    histo=(double *) calloc(Nhisto,sizeof(double));


    extrema_lista(dimx,datos,&mm[0]);
    for(ix=0;ix<dimx;ix++)
    {
	ip=(int)(((double)Nhisto)*(datos[ix]-mm[0])/(mm[1]-mm[0]));
	if(ip>Nhisto-1) ip=Nhisto-1;
	histo[ip]+=1.;
    }

    for(ip0=0,ip=0;ip<Nhisto;ip++) if(histo[ip]>histo[ip0]) ip0=ip;
    out=mm[0]+(mm[1]-mm[0])*(0.5+(double)ip0)/((double)Nhisto);



    free(histo);
    return(out);
}

double moda_2D(int dimx, int dimy, double **datos)
{
    double *histo;
    double mm[2];
    double out;
    int Nhisto;
    int ix,iy,ip,ip0;

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

    for(ip0=0,ip=0;ip<Nhisto;ip++) if(histo[ip]>histo[ip0]) ip0=ip;
    out=mm[0]+(mm[1]-mm[0])*(0.5+(double)ip0)/((double)Nhisto);



    free(histo);
    return(out);
}
