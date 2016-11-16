/*        multifractal.c. Version del 28 de Septiembre, 2005   */

#define MULTIFRACTAL_C

#ifndef STRUCT_DEF
#include <struct_def.c>
#endif

#ifndef DERIVACION_C
#include <derivacion.c> 
#endif


#ifndef GRAFICOS_C
#include <graficos.c>
#endif

/*           Run-time variables           */

/*       Common to 1 and 2D analyses    */

int NPOINTS=6;      // Number of scales in the multifractal evaluation
int WAV=0;          // wavelet of choice:
                    // 0: gaussian
                    // 1: Lorentzian, order 0.5
                    // 2: Lorentzian, order 1
                    // 3: Lorentzian, order 1.5

int ORD_DER=0;      // Differentiation order
int HOLDER=0;       // Analizes de measure (0) or the own function
float S0=1.;        // Initial scale for the wavelets
float WAV_RANGE=64.; // Range of scales in wavelet analysis

/*          Specific to 2D analysis                        */

int DIM1=0;         // Flag. If set, the analysis is performed in 1D
float THETAU=0;     // Direction for the 1D analysis

/*  
minimum scales (SC0) and scale steps (QS). They depend on the derivative 
order (first argument) and on the wavelet, defined by WAV (second argument) 
*/

/* wavelet exponents 1D (negative for gaussian) */

const double EXP_WL_1D[4] = { -1.  , 0.5  , 1.,  1.5   }; 

const double SC0_1D[3][4] =  
{ { 1.000, 1.000, 1.000, 1.000}, // Derivative order 0
  { 2.000, 2.000, 2.000, 2.000}, // Derivative order 1
  { 4.000, 4.000, 4.000, 4.000}  // Derivative order 2
}; 


/* wavelet exponents 2D  (negative for gaussian) */

const double EXP_WL_2D[4] = { -1.  , 1.,  1.5  , 2.  }; 

const double SC0_2D[3][4] =  
{ { 1.000, 0.500, 0.500, 1.0000}, // Derivative order 0
  { 2.000, 2.000, 2.000, 2.0000}, // Derivative order 1
  { 4.000, 4.000, 4.000, 4.0000}  // Derivative order 2
}; 


/*          Function prototypes      */

int parsing_multifractal_1D( int in0, int siflag, int *deflag, 
			     char **olarg, char **olval, char **olexp,
			     float **ptrvar_f, float **ptrval_f, 
			     int **ptrvar_i, int **ptrval_i,
			     int **ptrflag, int *type);


double calcula_multifractal_1D( int dimx, int verb, double *signal, 
				double *expon, double *mm_h);
double escala_wavelet_1D( int dimx);
void genera_wavelet_1D( int dimx, double sc, int ord_der, double expon, 
			double *wave);
void define_wavelet_1D( int dimx, double sc, int ord_der, double expon, 
		     double *wave);
void normaliza_wavelet_1D( int dimx, double sc, int ord_der, double *wave);




int parsing_multifractal_2D( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type);

double escala_lineal_2D( int dimx, int dimy);
double escala_wavelet_2D( int dimx, int dimy);
double escala_wavelet_2D_lineas( int dimx, int dimy);

void modderiva( int dimx, int dimy, double **modg);
void modderiva_1D( int dimx, int dimy, double **modg);


double calcula_multifractal_2D( int dimx, int dimy, int verb, double **signal,
				double **expon, double *mm_h);
void genera_wavelet_2D( int dimx, int dimy, double sc, double **wave);
void define_wavelet_2D( int dimx, int dimy, double sc, double **wave);
void genera_wavelet_1D_2D( int dimx, int dimy, double sc, double **wave);
void define_wavelet_1D_2D( int dimx, int dimy, double sc, double **wave);
void normaliza_wavelet_2D( int dimx, int dimy, double sc, double **wave);

void anorma1_lineas( int dimx, int dimy, double **data);



void visualiza_color( int dimx, int dimy, double blout, char *dext, 
		      double h0, double h1, double **expon);
void paleta( double rep, char *Red, char *Green, char *Blue);
void visualiza_gris_trozos( int dimx, int dimy, double blout, char *dext,
			    double h0, double h1, double delta_h, 
			    double **expon);


/*     Function declarations     */


int parsing_multifractal_1D( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type)
{
        int in;

	in=in0;


//   Argument WAV. Type 1: Integer

	strcpy(olarg[in],"-wav");
	strcpy(olval[in],"wav_index");
	sprintf(olexp[in],"%s\n %s : %s %d\n",
		"\nWAVELET VARIABLES\n=================",
		olarg[in],
		"Wavelet to be used.\n   0: Gaussian\n   1: Lorentzian at exponent 0.5\n   2: Lorentzian\n   3: Lorentzian at 1.5.\nDefault: "
		,WAV);
	if(siflag)
	{
	  type[in]=4;
	  ptrflag[in]=deflag;
	}
	else type[in]=1;
	ptrvar_i[in]=&WAV;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=4;
	in++;

//   Argument ORD_DER. Type 1: Integer

	strcpy(olarg[in],"-der");
	strcpy(olval[in],"deriv_order");
	sprintf(olexp[in]," %s : %s %d\n",olarg[in],
		"Derivative order for the wavelet. Default: ",ORD_DER);
	if(siflag)
	{
	  type[in]=4;
	  ptrflag[in]=deflag;
	}
	else type[in]=1;
	ptrvar_i[in]=&ORD_DER;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=2;
	in++;

//   Argument HOLDER (<multifractal_2D.c>). Type 0: Flag

	strcpy(olarg[in],"-hold");
	sprintf(olexp[in]," %s : %s\n",olarg[in],
		"Flag. It enables function processing (in opposition to\nmeasure analysis). Default: DISABLED");
	if(siflag)
	{
	  type[in]=3;
	  ptrflag[in]=deflag;
	  ptrvar_i[in]=&HOLDER;
	}
	else
	{
	  type[in]=0;
	  ptrflag[in]=&HOLDER;
	}
	in++;

//   Argument S0 (<multifractal_2D.c>). Type 2: Float

	strcpy(olarg[in],"-s0");
	strcpy(olval[in],"scale_0");
	sprintf(olexp[in]," %s : %s %f\n",olarg[in],
		"Initial scale for wavelets. Default: ",S0);
	if(siflag)
	{
	  type[in]=5;
	  ptrflag[in]=deflag;
	}
	else type[in]=2;
	ptrvar_f[in]=&S0;
	ptrval_f[in][0]=.33;
	ptrval_f[in][1]=100.;
	in++;

//   Argument WAV_RANGE . Type 2: Float

	strcpy(olarg[in],"-range");
	strcpy(olval[in],"scale_rng");
	sprintf(olexp[in]," %s : %s %f\n",olarg[in],
		"Range of scales in wavelet analysis. Default: ",S0);
	if(siflag)
	{
	  type[in]=5;
	  ptrflag[in]=deflag;
	}
	else type[in]=2;
	ptrvar_f[in]=&WAV_RANGE;
	ptrval_f[in][0]=1.5;
	ptrval_f[in][1]=1000.;
	in++;

	return(in);
}



double calcula_multifractal_1D( int dimx, int verb, double *signal, 
				double *expon, double *mm_h)
{
    double **wave;
    double *dmasa;
    double *xx,*yy;
    double sc,sc0;
    double a,b,corr;
    double qs_1D;
	
    int xeff;
    int Nbuen;
    int ix,ip;

    xeff=mdimensiona(dimx);
    
    wave=reservar_matriz(NPOINTS+1,xeff);
    dmasa=(double *) calloc(xeff,sizeof(double));
    xx=(double *) calloc(NPOINTS+1,sizeof(double));
    yy=(double *) calloc(NPOINTS+1,sizeof(double));
    
    asigna_lista(dimx,signal,dmasa);
    if(HOLDER==0)
    {
	gradiente_1D(xeff,dmasa);
	dmasa[dimx-1]=dmasa[xeff-1]=0.; // to avoid boundary jumps
	for(ix=0;ix<xeff;ix++) dmasa[ix]=fabs(dmasa[ix]);
	}
    
/*   Introducing the derivative order  */

    for(ip=0;ip<ORD_DER;ip++)
    {
	if(verb) printf("Applying the %d-th derivative\n",ip+1);
	gradiente_1D(xeff,dmasa);
    }
	
/*    Calculating the wavelet projections   */

    sc0=escala_wavelet_1D(xeff);
    qs_1D=pow(WAV_RANGE,1./((double)NPOINTS));
    for(ip=0,sc=S0*SC0_1D[ORD_DER][WAV];ip<=NPOINTS;ip++,sc*=qs_1D)
    {
	xx[ip]=1./log(sc0*sc);
	genera_wavelet_1D(xeff,sc,ORD_DER,EXP_WL_1D[WAV],wave[ip]);
	convuelto_1D(xeff,dmasa,wave[ip]);
    }

    Nbuen=0;
    for(ix=0;ix<dimx;ix++)
    {
	for(ip=0;ip<=NPOINTS;ip++)
	{
	    a=fabs(wave[ip][ix]);
	    if(a>1e-30) yy[ip]=log(a);
	    else yy[ip]=-30.*log(10.);
	    yy[ip]*=xx[ip];
	}
	fit(xx,yy,NPOINTS+1,&a,&b,&corr);
	if(fabs(corr)>0.9) Nbuen++;
	expon[ix]=b;
	mm_h[0]=fMin(mm_h[0],b);
	mm_h[1]=fMax(mm_h[1],b);
    }
    free(xx);
    free(yy);
    
    if(verb)
    {
	printf("Multifractal analysis\n");
	printf("=====================\n");
	switch(WAV)
	{
	    case 0:
		printf("Gaussian wavelet\n");
		break;
	    default:
		printf("Lorentzian wavelet of exponent %f\n",EXP_WL_1D[WAV]);
		break;
		
	}
	if(ORD_DER) printf("in %d-th derivative\n",ORD_DER);
	printf("Singularities: minimum: %f and maximum: %f\n",
	       mm_h[0],mm_h[1]);
	printf("Percentage of good regression points: %f %%\n",
	       100.*((double)Nbuen)/((double)dimx));
    }
    
    liberar_matriz(wave,NPOINTS+1);
    free(dmasa);
    
    
    return(((double)Nbuen)/((double)dimx));
}

double escala_wavelet_1D( int dimx)
{
	double *wave;
        double sc;
	double norma;
	int xeff;
	int ix;

	xeff=mdimensiona(dimx);

	wave=(double *) calloc(xeff,sizeof(double));
	define_wavelet_1D(xeff,SC0_1D[ORD_DER][WAV],ORD_DER,EXP_WL_1D[WAV],wave);
	for(ix=0,norma=0.;ix<xeff;ix++) norma+=wave[ix]*wave[ix];

// Like that, norma is assumed to be a diameter. It follows

	sc=sqrt(norma/dimx);

	free(wave);

	return(sc);
}


void genera_wavelet_1D( int dimx, double sc, int ord_der, double expon, 
	double *wave)
{

        define_wavelet_1D(dimx,sc,ord_der,expon,wave);
	normaliza_wavelet_1D(dimx,sc,ord_der,wave);
}

void define_wavelet_1D( int dimx, double sc, int ord_der, double expon, 
	double *wave)
{
	double x,valor;
	int ix;

	for(ix=0;ix<dimx;ix++)
	{
	  x=(double)ix;
	  if(ix>=dimx/2) x-=(double)dimx;
	  x=x/sc;
	  if(expon>0.) valor=pow(1.+x*x,-expon);
	  else  valor=exp(-0.5*x*x);
	  wave[ix]=valor;		
	}
	
}

void normaliza_wavelet_1D( int dimx, double sc, int ord_der, double *wave)
{
    const double D_space=1;
    double norma;
    int ix;
    
    //for(ix=0,norma=0.;ix<dimx;ix++) norma+=fabs(wave[ix]);
    //norma*=pow(sc,D_space-ord_der);
    norma=pow(sc,D_space-ord_der);
    for(ix=0;ix<dimx;ix++) wave[ix]/=norma;

}

int parsing_multifractal_2D( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type)
{
        int in;

	in=in0;

//   Argument DIM1. Type 0: Flag

	strcpy(olarg[in],"-dim1");
	sprintf(olexp[in],"%s\n %s : %s\n",
		"\nMULTIFRACTAL VARIABLES\n======================",olarg[in],
		"Flag. If enabled, analysis is carried out in 1D.\nDefault: DISABLED");
	if(siflag)
        {
	  type[in]=3;
	  ptrvar_i[in]=&DIM1;
	  ptrflag[in]=deflag;
	}
	else
	{
	  type[in]=0;
	  ptrflag[in]=&DIM1;
	}
	in++;

	
//   Argument THETAU. Type 2: Float

	strcpy(olarg[in],"-theta");
	strcpy(olval[in],"angle_1D");
	sprintf(olexp[in]," %s : %s %f\n",olarg[in],
		"Multifractal variable. Direction for 1D analysis. Default:",THETAU);
	if(siflag)
	{
	  type[in]=5;
	  ptrflag[in]=deflag;
	}
	else type[in]=2;
	ptrvar_f[in]=&THETAU;
	ptrval_f[in][0]=0;
	ptrval_f[in][1]=PI/2.;
	in++;

//   Argument WAV. Type 1: Integer

	strcpy(olarg[in],"-wav");
	strcpy(olval[in],"wav_index");
	sprintf(olexp[in],"%s\n %s : %s %d\n",
		"\nWAVELET VARIABLES\n=================",
		olarg[in],
		"Wavelet to be used.\n   0: Gaussian\n   1: Lorentzian at exponent 0.5\n   2: Lorentzian\n   3: Lorentzian at 1.5.\nDefault: "
		,WAV);
	if(siflag)
	{
	  type[in]=4;
	  ptrflag[in]=deflag;
	}
	else type[in]=1;
	ptrvar_i[in]=&WAV;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=4;
	in++;

//   Argument ORD_DER. Type 1: Integer

	strcpy(olarg[in],"-der");
	strcpy(olval[in],"deriv_order");
	sprintf(olexp[in]," %s : %s %d\n",olarg[in],
		"Derivative order for the wavelet. Default: ",ORD_DER);
	if(siflag)
	{
	  type[in]=4;
	  ptrflag[in]=deflag;
	}
	else type[in]=1;
	ptrvar_i[in]=&ORD_DER;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=2;
	in++;

//   Argument HOLDER (<multifractal_2D.c>). Type 0: Flag

	strcpy(olarg[in],"-hold");
	sprintf(olexp[in]," %s : %s\n",olarg[in],
		"Flag. It enables function processing (in opposition to\nmeasure analysis). Default: DISABLED");
	if(siflag)
	{
	  type[in]=3;
	  ptrflag[in]=deflag;
	  ptrvar_i[in]=&HOLDER;
	}
	else
	{
	  type[in]=0;
	  ptrflag[in]=&HOLDER;
	}
	in++;

//   Argument S0 (<multifractal_2D.c>). Type 2: Float

	strcpy(olarg[in],"-s0");
	strcpy(olval[in],"scale_0");
	sprintf(olexp[in]," %s : %s %f\n",olarg[in],
		"Initial scale for wavelets. Default: ",S0);
	if(siflag)
	{
	  type[in]=5;
	  ptrflag[in]=deflag;
	}
	else type[in]=2;
	ptrvar_f[in]=&S0;
	ptrval_f[in][0]=.33;
	ptrval_f[in][1]=100.;
	in++;

	return(in);
}


double escala_lineal_2D( int dimx, int dimy)
{
        return(1./sqrt((double)dimx*dimy));
}


double escala_wavelet_2D( int dimx, int dimy )
{
	double **wave;
        double sc;
	double norma;
	double x,y,r;
	int ix,iy;


	wave=reservar_matriz(dimy,dimx);
	define_wavelet_2D(dimx,dimy,SC0_2D[ORD_DER][WAV],wave);
	norma=media(dimx,dimy,wave);
	sc=sqrt(norma/PI);

	liberar_matriz(wave,dimy);

	return(sc);
}


double escala_wavelet_2D_lineas( int dimx, int dimy)
{
        double **wave;
	double x,y,mod,proy;
	double ux,uy;
        double sc;
	int lmax;
	int ix,iy;
	      
	wave=reservar_matriz(dimy,dimx);
	define_wavelet_1D_2D(dimx,dimy,SC0_2D[ORD_DER][WAV],wave);


	ux=cos(THETAU);
	uy=sin(THETAU);

	sc=0.;
	lmax=0;
	for(iy=0;iy<dimy;iy++)
	{
	  y=(double)iy;
	  if((iy>0)&&(iy>=dimy/2)) y-=(double)dimy;
	  for(ix=0;ix<dimx;ix++)
	  {
	    x=(double)ix;
	    if((ix>0)&&(ix>=dimx/2)) x-=(double)dimx;
	    mod=sqrt(x*x+y*y);
	    if(mod<1e-30)
	    {
	      lmax++;
	      sc+=1.;
	    } 
	    else
	    {
	      proy=fabs(uy*x-ux*y)/mod;
	      if(proy<1e-10)
	      {
		lmax++;
		sc+=fabs(wave[iy][ix]);
	      }
	    }
	}
	}
	sc/=(double)lmax;

	liberar_matriz(wave,dimy);
	return(sc);
}

void modderiva( int dimx, int dimy, double **modg)
{
        double **aux;
        int ix,iy;

	aux=reservar_matriz(dimy,dimx);

	gradiente_2D(dimx,dimy,modg,aux);
	for(iy=0;iy<dimy;iy++)
	{ 
	for(ix=0;ix<dimx;ix++)
	{ 
	  modg[iy][ix]=sqrt(modg[iy][ix]*modg[iy][ix]+aux[iy][ix]*aux[iy][ix]); 
	}
	}
	liberar_matriz(aux,dimy);
}

void modderiva_1D( int dimx, int dimy, double **modg)
{
        double **aux;
	double ux,uy;
	int ix,iy;

	ux=cos(THETAU);
	uy=sin(THETAU);

	aux=reservar_matriz(dimy,dimx);

	gradiente_2D(dimx,dimy,modg,aux);
	for(iy=0;iy<dimy;iy++)
	{ 
	for(ix=0;ix<dimx;ix++)
	{ 
	  modg[iy][ix]=fabs(ux*modg[iy][ix]+uy*aux[iy][ix]); 
	}
	}
	liberar_matriz(aux,dimy);
}

double calcula_multifractal_2D( int dimx, int dimy, int verb, double **signal,
				double **expon, double *mm_h)
{
    double ***wave;
    double **dmasa;
    double *xx,*yy;
    double sc,sc0;
    double a,b,corr;
    double qs_2D;

    int xeff,yeff;
    int Nbuen;
    int ix,iy,ip;
    
    xeff=mdimensiona(dimx);
    yeff=mdimensiona(dimy);
    

    wave=reservar_tritensor(NPOINTS+1,yeff,xeff);
    dmasa=reservar_matriz(yeff,xeff);
    xx=(double *) calloc(NPOINTS+1,sizeof(double));
    yy=(double *) calloc(NPOINTS+1,sizeof(double));
    
    
    asigna(dimx,dimy,signal,dmasa);
    if(HOLDER==0)
    {
	if(DIM1) modderiva_1D(xeff,yeff,dmasa);
	else modderiva(xeff,yeff,dmasa);
    }
    if(ORD_DER) filtro_2D(xeff,yeff,(double)ORD_DER,dmasa);
    
/*    Calculating the wavelet projections   */

    if(DIM1) sc0=escala_wavelet_2D_lineas(dimx,dimy);
    else sc0=escala_wavelet_2D(dimx,dimy);

    qs_2D=pow(WAV_RANGE,1./((double)NPOINTS));
    for(ip=0,sc=S0*SC0_2D[ORD_DER][WAV];ip<=NPOINTS;ip++,sc*=qs_2D)
    {
	xx[ip]=1./log(sc*sc0);
	if(DIM1) genera_wavelet_1D_2D(xeff,yeff,sc,wave[ip]);
	else genera_wavelet_2D(xeff,yeff,sc,wave[ip]);
	
	convuelto_2D(xeff,yeff,dmasa,wave[ip]);
	anorma1_lineas(xeff,yeff,wave[ip]);
    }

    Nbuen=0;
    for(iy=0;iy<dimy;iy++)
    {
	for(ix=0;ix<dimx;ix++)
	{
	    for(ip=0;ip<=NPOINTS;ip++)
	    {
		a=fabs(wave[ip][iy][ix]);
		if(a>1e-30) yy[ip]=log(a);
		else yy[ip]=-30.*log(10.);
		yy[ip]*=xx[ip];
	    }
	    fit(xx,yy,NPOINTS+1,&a,&b,&corr);
	    if(fabs(corr)>0.9) Nbuen++;
	    expon[iy][ix]=b;
	    mm_h[0]=fMin(mm_h[0],b);
	    mm_h[1]=fMax(mm_h[1],b);
	}
    }
    free(xx);
    free(yy);
    
    if(verb)
    {
	printf("Multifractal analysis\n");
	printf("=====================\n");
	switch(WAV)
	{
	    case 0:
		printf("Gaussian wavelet\n");
		break;
	    case 1:
		printf("Lorentzian wavelet\n");
		break;
	    case 2:
		printf("Lorentzian wavelet of exponent 1.5\n");
	    default:
		printf("Lorentzian wavelet of exponent 2\n");
		break;
		
	}
	if(ORD_DER) printf("in %d-th derivative\n",ORD_DER);
	printf("Singularities: minimum: %f and maximum: %f\n",
	       mm_h[0],mm_h[1]);
	printf("Percentage of good regression points: %f %%\n",
	       100.*((double)Nbuen)/((double)dimx*dimy));
    }
    
    liberar_tritensor(wave,NPOINTS+1,yeff);
    liberar_matriz(dmasa,yeff);
    
    
    return(((double)Nbuen)/((double)dimx*dimy));
}

void genera_wavelet_2D( int dimx, int dimy, double sc, double **wave)
{
        define_wavelet_2D(dimx,dimy,sc,wave);
	normaliza_wavelet_2D(dimx,dimy,sc,wave);
}


void define_wavelet_2D( int dimx, int dimy, double sc, double **wave)
{
	double x,y,valor;
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	  y=(double)iy;
	  if((iy>0)&&(iy>=dimy/2)) y-=(double)dimy;
	  y/=sc;
	  for(ix=0;ix<dimx;ix++)
	  {
	    x=(double)ix;
	    if((ix>0)&&(ix>=dimx/2)) x-=(double)dimx;
	    x/=sc;
	    if(EXP_WL_2D[WAV]>0.) valor=pow(1.+x*x+y*y,-EXP_WL_2D[WAV]);
	    else  valor=exp(-0.5*(x*x+y*y));
	    wave[iy][ix]=valor;		
	  }
	}
}

void genera_wavelet_1D_2D( int dimx, int dimy, double sc, double **wave)
{
        define_wavelet_1D_2D(dimx,dimy,sc,wave);
	normaliza_wavelet_2D(dimx,dimy,sc,wave);
}

void define_wavelet_1D_2D( int dimx, int dimy, double sc, double **wave)
{
	double x,y,mod,proy,valor;
	double ux,uy;
	int ix,iy;

	ux=cos(THETAU);
	uy=sin(THETAU);

		
	for(iy=0;iy<dimy;iy++)
	{
	  y=(double)iy;
	  if((iy>0)&&(iy>=dimy/2)) y-=(double)dimy;
	  y/=sc;
	  for(ix=0;ix<dimx;ix++)
	  {
	    x=(double)ix;
	    if((ix>0)&&(ix>=dimx/2)) x-=(double)dimx;
	    x/=sc;
	    mod=sqrt(x*x+y*y);
	    if(mod<1e-30) wave[iy][ix]=1.;
	    else
	    {
	      proy=fabs(uy*x-ux*y)/mod;
	      if(proy<1e-10)
	      {
		if(EXP_WL_2D[WAV]>0.) valor=pow(1.+mod*mod,-EXP_WL_2D[WAV]);
		else  valor=exp(-0.5*(mod*mod));
		wave[iy][ix]=valor;
	      }
	      else wave[iy][ix]=0.;
	    }
	  }
	}
}


void normaliza_wavelet_2D( int dimx, int dimy, double sc, double **wave)
{
        double norma;
	int ix,iy;

	norma=0.;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  norma+=fabs(wave[iy][ix]);
	}
	}
	norma*=pow(sc,-ORD_DER);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  wave[iy][ix]/=norma;
	}
	}

}

void anorma1_lineas( int dimx, int dimy, double **data)
{
	double norma;
	int ix,iy;

	
/*          Version incompleta: solo vale para lineas horizontales
   	                y verticales                                      */

	if(fabs(THETAU)<0.01)
	{
	  for(iy=0;iy<dimy;iy++)
	  {
	    norma=0.;
	    for(ix=0;ix<dimx;ix++) norma+=fabs(data[iy][ix]);
	    norma/=(double)dimx;
	    for(ix=0;ix<dimx;ix++) data[iy][ix]/=norma;
	  }
	}
	else if(fabs(THETAU-PI/2)<0.01)
	{		
	  for(ix=0;ix<dimx;ix++) 
	  {
	    norma=0.;
	    for(iy=0;iy<dimy;iy++) norma+=fabs(data[iy][ix]);
	    norma/=(double)dimy;
	    for(iy=0;iy<dimy;iy++) data[iy][ix]/=norma;
	  }
	}
	else anorma1(dimx,dimy,data);

}


void visualiza_color( int dimx, int dimy, double blout, char *dext, 
		      double h0, double h1, double **expon)
{
        char nombre[90];
	char **Red,**Green,**Blue;
	double **exponeff;
	double rep;
	int dimxf,dimyf;
	int ix,iy;

	dimxf=blout*dimx;
	dimyf=blout*dimy;


	Red=reservar_matriz_char(dimyf,dimxf);
	Green=reservar_matriz_char(dimyf,dimxf);
	Blue=reservar_matriz_char(dimyf,dimxf);
	exponeff=reservar_matriz(dimyf,dimxf);

	coarse_resolution(dimx,dimy,1./blout,expon,exponeff);

	for(iy=0;iy<dimyf;iy++)
	{
	for(ix=0;ix<dimxf;ix++)
	{
	  rep=(exponeff[iy][ix]-h0)/(h1-h0);
	  if(rep<0.) rep=0.;
	  if(rep>1.) rep=1.;
	  paleta(rep,&(Red[iy][ix]),&(Green[iy][ix]),&(Blue[iy][ix]));
	}
	}
	sprintf(nombre,"sing_total.%s.ppm",dext);
	graba_ppm(dimxf,dimyf,1,nombre,Red,Green,Blue);


	liberar_matriz_char(Red,dimyf);
	liberar_matriz_char(Green,dimyf);
	liberar_matriz_char(Blue,dimyf);
	liberar_matriz(exponeff,dimyf);

}

void paleta( double rep, char *Red, char *Green, char *Blue)
{
  double rep2;
  int gamma;
  int scR,scG,scB;
  int ic;

  gamma=0;

  rep2=5.*rep;
  ic=(int) rep2;
  if(ic==5) ic=4;
  rep2-=(double)ic;
  rep2=1.-rep2;


  switch(ic)
  {
      case 0:
	scR=255-gamma;
	scG=255-gamma;
	scB=(255-gamma)*rep2;
	break;
      case 1:
	scR=(255-gamma)*rep2;
	scG=(255-gamma)/2;
	scB=0.;
	break;
      case 2:
	scR=(255-gamma)*(1-rep2);
	scG=0.;//(255-gamma)*rep2;
	scB=0.;
	break;
      case 3:
	scR=(255-gamma)*rep2;
	scG=0.;
	scB=(255-gamma)*(1.-rep2);
	break; 
      case 4:
	scR=0.;
	scG=0.;
	scB=(255-gamma)*rep2;
	break;
  }
  *Red=(char)(gamma+scR);
  *Green=(char)(gamma+scG);
  *Blue=(char)(gamma+scB);


}


void visualiza_gris_trozos( int dimx, int dimy, double blout, char *dext,
			    double h0, double h1, double delta_h, 
			    double **expon)
{
        char nombre[90];
	char ***manifold,***excluded;
	double singmax;
	int Nsing;
	int ix,iy,ip;


	Nsing=(h1-h0)/delta_h;

	manifold=reservar_tritensor_char(Nsing,dimy,dimx);
	excluded=reservar_tritensor_char(2,dimy,dimx);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
        {
	      for(ip=0;ip<Nsing;ip++) manifold[ip][iy][ix]=(char)255;
	      excluded[0][iy][ix]=(char)255; 
	      excluded[1][iy][ix]=(char)255; 
	      
	      ip=(int) floor((expon[iy][ix]-h0)/delta_h);
	      if((ip>=0)&&(ip<Nsing)) manifold[ip][iy][ix]=(char)0;
	      if(expon[iy][ix]<h0-delta_h) excluded[0][iy][ix]=(char)0;
	      if(expon[iy][ix]>singmax-delta_h) excluded[1][iy][ix]=(char)0;

	}
	}


/*	The manifolds are represented	*/

	for(ip=0;ip<Nsing;ip++)
	{
	  sprintf(nombre,"manifold_%02d.%s.gif",ip,dext);
	  graba_mask_block(dimx,dimy,blout,nombre,manifold[ip]);
	}

	sprintf(nombre,"excluded-.%s.gif",dext);
	graba_mask_block(dimx,dimy,blout,nombre,excluded[0]);
	sprintf(nombre,"excluded+.%s.gif",dext);
	graba_mask_block(dimx,dimy,blout,nombre,excluded[1]);

/*     Freeing memory before closign the loop    */

	liberar_tritensor_char(excluded,2,dimy);
	liberar_tritensor_char(manifold,Nsing,dimy);
}

