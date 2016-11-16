/********************************************************/
/*                                                      */
/*                multifractal_wtmm.c                   */
/*           Version del 2 de Dekembriou, 2004          */
/*                                                      */
/*  Main calls to the WTMM functions.                   */
/*                                                      */
/* List of functions:                                   */
/* analiza_series_WTMM: launches the process            */
/* estima_Dh_WTMM: calls the function                   */
/* checkWTMMparameters                                  */
/* checkdimension                                       */
/*                                                      */
/********************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

/* Personal libraries */
#include <const_wtmm.h>
#include <utils.h>
#include <extrema.h>
#include <approx_legendre.h>

#include <multifractal_wtmm.h>


extern int VERBOSE;
extern int NSERIES;

/* Variables related to the WTMM method */
extern int WTMMMethod;
extern int IsBinary;

extern int ord_der_wtmm;
extern int wav_wtmm; 
extern char wav_name[][20];

extern float MinScale, MaxScale;
extern float ScStep;
extern int NVoices;
extern int ScTime; 
extern float ScRatio; 

extern float MinMoment;
extern float MaxMoment;
extern float QStep;
extern float qDefArray[];
extern int nqDef;

extern  float SCmin, SCmax;

extern float ShiftSpectrum;

extern int flagML, flagPF, flagTauq;

extern FILE *mlfile;

#ifdef FLAG_MSM
extern void limpia(int dimx, int dimy, double **matriz);
extern void lee_datos( int leff, char *nombre_in, double **signal);
extern double **reservar_matriz( int ydim, int xdim);
extern void liberar_matriz( double **m, int ydim);
extern float **reservar_matriz_float( int ydim, int xdim);
extern void liberar_matriz_float( float **m, int ydim);
extern int **reservar_matriz_int( int ydim, int xdim);
extern void liberar_matriz_int( int **m, int ydim);
extern char **reservar_matriz_char( int ydim, int xdim);
extern void liberar_matriz_char( char **m, int ydim);
#endif


/* ===================================================================== */
void analiza_series_WTMM( int leff, int dimy, char *base) {
  /* ===================================================================== 
   *
   * Launches the WTMM analysis on file which basename is base.
   * 
   * ===================================================================== */
  
  char nombre[90];
  
  double *tauq, *h, *Dh;
  double *qArray;
  int nq;

  if(VERBOSE) printf("\n\nWTMM Estimation - Main Parameters"
		     "\n=================================");
	   
  /* Check the dimension of the file */
  leff=checkdimension(base,leff);
 
  /* Check/update the parameters passed to the program */
  checkWTMMparameters(leff);

   /** Internal parameter **/
  /* Number of moments */
  if( QStep != 0.) 
    nq = (int)((MaxMoment-MinMoment)/QStep)+1;
  else
    nq = nqDef;
  /* Note that (nq-2) will also be the number of estimated values for
   * the spectrum */

   /** Memory allocations **/
  if(/* List of computed moments  */
     (qArray=(double*)calloc(nq, sizeof(double))) == NULL ||
     /* Matrices for the multifractal exponents and spectra */  
     (tauq = (double*)calloc(nq,sizeof(double))) == NULL ||
     (h = (double*)calloc(nq,sizeof(double))) == NULL ||
     (Dh = (double*)calloc(nq,sizeof(double))) == NULL ) {
    fprintf(stderr,"\n Error while allocating memory in analiza_series_WTMM");
    exit(-1);
  }

  if(VERBOSE) printf("\n\nWTMM Estimation - Main Computation"
		     "\n=================================");
  /** Main computation related to CWT and WTMM **/
  estima_Dh_WTMM( leff, dimy, base, nq, qArray, tauq, h, Dh );

  /** Registration **/
  if(VERBOSE) printf("\n\nWTMM Estimation - Registration"
		     "\n=================================");

  sprintf(nombre,"Dh_wtmm_%s",base);
  if(VERBOSE) printf("\n Writing in %s...", nombre);
  /* Eventually shift to compare it with other spectrum */
  if(ShiftSpectrum)  asigna_desplaza_lista(nq-1,  (double)ShiftSpectrum, Dh);
  registra_WTMM( 1, nq-1, nombre, h, Dh );

  if(flagTauq) {
    sprintf(nombre,"Tq_wtmm_%s",base);
    if(VERBOSE) printf("\n Writing in %s...", nombre);
    registra_WTMM( 0, nq, nombre, qArray, tauq );
  }

  /** Free memory */
  if(tauq) free(tauq);
  if(h) free(h);
  if(Dh) free(Dh);
  if(qArray) free(qArray);
}
  

/* ===================================================================== */
void estima_Dh_WTMM( int dimx, int dimy, char *base,
		     int nq, double *qArray, 
		     double *tauq, double*h, double*Dh ) {
  /* ===================================================================== 
   *
   * Main computation - See file extrema.c and approx_legendre.c for a full
   *                    description of parameters
   * Main results are stored in the tabular tauq, h and Dh of multifractal
   * features.
   *
   * ===================================================================== */
  
  char nombre[90];
  double **signal; /* note: we are not yet able to deal with 2D signals */
  int ind_q, ind_sc;
  double **ExtWTlis, **Z;
  double *scArray;
  int *n_ext, nsc=0; 
  int  icolor, nn;
  double** sTq, **sTqLogT /* **logSTq */;
  int in,ix;
  double sc;

  FILE *fp;
  double temp;

   /** Useful parameters **/
  /* Number of scales */
  for (sc=MinScale; sc<=MaxScale; sc*=ScStep)    nsc++;
  /* Maximum size of wavelet filter */
  nn = (int)(MaxScale* ScTime);

  /** Memory allocations **/
  if( /* original signal */
     (signal = reservar_matriz( dimy, dimx+1 )) == NULL ||
     /* List of scales under consideration */
     (scArray=(double*)calloc(nsc, sizeof(double))) == NULL || 
     /* Matrice of wavelet coefficients */
     (ExtWTlis = reservar_matriz( nsc, dimx+1 )) == NULL || 
     /* Matrice of the coefficients of the partition function 
      *   Z(a,q) = \sum_{L \in MaxLines} 
      *           [ sup_{(x,a') \in L} || */
     (Z = reservar_matriz( nsc, nq+1 )) == NULL ||
     /* Number of maxima for each scale */
     (n_ext = (int*)calloc(nsc+1,sizeof(int))) == NULL) {
    fprintf(stderr,"\n Error while allocating memory in estima_Dh_WTMM");
    exit(-1);
  }

  if( WTMMMethod == CANONICAL ) /* Intermedary internal parameters used for 
				 * the approximation of Legendre transform in 
				 * canonical method */
    if( (sTq = reservar_matriz( nsc,  nq)) == NULL ||
	(sTqLogT = reservar_matriz( nsc, nq)) == NULL ) {
      fprintf(stderr,"\n Error while allocating memory in estima_Dh_WTMM");   
      exit(-1); 
    } 

  /* Allocation for representation of Wavelet Transform and Maxima Lines */
  /* if(flagWT) wtcolor=(double*)calloc((dimx-2*nn-1)*nsc+1, sizeof(double)); */
  if(flagML) mlfile=fopen(MLFILE, "wt");

  /** Initialize the scales and the moments for which estimantions will
   ** be made **/
  if( QStep != 0. ) /* Range of moments */
    for (ind_q=0; ind_q<nq; ind_q++) 
      qArray[ind_q] = MinMoment + (double)ind_q*QStep;
  else /* List of moments by default */
    for (ind_q=0; ind_q<nq; ind_q++) 
      qArray[ind_q] = qDefArray[ind_q];
  for (ind_sc=1, scArray[0]=MinScale; ind_sc<nsc; ind_sc++) 
    scArray[ind_sc] = scArray[ind_sc-1] * ScStep;

  if(VERBOSE) {
    printf("\n Compute Wavelet Transform, Wavelet Extrema and Partition Function"
	   "\n for all signals and over all scales.                             ");    
    if(NSERIES > 1) 
      printf("\n The statistics of the different signals are accumulated.");
  }

  /** Compute the multifractal exponents **/
  for(in=0;in<NSERIES;in++)  {
    
    /* Note: if the lists of moments are strickly equal for the different
     * series, then the result is the addition of the partition function 
     * for each q: we are simply adding statistics (they just have to be 
     * the result of the same type of analysis - same WT on the same range 
     * of scale on signals of the same size).
     */

    if(NSERIES == 1 && IsBinary != BINARY) 
      sprintf(nombre,"%s",base); /* 1 file submitted only */
    else             sprintf(nombre,"%s-N%03d",base,in);
    
    if(IsBinary == BINARY) lee_datos(dimx,nombre,signal); /* binary files */
    else                   lee_serie(dimx,nombre,signal); /* text files */

   /* clean here the WTMM first : ExtWTlis = 0 everywhere, even if we
     * rewrite on it  */
    if(in>0) limpia( dimx+1, nsc, ExtWTlis );

   /* Generate and track maxima lines for extrema extraction */
    WTMM_PF_compute( signal[0], dimx, nn, 
		     (double)wav_wtmm, (double)ord_der_wtmm, 
		     nq, qArray, nsc, scArray, ExtWTlis, Z, n_ext );
    
    if ( WTMMMethod == CANONICAL ) 
      canonMeanCompute( ExtWTlis, dimx, n_ext, nq, qArray, nsc, scArray,
			sTq, sTqLogT /* logSTq */ );
    
    /* Representation of the Wavelet Transform */
    /* if(flagWT) registra_WTcolor(WTFILE, wtcolor, icolor, dimx-2*nn-1, nsc ); */
  } /* end of the loop "for(in=0;in<NSERIES;..." */
  
  
  /** Finally compute the exponents tau(q) and the spectrum (h,D(h)) from the
   ** values of the WT over the detected extrema **/
  if ( WTMMMethod == DIRECT ) {
    if(VERBOSE) 
      printf("\n Approximate Legendre Transform for estimating multifractal exponents"
	     "\n with DIRECT method");
    directSpecCompute( Z, dimx, nq, qArray, nsc, scArray, 
		       LOG(SCmin), LOG(SCmax), tauq, h, Dh );
    
  } else if ( WTMMMethod == CANONICAL ) {
  if(VERBOSE)
      printf("\n Approximate Legendre Transform for estimating multifractal exponents"
	     "\n with CANONICAL method");
    canonSpecCompute( sTq, sTqLogT, /* logSTq */
		      dimx, nq, qArray, nsc, scArray, LOG(SCmin), LOG(SCmax), 
		      tauq, h, Dh );
  }

  if(flagPF)
    registra_PF(PFFILE, Z, nq, qArray, nsc, scArray); 
     /* graba_PF( PFFILE, Z, nq, nsc, scArray );*/

  /** Free memory */
  if(signal) liberar_matriz(signal,dimy);
  if(Z) liberar_matriz(Z,nsc);
  if(ExtWTlis) liberar_matriz(ExtWTlis,nsc);
  if(n_ext) free(n_ext);
  if( WTMMMethod == CANONICAL ) {
    liberar_matriz(sTq,nsc);
    liberar_matriz(sTqLogT,nsc);
  }
  if(scArray) free(scArray);

  if(flagML) fclose(mlfile);
  /* if(flagWT && wtcolor) free(wtcolor);*/

}


/* ===================================================================== */
void checkWTMMparameters(int dimx) {
  /* =====================================================================
   *
   * Put the default parameters of WTMM analysis scheme if necessary
   * 
   * ===================================================================== */
  
  int flagratio; // useless temporary flag

  if(wav_wtmm == CRAZY) wav_wtmm = WAV_WTMM;
  if(ord_der_wtmm == CRAZY) ord_der_wtmm = ORD_DER_WTMM;

  if(ScTime == CRAZY) ScTime = SCTIME;

  if(ScRatio == CRAZY) ScRatio = RATIO;
  
  if(MinScale == CRAZY) MinScale = MINSCALE;
  
  if (MaxScale == CRAZY || MinScale>=MaxScale) {
    MaxScale = ScRatio*(((float)dimx)-1.0)/(float)ScTime;
    flagratio=YES;
  } else 
    flagratio=NO;

  if(ShiftSpectrum == CRAZY) ShiftSpectrum=SHIFTSPEC;

  if(NVoices == CRAZY)    NVoices = NVOICES;
  ScStep = pow(2.,(1./(double)NVoices));

  if (MinMoment >= MaxMoment) {
    MinMoment = MINMOMENT;
    MaxMoment = MAXMOMENT;
  }

  if(QStep == CRAZY) QStep = DQ;

  if(SCmin >= SCmax) {
    SCmin = SCMIN; 
    SCmax = SCMAX;
  }
 
  if(NSERIES > 1) flagML = NO;
  /* no representation at all: it's too huge and unsignificant */

  if(VERBOSE) { 
    printf("\n Wavelet to be used in WTMM scheme:   wav_wtmm=%d => %s",
	   wav_wtmm,wav_name[wav_wtmm+1]);
    if(wav_wtmm>0) 
      printf("\n Derivative order of the gaussian wavelet:"
	     "      deriv_order=%d ",ord_der_wtmm);
    printf("\n Number of voices: no_voices=%d", NVoices);
    printf("\n Range of scales of analysis in WTMM scheme:"
	   "\n      [scale_0=%g, scale_max=%g]", MinScale, MaxScale);
    if(flagratio == YES)  
      printf("\n Factor of convolution range in WTMM scheme:"
	     "      time_scale=%d" 
	     "\n Ratio between maximum scale and signal length:"
	     "      ratio=%g", ScTime, ScRatio);
    if(QStep != 0.)
      printf("\n Range of moments and moment step used to perform WTMM estimation:"
	     "\n      [min_moment=%g, max_moment=%g] - moment_step=%g",
	     MinMoment, MaxMoment, QStep);
    else {
      int i;
      printf("\n List of moments used to perform WTMM estimation:\n [");
      for( i=0; i<nqDef; i++ ) printf(" %g",qDefArray[i]);
      printf(" ]");
    }
    printf("\n Range of (log)scales used to run regression for WTMM estimation:"
	   "\n     [log(sc_min)=%g, log(sc_max)=%g]", LOG(SCmin), LOG(SCmax));
  }
  
}


/* ===================================================================== */
int checkdimension( char *base, int dimx ) {
/* ===================================================================== 
 *
 * Check the dimension of the input file
 *
 * ===================================================================== */
  char nombre[90];

  if(dimx <= 0) {
    if(NSERIES == 1 && IsBinary != BINARY) sprintf(nombre,"%s",base);
    else             sprintf(nombre,"%s-N000",base);
    if(IsBinary == BINARY) dimx=lee_dimension_datos( nombre );
    else                   dimx=lee_dimension_serie( nombre );

    if(VERBOSE) printf("\n Supposed size of the original signal: %d",dimx);
  }

  return dimx;
}
