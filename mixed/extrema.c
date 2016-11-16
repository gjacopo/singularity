/********************************************************/
/*                                                      */
/*    extrema.c - Version del 2 de Dekembriou, 2004     */
/*                                                      */
/*  Implementation of the Extrema Extraction from the   */ 
/*  Continuous Wavelet Transform (CWT).                 */ 
/*                                                      */
/* List of functions:                                   */
/* WTMM_PF_compute: computes the CWT and                */
/*  - store the WT computed on selected extrema (the    */
/*    so-called WTMM) for further use,                  */
/*  - computes the factors of partition function from   */
/*    these WTMM.                                       */
/* WTMM_find_extrema                                    */
/* PF_extrema_track                                     */
/* WT_transform: computes the wavelet used in CWT.      */
/* simple_convolve                                      */
/* define_wavelet                                       */
/* comp_coefficient                                     */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

/* Personnal libraries */
#include <const_wtmm.h>
#include <extrema.h>

/* Global external variables to be defined in main program */

extern int ScTime;

extern int flagML; 
extern FILE *mlfile;


/* ===================================================================== */
void WTMM_PF_compute( double *signal, int dimx, int maxfilt, 
		     double wav, int order,
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
  int icolor=1;
  int *wtrans_ind, *maxsig_ind, *ptl;
  int ind_sc, ind_q, n_max;

   /* Allocations of memory
    * Note: making the allocations here (once for all, instead of making
    * it for each scale) enables to reduce time consuming and to make
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
  
  for( ind_sc=0; ind_sc<nsc; ind_sc++ ) {
    sc=scArray[ind_sc]; /* current scale of analysis */
    
    /** Wavelet transform of the signal for increasing wavelet scale. 
     **/
    /* Direct computation (in temporal space, not in Fourier space)
     * of the WT at scale sc */
    WT_transform( signal, dimx, wav, order, sc, /* maxfilt, */ wtrans );
        
    /** Find the local maxima of the wavelet coefficient for each scale **/
    n_max = WTMM_find_extrema( wtrans, wtrans_ind, dimx, sc, maxfilt ); 
    /* n_max : number of extrema at current scale */
    
    /** Facultative representation of the Maxima Lines **/
    if(flagML && mlfile!=NULL)
      for( i=0; i<n_max; i++ ) 
	fprintf(mlfile,"%g %d\n", LOG(sc), wtrans_ind[i]+1);
    /* the indexes of the selected extrema at each scale are stored */
    
    /** Tracking the maxima lines: test for supremum 
     ** Generate maxima lines and compute the partition function. */
    if( ind_sc > 0)    /* i.e.: if sc > min_sc */
      n_ext[ind_sc] =
	PF_extrema_track( wtrans, wtrans_ind, maxsig, maxsig_ind, 
			  jj, n_max, nq, qArray, nsc, scArray, 
			  ExtWTlis[ind_sc], Z[ind_sc] );
    
    /* DEBUG    {
      FILE *f;
      if(!ind_sc) f=fopen("update-wtmm.txt","w"); 
      else        f=fopen("update-wtmm.txt","a");
      for ( i=0; i<n_ext[ind_sc]; i++ ) fprintf(f,"%g ", ExtWTlis[ind_sc][i]);
      fprintf(f,"\n");
      fclose(f);
      } */

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
			double sc, int maxfilt ) {
  /* =====================================================================
   *
   * Find the local maxima of the wavelet coefficient for each scale.
   *
   * Parameters:
   *    - sc : current scale of analysis,
   *    - wtrans : 1d tabular storing  the values of the WT at scale sc,
   *    - maxfilt : maximal size of the wavelet filter,
   *    - dimx : common size of both tabulars wtrans and wtrans_ind.
   * Outputs:
   *    - wtrans : it is modified to finally store the values of the extrema
   *      only,
   *    - wtrans_ind : 1d tabular storing the indexes of the selected 
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
  
  sign = SIGN(wtrans[maxfilt+1] - wtrans[maxfilt]);
  temp1 = 0.; 
  for( j=0, i=maxfilt+2; i<(dimx-maxfilt-1); i++ ) {  
    if ((fabs(wtrans[i] - wtrans[i-1]) > 0.0) && 
	((sign == 1) && ((SIGN(wtrans[i] - wtrans[i-1])) == -1))) {
      temp = wtrans[i-1];
            
      n_max++;
      /* affectation du max */
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
   *    - wtrans, wtrans_ind : 1d tabulars storing resp. the values and the
   *      indexes of the WT at the current scale,
   *    - maxsig, maxsig_ind : 1d tabulars storing resp. the values and the
   *      indexes of the WT at the previous scale,
   *    - n_max : number of local extrema,
   *    - jj : index of the last detected extrema in the list of extremas
   *      already selected.
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nsc, scArray : resp. no of sceles and list of scales,
   * Outputs:
   *    - Z : 1d tabular (indexed by moment) storing the factors of the 
   *      partition function at current scale,
   *    - ExtWTlis : 1d tabular storing the WTMM over selected extrema
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
  while (((i1-1) < n_max) && ((i2-1) < jj)) {
    /* Choose the supremum between all previous values at smaller scales. */
    if ((wtrans_ind[i1] - maxsig_ind[i2]) <= 
	(maxsig_ind[i2+1] - wtrans_ind[i1]))
      wtrans[i1] = MAX( wtrans[i1], maxsig[i2] );
    else
      wtrans[i1] = MAX( wtrans[i1], maxsig[i2+1] );
    
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
   * Note : unefficient but safe implementation (no FFT).
   * Called by : WTMM_PF_compute  
   * ===================================================================== */
 
  int i, j;
  int tempi;
  double *wavesig;
  
  /* size of the wavelet filter at scale sc (support of the wavelet) */
  tempi = (int)((double)ScTime*sc);

  /* allocate memory for the wavelet signal */
  if( (wavesig = (double*)calloc(2*tempi+1,sizeof(double))) == NULL ) {
    fprintf(stderr,"\n Error allocation in WT_transform");
    exit(-1);
  }
  
  /* useless initializations:
     for (i=0; i<2*maxfilt+1; i++)  wavesig[i] = 0.;  
     for (i=0; i<dimx+1; i++) wtrans[i] = 0.;
  */

  /* Create the wavelet filter */
   define_wavelet( tempi, sc, wav, order, wavesig ); 

  /* Then do the wavelet convolution of the signal
   * !!! it would be faster with FFT !!! 
   * but i reduce numerical uncertainty */
   simple_convolve( signal, dimx, wavesig, tempi, sc, wtrans );

  if(wavesig)  free(wavesig);
}


/* ===================================================================== */
void simple_convolve( double *signal, int dimx, double *filt, int sfilt, 
		      double sc, double *wtrans ) {
  /* ===================================================================== 
   *
   * Naive convolution + normalization to compute the wavelet coefficient
   *
   * Parameters:
   *    - signal : original signal,
   *    - dimx : size of the signal,
   *    - sc : current scale of analysis,
   *    - filt : wavelet filter used to convolve the signal at scale sc,
   *    - sfilt : size of the wavelet filter.
   * Output:
   *    - wtrans : values of the WT at scale sc,
   *
   * ===================================================================== 
   * Called by : WT_transform
   * ===================================================================== */
  
  double conv;
  int i, j;
  
  for (i=sfilt; i<dimx-sfilt-1; i++) {
    for( j=i-sfilt, conv=0.; j<=i+sfilt; j++ )
      conv += signal[j] * filt[j-i+sfilt];
    /* absolute value of wavelet coefficient */
    wtrans[i] = fabs(conv/sc);  
  }
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
   *    - dimx : dimension of the output wavelet function (this will be
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
  
  for( ix=0; ix<=2*dimx; ix++ ){
    t=((double)ix-(double)dimx)/sc;
    /* implementation turiel
    t = (double)ix;
    if(ix>=dimx/2) t -= (double)dimx;
    t = t/sc;
    */
    if(wav<0.) wavesig[ix] = pow(1.+SQR(t),-1); /* lorentzian */
    else if (wav == 0.) wavesig[ix] = exp(-SQR(t)/2.) * cos(5.*t); /* morlet */
    else /*if(wav>0.)*/  wavesig[ix] = comp_coefficient(order,t); /* gaussian */
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
  case 0:  return exp(-0.5*SQR(t));
  case 1:  return -t * exp(-0.5*SQR(t));
  case 3:  return t * exp(-0.5*SQR(t)) * (3-SQR(t));
  case 4:  return exp(-0.5*SQR(t)) * (pow(t,4.0) - 6*SQR(t) + 3);
  case 5:  return -t * exp(-0.5*SQR(t))* (pow(t,4.0) - 10*SQR(t) + 15);
  default:
  case 2:  return (1-SQR(t)) * exp(-0.5*SQR(t));
    // case 2:  return (1-SQR(t)) * exp(-0.5*SQR(t)) * 2./ (SQRT(3)* pow(M_PI,0.25));
  }
}



/*************************************************************************/


/* ===================================================================== */
static int FFT( double *x, double *y, int m, int dir) {
  /* ===================================================================== 
   * This computes an in-place complex-to-complex FFT
   *
   * Parameters :
   *    - x, y : real and imaginary arrays of 2^m points resp.
   *    - dir : dir=1 gives forward transform, dir=-1 gives reverse transform  
   *
   * ===================================================================== 
   * Not used here
   * ===================================================================== */
  
   int nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;
   
   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)  nn *= 2;
   
   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
     if (i < j) {
       tx = x[i];       ty = y[i];
       x[i] = x[j];     y[i] = y[j];
       x[j] = tx;       y[j] = ty;
     }
     k = i2;
     while (k <= j) {
       j -= k;
       k >>= 1;
     }
     j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
     l1 = l2;
     l2 <<= 1;
     u1 = 1.0;
     u2 = 0.0;
     for (j=0;j<l1;j++) {
       for (i=j;i<nn;i+=l2) {
	 i1 = i + l1;
	 t1 = u1 * x[i1] - u2 * y[i1];
	 t2 = u1 * y[i1] + u2 * x[i1];
	 x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
       }
       z =  u1 * c1 - u2 * c2;
       u2 = u1 * c2 + u2 * c1;
       u1 = z;
     }
     c2 = sqrt((1.0 - c1) / 2.0);
     if (dir == 1)   c2 = -c2;
     c1 = sqrt((1.0 + c1) / 2.0);
   }
   
   /* Scaling for forward transform */
   if (dir == 1) {
     for (i=0;i<nn;i++) {
       x[i] /= (double)nn;
       y[i] /= (double)nn;
     }
   }

   return 0;
}


/* ===================================================================== */
static void convolve( int dimx, int nb2, double *R1, double *R2) {
  /* ===================================================================== 
   * Not used here
   * ===================================================================== */ 
  
  double *I1,*I2;
  double tmp;
  int ix;

  if( (I1=(double *)calloc(dimx,sizeof(double))) == NULL ||
      (I2=(double *)calloc(dimx,sizeof(double))) == NULL ) {
    fprintf(stderr,"\n Error allocation in convolve");
    exit(-1);
  }

  FFT( R2, I2, nb2, 1 );
  FFT( R1, I1, nb2, 1 );
  for( ix=0; ix<dimx; ix++ )
    {
      tmp = R1[ix]*R2[ix] - I1[ix]*I2[ix];
      I2[ix] = R1[ix]*I2[ix] + I1[ix]*R2[ix];
      R2[ix] = tmp;
    }
  FFT( R2, I2, nb2, -1 );

  if( I1 ) free( I1 );
  if( I2 ) free( I2 );
}

