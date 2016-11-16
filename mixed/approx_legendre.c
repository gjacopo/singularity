/********************************************************/
/*                                                      */
/*                 approx_legendre.c                    */
/*           Version del 2 de Dekembriou, 2004          */
/*                                                      */
/*  Implementation of the Legendre transform.           */
/*                                                      */
/* List of functions:                                   */
/* directSpecCompute: classical method                  */
/* canonMeanCompute, canonSpecCompute: canonical method */
/* (used in LastWave)                                   */
/* linefit                                              */
/* compare                                              */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

/* Personnal libraries */
#include <const_wtmm.h>
#include <approx_legendre.h>


extern int VERBOSE;

/**
 ** ===================================================================== **
 ** METHOD I : direct computation of the density spectrum through the 
 **            exponent tau(q) and the partition function
 ** ===================================================================== **
 **/

/* ===================================================================== */
void directSpecCompute( double **Z, int dimx, 
			int nq, double *qArray, int nws, double *wsArray,
			double logaMin, double logaMax,
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
   *    - logaMin, logaMax : borns of the range of (log)scales used in the
   *      estimation of the exponents.
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
    Lws = LOG( wsArray[ind_ws] );

    if(Lws>=logaMin && Lws<=logaMax) {
      /* number of scales considered till now */
      count ++;
      /* sum of scales */
      sumsc += Lws;
      /* sum of squared scales */
      sumsqsc += SQR(Lws);
      for( ind_q=0; ind_q<nq; ind_q++ ) {
	/* log of the partition function */
	LZ = LOG( Z[ind_ws][ind_q] );
	/* sum of the partition function */
	sumpf[ind_q] += LZ;
	/* sum of partition function weightened by the scale  */
	sumscpf[ind_q] += Lws * LZ;	
      }
    }
  } /* end loop over the scales: "for (ind_ws=0;..." */
  
   /** Compute the tauq */
  for(ind_q=0; ind_q<nq; ind_q++) {
    /** Compute the tauq exponent corresponding to qArray[ind_q]
     * */
    Tauq[ind_q] = (sumscpf[ind_q]*count - sumsc*sumpf[ind_q]) 
      / (count*sumsqsc - sumsc*sumsc);
    /* Note: this is probably nothing else than : 
     *        Tauq[ind_q] = sumscpf[ind_q] / sumsqsc;
     * to be checked...     
     */

    /* MODIF
       Tauq[ind_q] -= qArray[ind_q]/2.;
       END MODIF */


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
  for( i=0; i<nws; i++ )   maxsize = MAX( maxsize, n_ext[i] );
  
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
	tempLogT[i] = LOG(tempWT[i] / tm);
      }
      
      /** Note that the sum below allow to compute statistical variables 
       ** over several signals (NSERIES>1) by simply adding their values **/

      /* We compute sTq and sTqLogT */
      if(q >= 0.)      
	for( i=imin; i<size; i++ )	{
	  sTq[ind_ws][ind_q] += tempTq[i]; /* pow( tempWT[i], q ); */
	  sTqLogT[ind_ws][ind_q] += tempTq[i] * tempLogT[i];
	  /* pow( tempWT[i], q ) * LOG( tempWT[i]/tm ); */
	}
      
      else 
	for( i=size-1; i>=imin; i-- )	{
	  sTq[ind_ws][ind_q] += tempTq[i];/* pow( tempWT[i], q ); */
	  sTqLogT[ind_ws][ind_q] += tempTq[i] * tempLogT[i]; 
	  /* pow( tempWT[i], q ) * LOG( tempWT[i]/tm ); */
	}
      
      /* sTqLogT = sTqLogT + log(tm)*sTq */
      sTqLogT[ind_ws][ind_q] += LOG(tm)*sTq[ind_ws][ind_q];
      
      /* We compute LogSTq
	 if(sTq[ind_ws][ind_q] != 0.0)	{
	 logSTq[ind_ws][ind_q] = 
	 LOG(sTq[ind_ws][ind_q]/((double) (size-imin)));
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
		       double logaMin, double logaMax,
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
   *    - nws, wsArray : resp. no of sceles and list of scales,
   *    - logaMin, logaMax : borns of the range of (log)scales used for the
   *      estimation of the multifractal exponents.
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
    LwsArray[ind_ws] = LOG( wsArray[ind_ws] );
  
  for(ind_q=0; ind_q<nq; ind_q++) { 
    
    /* consider the current moment for computing the variables */
    q=qArray[ind_q];
    
    for( ind_ws=0; ind_ws<nws; ind_ws++ ) { 
      stq = sTq[ind_ws][ind_q];
      
      /** Compute the canonical multifractal exponents tq */ 
      if(stq == 0.)	tq[ind_ws] = 0.;
      else   tq[ind_ws] = (double) LOG(stq);
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
				   - LOG(stq));
      /* equivalent to the function PFAccessDQFloat */
      
    } /* end loop over the scales: "for (ind_ws=0;..." */
    
    
    /* Regression over the scales to get the different exponents.
     * For each moment, the tq, hq and dq are tabulars of nws values
     * whose slopes give the corresponding exponent Tauq, H and Dh. */
    ss = linefit( tq, LwsArray, nws, logaMin, logaMax, &(Tauq[ind_q]) );
    
    linefit( hq, LwsArray, nws, logaMin, logaMax, &(H[ind_q]) );
    
    linefit( dq, LwsArray, nws, logaMin, logaMax, &(Dh[ind_q]) );
    
  } /* end loop over the moments: "for (ind_q=0;..." */
  
  if(VERBOSE) 
    printf("\n     %d scales used for the estimation of multifractal exponents", ss);
  
  /* Free memory */ 
  if(tq) free(tq);
  if(LwsArray) free(LwsArray);
}



/* ===================================================================== */
int
linefit( double *yValues, double *xValues, int dimx, 
	      double xMin, double xMax, double *slope )  {
  /* ===================================================================== 
   *
   * Fit a signal with a straight line by regression, i.e. it computes the
   * slope a of the line y=ax+b that best fits the data [xValues,yValues].
   *
   * Parameters: 
   *      - xValues : the list of indexes of the signal,
   *      - yValues : the values of the signal to fit, 
   *      - dimx : lenght of the signal,
   *      - xMin, xMax : borns of the range of values of xValues where
   *        the signal will be fitted,
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
    if ( xValues[i]>=xMin && xValues[i]<=xMax ) {
      sx += xValues[i];
      sy += yValues[i];
      ss++;
    }

  sxoss = sx/(double)ss;

  for( i=0; i<dimx; i++ ) 
    if ( xValues[i]>=xMin && xValues[i]<=xMax ) {
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


