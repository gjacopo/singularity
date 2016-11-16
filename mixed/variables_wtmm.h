/********************************************************/
/*                                                      */
/*                  variables-wtmm.h                    */
/*           Version del 2 de Dekembriou, 2004          */
/*                                                      */
/*  Global variables used in the WTMM analyzis scheme.  */
/*                                                      */
/********************************************************/

#ifndef   	VARIABLES_WTMM_H
# define   	VARIABLES_WTMM_H

/** Define the FLAG_WTMM variable **/
#ifndef FLAG_WTMM
#define FLAG_WTMM
#endif

/* dummy variable */
#ifndef CRAZY
#define CRAZY -10.
#endif

int IsBinary=BINARY;

int WTMMMethod=DIRECT; /* Choice of the method used for Legendre transform 
			   * approximation: direct or canonical */

int wav_wtmm=-10;     /* Choice of the analyzing wavelet. 
			      * Default will be Gaussian (1) */
char wav_name[][20]={"lorentzian","morlet","gaussian"};

int ord_der_wtmm=-10; /* Order of the gaussian derivative used for wavelet
			      * analysis */

float MinScale=CRAZY; /* Minimum scale of analysis */
float MaxScale=CRAZY; /* Maximum scale of analysis */
int NVoices=CRAZY;  /* Number of voices per octave. The scale resolution ScStep
		     * for which the partition function is calculated is 
		     * 2^(1/NVoices) */
float ScStep;

int ScTime=-10; /* Factor of convolution range; i.e.: 
		        *    [wavelet box] = #{-ScTime*nsc,...,ScTime*nsc} */
float ScRatio=CRAZY;    /* Ratio between maximum scale (wavelet box) and 
		        * signal length */


float MinMoment=CRAZY; /* Minimum moment for estimation */ 
float MaxMoment=CRAZY; /* Maximum moment for estimation */ 
float QStep=CRAZY;     /* Resolution of moments q */

/* Variable used in this file only: default list of moments q, 
 * where they are not regularly spaced */
float qDefArray[]={-4.,-3.6,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.3,2.6,3.,3.5,4.,5.,6.,7.,8.};
int nqDef=65;

float SCmin=CRAZY; /* Log range of scales to be used in the estimation */
float SCmax=CRAZY; /* of multifractal exponents through regression */

float ShiftSpectrum=CRAZY; /* Constant to add to the values of the spectrum when
			    * representing it - Useful if we want to compare
			    * different spectra */ 

int flagML=NO; 
int flagTauq=NO;
int flagPF=NO;

double *wtcolor;
FILE *mlfile;

#ifndef FLAG_MSM 
/* Constant of the MSM lib used in some functions of the WTMM lib
 * If they are not declared elsewhere, wee do it here */
int D_space=1;
int NSERIES=1;
int VERBOSE=YES;
#endif



#endif 	    /* !VARIABLES_WTMM_H */
