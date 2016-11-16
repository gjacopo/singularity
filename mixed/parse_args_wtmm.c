/********************************************************/
/*                                                      */
/*  parse_args_wtmm.c. - Version tou 2 Dekembriou 2004  */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Personnal libraries */
#include <const_wtmm.h>
#include <utils.h>
#include <parse_args_wtmm.h>

/* Global variables passed to the functions */

/* parsing all other parameters of the analysis
 * See file const-msm.h for description  */
extern int VERBOSE;
extern int NSERIES;
extern int D_space;

extern int WTMMMethod;
extern int IsBinary;

extern int wav_wtmm;
extern int ord_der_wtmm;

extern float MinScale;
extern float MaxScale;
extern int   NVoices;
extern int   ScTime; 
extern float ScRatio; 

extern float MinMoment;
extern float MaxMoment;
extern float QStep;

extern float SCmin, SCmax;

extern float ShiftSpectrum;

extern int flagTauq;

#ifndef FLAG_MSM 
extern char infilename[];
extern char **ptrvar_c;
#endif


/* ===================================================================== */
int parsing_wtmm( int in0, int siflag, int *deflag, 
		  char **olarg, char **olval, char **olexp,
		  float **ptrvar_f, float **ptrval_f, 
		  int **ptrvar_i, int **ptrval_i,
		  int **ptrflag, int *type, int*olnumb)  {
/* ===================================================================== */
  int in;

  in=in0;

  //    Argument WTMMMethod. Type 0: flag

  sprintf(olarg[in],"%s","-canon"); 
  sprintf(olexp[in],"%s\n %s : %s\n",
	  "\nWTMM PARAMETERS\n===============",
	  olarg[in],
	  "Flag. If enabled, the method used in WTMM scheme to approximate the"
	  "\n  Legendre Transform is the canonical approach."
	  " Default:  DISABLED (direct method)"); 
  type[in]=0;
  olnumb[in]=1;
  ptrflag[in]=&WTMMMethod;
  in++;
 
  //   Argument wav_wtmm. Type 1: Integer

  strcpy(olarg[in],"-wav");
  DEFFLAG(olarg[in]); /* Add the string "wtmm" to the flag when
		       * compiling with the MSM library to avoid confusion
		       * with the already existing flag -wav */
  strcpy(olval[in],"wav_index");
  sprintf(olexp[in]," %s : %s %d\n",
	  olarg[in],
	  "Wavelet to be used in WTMM scheme, either:"
	  "\n      -1: Lorentzian,"
	  "\n       0: Morlet,"
	  "\n       1: Gaussian."
	  "\n  Default: ", WAV_WTMM);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=1;
  olnumb[in]=1;
  ptrvar_i[in]=&wav_wtmm;
  ptrval_i[in][0]=-1;
  ptrval_i[in][1]=1;
  in++;

  //   Argument ord_der_wtmm. Type 1: Integer

  strcpy(olarg[in],"-der");
  DEFFLAG(olarg[in]);
  strcpy(olval[in],"deriv_order");
  sprintf(olexp[in]," %s : %s %d\n",olarg[in],
	  "Derivative order, only when the wavelet used in WTMM scheme is the gaussian."
	  "\n  Default: ",
	  ORD_DER_WTMM);
  if(siflag)
    {
      type[in]=4;
      ptrflag[in]=deflag;
    }
  else type[in]=1;
  olnumb[in]=1;
  ptrvar_i[in]=&ord_der_wtmm;
  ptrval_i[in][0]=0;
  ptrval_i[in][1]=5;
  in++;

  //   Argument MinScale. Type 2: Float

  strcpy(olarg[in],"-sc0");
  DEFFLAG(olarg[in]);
  strcpy(olval[in],"scale_0");
  sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	  "Initial scale of analysis in WTMM scheme. Default: ",MINSCALE);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=1;
  ptrvar_f[in]=&MinScale;
  ptrval_f[in][0]=.33; 
  ptrval_f[in][1]=100.;
  in++;

  //   Argument MaxScale. Type 2: Float

  strcpy(olarg[in],"-scmax");
  DEFFLAG(olarg[in]);
  strcpy(olval[in],"scale_max");
  sprintf(olexp[in]," %s : %s\n",olarg[in],
	  "Maximal scale of analysis in WTMM scheme. By default, if it is not given,"
	  "\n  it is computed automatically according to the lenght of the signal (see options"
	  "\n  '-tsc' and '-scrat' below).");
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=1;
  ptrvar_f[in]=&MaxScale;
  ptrval_f[in][0]=1; 
  ptrval_f[in][1]=1000.;
  in++;

  //   Argument NVoices. Type 1: Integer

  strcpy(olarg[in],"-nvoi");
  //  DEFFLAG(olarg[in]);
  strcpy(olval[in],"no_voices");
  sprintf(olexp[in]," %s : %s %d\n",olarg[in],
	  "Number of voices per octave. The (multiplicative) scale step in WTMM scheme"
	  "\n  will be 2^<1/no_voices>. Default: ", NVOICES);
  if(siflag)
    {
      type[in]=4;
      ptrflag[in]=deflag;
    }
  else type[in]=1;
  olnumb[in]=1;
 ptrvar_i[in]=&NVoices;
  ptrval_i[in][0]=1;  
  ptrval_i[in][1]=100.;
  in++;

  //   Argument ScTime. Type 1: Integer

  strcpy(olarg[in],"-tsc");
  strcpy(olval[in],"time_scale");
  sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	  "Factor of convolution range in WTMM scheme, i.e.:\n"
	  "        [wavelet box] = #{-time_scale*nsc,...,time_scale*nsc}."
	  "\n  Default: ", SCTIME);
  if(siflag)
    {
      type[in]=4;
      ptrflag[in]=deflag;
    }
  else type[in]=1;
  olnumb[in]=1;
  ptrvar_i[in]=&ScTime;
  ptrval_i[in][0]=1;  
  ptrval_i[in][1]=100.;
  in++;

  //   Argument ScRatio. Type 2: Float

  strcpy(olarg[in],"-rsc");
  strcpy(olval[in],"ratio");
  sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	  "Ratio between maximum scale and signal length. Default: ",
	  RATIO);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=1;
  ptrvar_f[in]=&ScRatio;
  ptrval_f[in][0]=.001;
  ptrval_f[in][1]=10.;
  in++;

  //   Arguments MinMoment & MaxMoment. Type 2: Float

  strcpy(olarg[in],"-q");
  DEFFLAG(olarg[in]);
  strcpy(olval[in],"min_moment max_moment");
  sprintf(olexp[in]," %s : %s [%g,%g]\n",olarg[in],
	  "Range of moments used to perform the estimation of the multifractal exponents"
	  "\n  in WTMM scheme. Default: ",
	  MINMOMENT, MAXMOMENT);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=2;
  ptrvar_f[in]=&MinMoment;
  ptrvar_f[in+1]=&MaxMoment;
  ptrval_f[in][0]=ptrval_f[in+1][0]=-100.;
  ptrval_f[in][1]=ptrval_f[in+1][1]=100.;
  in+=2;


  //   Argument QStep. Type 2: Float

  strcpy(olarg[in],"-dq");
  DEFFLAG(olarg[in]);
  strcpy(olval[in],"moment_step");
  sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	  "Moment step of the WTMM analysis. If moment_step is choosed as 0, then a list"
	  "\n  of (irregularly spaced) moments will be used instead of the range of moments"
	  "\n  [min_moment, max_moment] (see variable qDefArray in file variables_wtmm.h)."
	  "\n  Default: ",DQ);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=1;
  ptrvar_f[in]=&QStep;
  ptrval_f[in][0]=0.;
  ptrval_f[in][1]=10.;
  in++;

  //   Argument SCmin & SCmax. Type 2: Float

  strcpy(olarg[in],"-regsc");
  //  DEFFLAG(olarg[in]);
  strcpy(olval[in],"sc_min sc_max");
  sprintf(olexp[in]," %s : %s [%g,%g]\n",olarg[in],
	  "Range of scales used to perform the estimation (through regression) of the"
	  "\n  multifractal exponents in WTMM scheme. Default: ",
	  SCMIN, SCMAX);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=2;
  ptrvar_f[in]=&SCmin;
  ptrvar_f[in+1]=&SCmax;
  ptrval_f[in][0]=ptrval_f[in+1][0]=0.1;
  ptrval_f[in][1]=ptrval_f[in+1][1]=1000.;
  in+=2;

  //   Argument ShiftSpectrum. Type 2: Float

  strcpy(olarg[in],"-cdh");
  strcpy(olval[in],"shift_spectrum");
  sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	  "Constant to add to the values of the WTMM spectrum when storing it."
	  "\n  Default: ",SHIFTSPEC);
  if(siflag)
    {
      type[in]=5;
      ptrflag[in]=deflag;
    }
  else type[in]=2;
  olnumb[in]=1;
  ptrvar_f[in]=&ShiftSpectrum;
  ptrval_f[in][0]=-2.;
  ptrval_f[in][1]=2.;
  in++;

  //    Argument flagTauq. Type 0: flag

  sprintf(olarg[in],"%s","-tauq"); 
  sprintf(olexp[in]," %s : %s\n",
	  olarg[in],
	  "Flag. If enabled, the multifractal exponent computed by the WTMM scheme"
	  "\n  will be saved. Default:  DISABLED"); 
  type[in]=0;
  olnumb[in]=1;
  ptrflag[in]=&flagTauq;
  in++;
 
  return(in);
}


#ifndef FLAG_MSM

/* ===================================================================== */
void parse_arguments(int argc, char *argv[]) {
/* ===================================================================== */

  char **olarg,**olval,**olexp;
  float **ptrvar_f,**ptrval_f;
  int **ptrvar_i,**ptrval_i;
  int **ptrflag;
  int *type, *olnumb;
  int lar;
  int flagv;
  int arglen;
  int Narg0=30;  // Initialization value; it should be greater than 
                 // (but not necessarily equal to) the number or arguments
  int in,Narg;
  int i, cur;

  /*         Defining the expected arguments            */
	
  olarg=reservar_matriz_char(Narg0,20);   // Name of the variable 
  olexp=reservar_matriz_char(Narg0,500);  // Explanation
  ptrvar_f=(float **) calloc(Narg0,sizeof(float *)); 
  // Pointer to the variable, if it is double
  // ptrval_f=reservar_matriz(Narg0,2); // Allowed min and max
  ptrval_f=reservar_matriz_float(Narg0,2); 
  ptrvar_i=(int **) calloc(Narg0,sizeof(int *)); // Pointer to the variable, if int
  ptrval_i=reservar_matriz_int(Narg0,2); // Allowed min and max
  ptrflag=(int **) calloc(Narg0,sizeof(int *)); // Pointer to the flag
  type=(int *) calloc(Narg0,sizeof(int)); // Type of variable
  /* Personal modification */ 
  olnumb=(int *) calloc(Narg0,sizeof(int)); // Number of arguments to follow 
  olval=reservar_matriz_char(Narg0,30);   // Name of the argument

  //  Pointer to the variable, if char
  ptrvar_c=( char**) calloc(Narg0,sizeof(char*));

  in=0;
	
  //    Argument NAME. Type 0: char list

  sprintf(olarg[in],"%s","-f");
  sprintf(olval[in],"%s","file_name");
  sprintf(olexp[in],"%s\n %s : %s\n",
	  "\nGENERAL INPUT PARAMETERS\n========================",
	  olarg[in],
	  "The name of the input file to be processed. It's a generic (base) name if the"
	  "\n  number NSERIES (option '-N') is >1 or if it's a binary file. It will be completed"
	  "\n  with the 'usual' strings '-N???', and also '.dat' or '.txt' according to the type"
	  "\n  of the file: see option '-txt' for flag BINARY. Otherwise, if it represents a"
	  "\n  unique text file to be processed, the string won't be completed and its size is"
	  "\n  automatically computed."); 
  type[in]=-1;
  olnumb[in]=1;
  ptrvar_c[in]=&(infilename[0]);
  in++;

  //    Argument NSERIES. Type 1: integer

  sprintf(olarg[in],"%s","-N");
  sprintf(olval[in],"%s","#series");
  sprintf(olexp[in]," %s : %s %d\n",
	  olarg[in],
	  "Number of series to be processed. Default:",NSERIES); 
  type[in]=1;
  olnumb[in]=1;
  ptrvar_i[in]=&NSERIES;
  ptrval_i[in][0]=1;
  ptrval_i[in][1]=10000;
  in++;

  //    Argument D_space. Type 1: integer

  /* NOT YET IMPLEMENTED FOR D_SPACE>1
     sprintf(olarg[in],"%s","-d_space");
     sprintf(olval[in],"%s","dimension");
     sprintf(olexp[in]," %s : %s %d\n",
     olarg[in],
     "Dimension of the embedding space. Default:",
     D_space); 
     type[in]=1;
     olnumb[in]=1;
     ptrvar_i[in]=&D_space;
     ptrval_i[in][0]=1;
     ptrval_i[in][1]=2;
     in++;
  */

  //    Argument IsBinary. Type 0: flag

  sprintf(olarg[in],"%s","-txt");
  sprintf(olexp[in]," %s : %s\n",
	  olarg[in],
	  "Flag. If enabled, the input data are supposed to be in text format,"
	  "\n  otherwise they are in binary format. Default: DISABLED"); 
  type[in]=0;
  olnumb[in]=1;
  ptrflag[in]=&IsBinary;
  in++;

 /*  multifractal parameters for WTMM analysis */

  in=parsing_wtmm(in,0,0,olarg,olval,olexp,ptrvar_f,ptrval_f,
			     ptrvar_i,ptrval_i,ptrflag,type, olnumb);

  //    Argument VERBOSE. Type 0: flag

  sprintf(olarg[in],"%s","-ver");
  sprintf(olexp[in],"%s\n %s : %s\n",
	  "\nGENERAL PARAMETER\n==================",
	  olarg[in],
	  "Flag. If enabled, the program shows lots of (verbose) information,"
	  "\n  specially for multifractal analysis. Default: DISABLED"); 
  type[in]=0;
  olnumb[in]=1;
  ptrflag[in]=&VERBOSE;
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
	  /**/
	  for ( i=0; i<olnumb[in]; i++ ) {
	  /**/
	    if (type[in] == -1) {
	      lar++;
	      strcpy(ptrvar_c[in+i],argv[lar]);
	      
	    }  else if(type[in]%3==2)
	      {
		lar++;
		sscanf(argv[lar],"%f",ptrvar_f[in+i]);
		if((*ptrvar_f[in+i]<ptrval_f[in+i][0])
		   ||(*ptrvar_f[in+i]>ptrval_f[in+i][1]))
		  {
		    printf("\nValue out of range (%g - %g) for argument %s\n",
			   ptrval_f[in+i][0],ptrval_f[in+i][1],olarg[in]);
		    flagv=0;
		  }
		if(type[in]==5) *ptrflag[in+i]=YES;
	      }
	    else if(type[in]%3==1)
	      {
		lar++;
		sscanf(argv[lar],"%d",ptrvar_i[in+i]);
		if((*ptrvar_i[in+i]<ptrval_i[in+i][0])
		   ||(*ptrvar_i[in+i]>ptrval_i[in+i][1]))
		  {
		    printf("\nValue out of range (%d - %d) for argument %s\n\n",
			   ptrval_i[in+i][0],ptrval_i[in+i][1],olarg[in]);
		    flagv=0;
		  }
		if(type[in]==4) *ptrflag[in+i]=YES;
	      }
	    else
	      {
		*ptrflag[in+i]=YES;
		if(type[in]==3) *ptrvar_i[in+i]=1;
	      }
	  }
	  
	  /**/
	}
      /**/
      
      /*            End of parsing loop             */
    }	
  
  if(argc==1 || flagv!=1)
    {
      printf("Usage: %s",argv[0]);
      for(in=0;in<Narg;in++) {
	if(type[in]%3==0) printf(" [%s]",olarg[in]);
	else printf(" [%s %s]",olarg[in],olval[in]);
	for( i=0,cur=in; i<olnumb[cur]-1; i++ ) 
	  /* skip some indices */ in++; 
      }
      printf("\n");
    }

  if(flagv==2)
    for(in=0;in<Narg;in++) printf("%s",olexp[in]);
  
  /*            Freeing memory before terminating        */
  
  liberar_matriz_char(olarg,Narg0);
  liberar_matriz_char(olval,Narg0);
  free(ptrvar_f);
  liberar_matriz_float(ptrval_f,Narg0);
  // liberar_matriz(ptrval_f,Narg0);
  free(ptrvar_i);
  liberar_matriz_int(ptrval_i,Narg0);
  free(ptrflag);
  free(type);
  free(olnumb);
	
  free(ptrvar_c);

  /*                Termination              */

  if(argc==1 || flagv!=1) exit(-1);

}

#endif /* #ifndef FLAG_MSM */
