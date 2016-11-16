/********************************************************/
/*                                                      */
/*                  wtmm_processor.c                    */
/*           Version del 2 de Dekembriou, 2004          */
/*                                                      */
/* Main stand-alone program for launching the WTMM      */
/* analyzis scheme for multifractal estimations.        */
/*                                                      */
/********************************************************/


#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Personal libraries */
#include <const_wtmm.h>
#include <parse_args_wtmm.h>
#include <extrema.h>
#include <approx_legendre.h>
#include <multifractal_wtmm.h>

/* all the global variables are defined in this header */
#include <variables_wtmm.h>


 /* Global variables used by the stand alone program only */
char infilename[90];
char **ptrvar_c;

/* ===================================================================== */
int main(int argc, char *argv[]) {
/* ===================================================================== */

  int dimx=0; /* Let a non positive number here as dimx initialization */
  int dimy;

  parse_arguments(argc, argv);
  dimy=(D_space==1)?1:dimx;
  
  analiza_series_WTMM( dimx, dimy, infilename);
  

  return 0;
}

