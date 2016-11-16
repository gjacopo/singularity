/********************************************************/
/*                                                      */
/*    MF_processor.c - Version: 18 de Octubre, 2004     */
/*                                                      */
/* Main program for launching the multifractal process  */
/*   - multifractal signals generation,                 */
/*   - multifractal signals analysis with Histogram,    */
/*     Singularity analyzis and WTMM schemes.           */
/*                                                      */
/********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Personnal MSM libraries  */
#include <struct_def.h>
#include <tensor.h>
#include <parse_args.h>
#include <graficos.h>
#include <operaciones.h>
#include <FFT.h>
#include <derivacion.h>
#include <multifractal.h>
#include <multifractal_generator.h>
#include <multifractal_analysis.h>

/* Personal WTMM libraries */
#include <const_wtmm.h>
#include <parse_args_wtmm.h>
#include <extrema.h>
#include <approx_legendre.h>
#include <multifractal_wtmm.h>

/* all the global variables are defined in these headers */
#include <variables_msm.h>
#include <variables_wtmm.h>


int main(int argc, char *argv[])
{

  /*	PARAMETERS		*/
  char base[90];
  int leff,dimy;

  /*		Program		*/

  parse_arguments(argc,argv);

  /*             Genera multifractal     */

  leff = crea_genera_multifractal( base );
  dimy=(D_space==1)?1:leff;

  /*       Analyzing series with multifractal methods      */

  if(ANALYSIS) {
    analiza_series( leff, dimy, base);
    if( dimy == 1)    /* analiza_series_WTMM works for dimy=1 only */
      analiza_series_WTMM( leff, dimy, base );
  }

  return 0;
}
	
