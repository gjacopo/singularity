#ifndef   	CONST_WTMM_H_
# define   	CONST_WTMM_H_

/* Define the FLAG_WTMM variable */
#ifndef FLAG_WTMM
#define FLAG_WTMM
#endif

/* Define the method used for estimating the Legendre transform */
#define  DIRECT 0
#define  CANONICAL (1-DIRECT)

/** Some useful routines **/
#define SQR(x) ( (x) * (x) )
#define SQRT(x) ( sqrt(x) )
#define MIN(x,y) ( ((x) < (y)) ? (x) : (y) )
#define MAX(x,y) ( ((x) > (y)) ? (x) : (y) )
#define SIGN(x) (((x)>0.)?(1):(-1))

#define BASELOG ((double)10.)
/* #define BASELOG ((double)2.) */
#define LOG(x) (log((x))/log(BASELOG))
/* #define LOG(x) (log(x)) */

/** Some default pseudo constant **/

#define WAV_WTMM 1      /* Default choice for wavelet: gaussian */
#define ORD_DER_WTMM 2  

#define MINSCALE 2.0   
//  #define MAXSCALE 
#define SCTIME 5.       
#define RATIO (1./16.) 
#define NVOICES 10

#define QDEFLIST {-4.,-3.6,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.3,2.6,3.,3.5,4.,5.,6.,7.,8.}
#define NQDEF 65

#define MINMOMENT -5.
#define MAXMOMENT 5.
#define DQ 0.2         

#define SCMIN 2.
#define SCMAX 100.

#define SHIFTSPEC 0.

#define   NG 255 

/* Variable used for color representation of WT */
#define WTFILE   "wtcolor.ppm"
#define MLFILE   "maxlines.txt"
#define PFFILE   "partfunc.txt"
#define TAUQFILE "tauq.txt"
#define SPECFILE "spectrum.txt"

#define DEBUG 0

#define BINARY 0

/* Useful routine to modify the name of some option to be parsed
 * to the program when they have the same name in both methods */
#ifdef FLAG_MSM 
#define STRWTMM "wtmm"
#else
#define STRWTMM "" /* empty string */
#endif
#define DEFFLAG(x) ( strcat((x),STRWTMM) )


/* dummy variable */
#define CRAZY -10.

#endif 	    /* !CONST_WTMM_H */
