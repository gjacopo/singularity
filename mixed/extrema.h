#ifndef   	EXTREMA_H
# define   	EXTREMA_H

void WTMM_PF_compute( double *signal, int dimx, int nn, 
		     double expon, int order, 
		     int nq, double *qArray, int nsc, double *scArray, 
		     double **ExtWTlis, double **Z, int *n_ext );

void WT_transform( double *signal, int dimx, double expon, int order,
		   double sc,  double *wtrans );

void define_wavelet( int dimx, double sc, double expon, int order, 
			double *wavesig);

void simple_convolve( double *signal, int dimx, double *filt, int sfilt, 
		      double sc, double *wtrans );

double comp_coefficient(int n, double t);

int WTMM_find_extrema( double *wtrans, int *wtrans_ind, 
			int dimx, double sc, int tempi_max );
  
int PF_extrema_track( double *wtrans, int *wtrans_ind, 
		       double *maxsig, int *maxsig_ind,
		       int jj, int  n_max,
		       int nq, double *qArray, int nsc, double *scArray, 
		       double *ExtWTlis, double *Z );




#endif 	    /* !EXTREMA_H */
