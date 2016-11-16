#ifndef   	MULTIFRACTAL_WTMM_H
# define   	MULTIFRACTAL_WTMM_H

void analiza_series_WTMM( int leff, int dimy, char *base);

void estima_Dh_WTMM( int leff, int dimy, char *base,
			int nq, double *qArray, 
			double *tauq, double*h, double*Dh );

int checkdimension( char *base, int dimx );

void checkWTMMparameters(int dimx);

#endif 	    /* !MULTIFRACTAL_WTMM_H */
