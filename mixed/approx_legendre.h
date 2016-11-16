
#ifndef   	APPROXLEGENDRE_H_
# define   	APPROXLEGENDRE_H_



int //double
linefit( double *yValues, double *xValues, int N, 
	 double xMin, double xMax, double *slope );

int compare(const double *d1,const double *d2);

void directSpecCompute( double **Z, int N, 
			int nq, double *qArray, int nws, double *wsArray,
			double logaMin, double logaMax,
			double *tauq, double *h, double *Dh );

void canonMeanCompute( double **ExtWTlis, int N, int *n_ext,
		       int nq, double *qArray, int nws, double *wsArray,
		       double **sTq, double **sTqLogT /* double **logSTq */ );

void canonSpecCompute( double **sTq, double **sTqLogT, /* double **logSTq */
		       int N, int nq, double *qArray, int nws, double *wsArray,
		       double logaMin, double logaMax,
		       double *Tauq, double *H, double *Dh);



#endif 	    /* !APPROXLEGENDRE_H_ */
