/********************************************************/
/*                                                      */
/*    utils.h - Version del 2 de Dekembriou, 2004      */
/*                                                      */
/********************************************************/

#ifndef   	UTILS_H
# define   	UTILS_H

void registra_WTMM(  int i0, int n, char *nombre, double *x, double *y );

int registra_WTcolor( const char *colname, double *wtcolor, 
		      int i, int nx, int ny );

int registra_PF( char *pfname, double **Z, 
		 int nq, double *qArray, int nws, double *wsArray );
int graba_PF( const char *pfname, double **Z, 
	      int nq, int nws, double *wsArray );

void lee_serie( int dimx, char *nombre_in, double **datos);

int lee_dimension_serie( char *nombre_in );

void asigna_desplaza_lista(int dimx, double shift, double *funcion);

#ifndef FLAG_MSM

void limpia(int dimx, int dimy, double **matriz);

void asigna_lista( int dimx, double *fuente, double *destino);

void lee_datos( int leff, char *nombre_in, double **signal);

void graba_datos( int leff, char *nombre_in, double **signal);

void graba_serie( int leff, char *nombre_in, double *datos);

double **reservar_matriz( int ydim, int xdim);

void liberar_matriz( double **m, int ydim);

float **reservar_matriz_float( int ydim, int xdim);

void liberar_matriz_float( float **m, int ydim);

int **reservar_matriz_int( int ydim, int xdim);

void liberar_matriz_int( int **m, int ydim);

char **reservar_matriz_char( int ydim, int xdim);

void liberar_matriz_char( char **m, int ydim);

#endif


#endif 	    /* !UTILS_H */
