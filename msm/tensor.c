/*	tensor.c.  Version del 2 de Septiembre, 2004	*/

#define TENSOR_C

#ifndef OPERACIONES_C
#include <operaciones.c>
#endif

/*     Function prototypes    */

double ****reservar_cuadritensor( int wdim, int zdim, int ydim, int xdim);
void liberar_cuadritensor( double ****m, int wdim, int zdim, int ydim);
double ***reservar_tritensor( int zdim, int ydim, int xdim);
double ***redimensiona_tritensor( int zdim0, int ydim0, int xdim0,
	int zdim, int ydim, int xdim, double ***pointer);
void liberar_tritensor( double ***m, int zdim, int ydim);
float ***reservar_tritensor_float( int zdim, int ydim, int xdim);
void liberar_tritensor_float( float ***m, int zdim, int ydim);
int ***reservar_tritensor_int( int zdim, int ydim, int xdim);
void liberar_tritensor_int( int ***m, int zdim, int ydim);
char ***reservar_tritensor_char( int zdim, int ydim, int xdim);
void liberar_tritensor_char( char ***m, int zdim, int ydim);
double **reservar_matriz( int ydim, int xdim);
void liberar_matriz( double **m, int ydim);
float **reservar_matriz_float( int ydim, int xdim);
void liberar_matriz_float( float **m, int ydim);
int **reservar_matriz_int( int ydim, int xdim);
void liberar_matriz_int( int **m, int ydim);
char **reservar_matriz_char( int ydim, int xdim);
void liberar_matriz_char( char **m, int ydim);
size_t **reservar_matriz_size_t( int ydim, int xdim);
void liberar_matriz_size_t( size_t **m, int ydim);


void limpia(int dimx, int dimy, double **matriz);
void limpia_int(int dimx, int dimy, int **matriz);
void limpia_char(int dimx, int dimy, char val, char **matriz);
void limpia_lista( int dimx, double *vector);
void limpia_char_lista(int dimx, char val, char *vector);


void asigna(int dimx, int dimy, double **fuente, double **destino);
void asigna_char(int dimx, int dimy, char **fuente, char **destino);
void asigna_resta(int dimx, int dimy, double **fuente, double **destino);
void asigna_suma(int dimx, int dimy, double **fuente, double **destino);
void asigna_combina(int dimx, int dimy, double escalar, double **fuente,
	double **destino);
void asigna_desplaza(int dimx, int dimy, double shift, double **funcion);
void asigna_escala(int dimx, int dimy, double scale, double **funcion);
void asigna_divide(int dimx, int dimy, double **fuente, double **destino);
void asigna_multiplica(int dimx, int dimy, double **fuente, double **destino);
void asigna_divide_vec(int dimx, int dimy, double **fx, double **fy, 
	double **dx, double **dy);
void asigna_multiplica_vec(int dimx, int dimy, double **fx, double **fy, 
	double **dx, double **dy);
void asigna_modulo( int dimx, int dimy, double **vx, double **vy,
	double **mod);

void asigna_lista( int dimx, double *fuente, double *destino);
void asigna_resta_lista( int dimx, double *fuente, double *destino);
void asigna_divide_lista(int dimx, double *fuente, double *destino);
void asigna_multiplica_lista(int dimx, double *fuente, double *destino);


double media( int dimx, int dimy, double **data);
double dispersion( int dimx, int dimy, double **cont);
double anorma( int dimx, int dimy, double **data);
void denorma( int dimx, int dimy, double norma, double **data);
double anorma1( int dimx, int dimy, double **data);
double anorma2( int dimx, int dimy, double **data);
void extrema( int dimx, int dimy, double **data, double *extr);

double dispersion_vec( int dimx, int dimy, double **vx, double **vy);

void normaliza_vector( int dimx, int dimy, double **vx, double **vy);
void gira(int dimx, int dimy, double angle, double **vx, double **vy);

double media_lista( int dimx, double *cont);
double anorma_lista( int dimx,  double *cont);
void denorma_lista( int dimx, double norma, double *cont);
double dispersion_lista( int dimx, double *cont);
double covarianza_lista( int dimx, double *x, double *y);
void extrema_lista( int dimx, double *data, double *extr);
double anorma1_lista( int dimx, double *cont);

void coarse_resolution( int dimx, int dimy, double block, double **source, 
			double **target);
void coarse_resolution_char( int dimx, int dimy, double block, char **source, 
			     char **target);
void coarse_resolution_mask( int dimx, int dimy, double block, int ivabs,
			     char **source, char **target);
void blur_resolution( int dimx, int dimy, int block, double **source, 
	double **target);
void coarse_resolution_lista( int dimx, double block, double *source,
			      double *target);

void de4(int dimx, int dimy, double **signal);	
void de4_mirror(int dimx, int dimy, double **signal);
void recorta_ventana( int dimx, int dimy, int ix0, int iy0, int dx, 
	int dy, double **data, double **window);
void remete_ventana( int dimx, int dimy, int ix0, int iy0, int ix1, 
		     int iy1, int tx, int ty, double **window, double **total);
void recorta_ventana_char( int dimx, int dimy, int ix0, int iy0, int dx, 
			   int dy, char **data, char **window);


int extrae_extension( char *nombre, char separador, char *ext);



/*   Function declarations   */

double ****reservar_cuadritensor( int wdim, int zdim, int ydim, int xdim)
{
	double ****m;
	int j,k,l;

	m=(double ****) calloc(wdim,sizeof(double ***));
	for(l=0;l<wdim;++l)
	{
	  m[l]=(double ***) calloc(zdim,sizeof(double **));
	  for(k=0;k<zdim;++k) 
	  {
		m[l][k]=(double **) calloc( ydim,sizeof(double *));
		for(j=0;j<ydim;++j) m[l][k][j]=(double *) calloc(xdim,sizeof(double));
	}
	}
	return m;
}

void liberar_cuadritensor( double ****m, int wdim, int zdim, int ydim)
{
	int j,k,l;

	for(l=0; l<wdim; ++l)
	{
	  for(k=0; k<zdim; ++k)
	  {
		for(j=0; j<ydim; ++j) free(m[l][k][j]);
		free(m[l][k]);
	  }
	  free(m[l]);
	}
	free(m);
}

double ***reservar_tritensor( int zdim, int ydim, int xdim)
{
	double ***m;
	int j,k;

	m=(double ***) calloc(zdim,sizeof(double **));
	for(k=0;k<zdim;++k) 
	{
		m[k]=(double **) calloc( ydim,sizeof(double *));
		for(j=0;j<ydim;++j) m[k][j]=(double *) calloc(xdim,sizeof(double));
	}
	
	return m;
}

double ***redimensiona_tritensor( int zdim0, int ydim0, int xdim0,
	int  zdim, int ydim, int xdim, double ***pointer)
{
	double ***m;
	int i,j,k;

	m=reservar_tritensor(zdim,ydim,xdim);

	for(k=0;k<zdim;k++)
	{
	for(j=0;j<ydim;j++)
	{
	for(i=0;i<xdim;i++)
	{
		if((k<zdim0)&&(j<ydim0)&&(i<xdim0))
			m[k][j][i]=pointer[k][j][i];
		else m[k][j][i]=0.;
	}
	}
	}
	liberar_tritensor(pointer,zdim,ydim);

	return m;
}

void liberar_tritensor( double ***m, int zdim, int ydim)
{
	int j,k;

	for(k=0; k<zdim; ++k)
	{
		for(j=0; j<ydim; ++j) free(m[k][j]);
		free(m[k]);
	}
	free(m);
}

float ***reservar_tritensor_float( int zdim, int ydim, int xdim)
{
	float ***m;
	int j,k;

	m=(float ***) calloc(zdim,sizeof(float **));
	for(k=0;k<zdim;++k) 
	{
		m[k]=(float **) calloc( ydim,sizeof(float *));
		for(j=0;j<ydim;++j) m[k][j]=(float *) calloc(xdim,sizeof(float));
	}
	
	return m;
}

void liberar_tritensor_float( float ***m, int zdim, int ydim)
{
	int j,k;

	for(k=0; k<zdim; ++k)
	{
		for(j=0; j<ydim; ++j) free(m[k][j]);
		free(m[k]);
	}
	free(m);
}

int ***reservar_tritensor_int( int zdim, int ydim, int xdim)
{
	int ***m;
	int j,k;

	m=(int ***) calloc(zdim,sizeof(int **));
	for(k=0;k<zdim;++k) 
	{
		m[k]=(int **) calloc( ydim,sizeof(int *));
		for(j=0;j<ydim;++j) m[k][j]=(int *) calloc(xdim,sizeof(int));
	}
	
	return m;
}

void liberar_tritensor_int( int ***m, int zdim, int ydim)
{
	int j,k;

	for(k=0; k<zdim; ++k)
	{
		for(j=0; j<ydim; ++j) free(m[k][j]);
		free(m[k]);
	}
	free(m);
}

char ***reservar_tritensor_char( int zdim, int ydim, int xdim)
{
	char ***m;
	int j,k;

	m=(char ***) calloc(zdim,sizeof(char **));
	for(k=0;k<zdim;++k) 
	{
		m[k]=(char **) calloc( ydim,sizeof(char *));
		for(j=0;j<ydim;++j) m[k][j]=(char *) calloc(xdim,sizeof(char));
	}
	
	return m;
}

void liberar_tritensor_char( char ***m, int zdim, int ydim)
{
	int j,k;

	for(k=0; k<zdim; ++k)
	{
		for(j=0; j<ydim; ++j) free(m[k][j]);
		free(m[k]);
	}
	free(m);
}


double **reservar_matriz( int ydim, int xdim)
{
	double **m;
	int j;

	m=(double **) calloc(ydim,sizeof(double *));
	for(j=0;j<ydim;++j) m[j]=(double *) calloc( xdim,sizeof(double));
	
	return m;
}

void liberar_matriz( double **m, int ydim)
{
	int j;

	for(j=0; j<ydim; ++j) free(m[j]);
	free(m);

}

float **reservar_matriz_float( int ydim, int xdim)
{
	float **m;
	int j;

	m=(float **) calloc(ydim,sizeof(float *));
	for(j=0;j<ydim;++j) m[j]=(float *) calloc( xdim,sizeof(float));
	
	return m;
}

void liberar_matriz_float( float **m, int ydim)
{
	int j;

	for(j=0; j<ydim; ++j) free(m[j]);
	free(m);
}


int **reservar_matriz_int( int ydim, int xdim)
{
	int **m;
	int j;

	m=(int **) calloc(ydim,sizeof(int *));
	for(j=0;j<ydim;++j) m[j]=(int *) calloc( xdim,sizeof(int));
	
	return m;
}

void liberar_matriz_int( int **m, int ydim)
{
	int j;

	for(j=0; j<ydim; ++j) free(m[j]);
	free(m);
}


char **reservar_matriz_char( int ydim, int xdim)
{
	char **m;
	int j;

	m=(char **) calloc(ydim,sizeof(char *));
	for(j=0;j<ydim;++j) m[j]=(char *) calloc( xdim,sizeof(char));
	
	return m;
}

void liberar_matriz_char( char **m, int ydim)
{
	int j;

	for(j=0; j<ydim; ++j) free(m[j]);
	free(m);
}


size_t **reservar_matriz_size_t( int ydim, int xdim)
{
	size_t **m;
	int j;

	m=(size_t **) calloc(ydim,sizeof(size_t *));
	for(j=0;j<ydim;++j) m[j]=(size_t *) calloc( xdim,sizeof(size_t));
	
	return m;
}

void liberar_matriz_size_t( size_t **m, int ydim)
{
	int j;

	for(j=0; j<ydim; ++j) free(m[j]);
	free(m);
}

void limpia(int dimx, int dimy, double **matriz)
{
	int ix,iy;

	for(ix=0;ix<dimx;ix++)
	{
	for(iy=0;iy<dimy;iy++)
	{
		matriz[iy][ix]=0.;
	}
	}
}

void limpia_int(int dimx, int dimy, int **matriz)
{
	int ix,iy;

	for(ix=0;ix<dimx;ix++)
	{
	for(iy=0;iy<dimy;iy++)
	{
		matriz[iy][ix]=0;
	}
	}
}

void limpia_char(int dimx, int dimy, char val, char **matriz)
{
	int ix,iy;

	for(ix=0;ix<dimx;ix++)
	{
	for(iy=0;iy<dimy;iy++)
	{
		matriz[iy][ix]=val;
	}
	}
}

void limpia_lista( int dimx, double *vector)
{
	int ix;

	for(ix=0;ix<dimx;ix++) vector[ix]=0.;
}

void limpia_char_lista(int dimx, char val, char *vector)
{
	int ix;

	for(ix=0;ix<dimx;ix++) vector[ix]=val;
}


void asigna(int dimx, int dimy, double **fuente, double **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		destino[iy][ix]=fuente[iy][ix];
	}
	}
}

void asigna_char(int dimx, int dimy, char **fuente, char **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		destino[iy][ix]=fuente[iy][ix];
	}
	}
}

void asigna_resta(int dimx, int dimy, double **fuente, double **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		destino[iy][ix]-=fuente[iy][ix];
	}
	}
}

void asigna_suma(int dimx, int dimy, double **fuente, double **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		destino[iy][ix]+=fuente[iy][ix];
	}
	}
}


void asigna_combina(int dimx, int dimy, double escalar, double **fuente,
	double **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		destino[iy][ix]+=escalar*fuente[iy][ix];
	}
	}
}


void asigna_desplaza(int dimx, int dimy, double shift, double **funcion)
{

	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		funcion[iy][ix]+=shift;
	}
	}

}

void asigna_escala(int dimx, int dimy, double scale, double **funcion)
{

	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		funcion[iy][ix]=scale*funcion[iy][ix];
	}
	}

}


void asigna_divide(int dimx, int dimy, double **fuente, double **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		if(fabs(fuente[iy][ix])>1e-20)
			destino[iy][ix]=destino[iy][ix]/fuente[iy][ix];
		else destino[iy][ix]=0.;
	}
	}
}

void asigna_multiplica(int dimx, int dimy, double **fuente, double **destino)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		destino[iy][ix]=destino[iy][ix]*fuente[iy][ix];
	}
	}
}

void asigna_divide_vec(int dimx, int dimy, double **fx, double **fy, 
	double **dx, double **dy)
{
	double buffx,buffy,mod;
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod=fx[iy][ix]*fx[iy][ix]+fy[iy][ix]*fy[iy][ix];
		if(mod>1e-30)
		{
			buffx=(fx[iy][ix]*dx[iy][ix]+
				fy[iy][ix]*dy[iy][ix])/mod;
			buffy=(fx[iy][ix]*dy[iy][ix]-
				fy[iy][ix]*dx[iy][ix])/mod;
			dx[iy][ix]=buffx;
			dy[iy][ix]=buffy;
		}
		else
		{
			dx[iy][ix]=0.;
			dy[iy][ix]=0.;
		}
	}
	}
}

void asigna_multiplica_vec(int dimx, int dimy, double **fx, double **fy, 
	double **dx, double **dy)
{
	double buffx,buffy;
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		buffx=fx[iy][ix]*dx[iy][ix]-fy[iy][ix]*dy[iy][ix];
		buffy=fx[iy][ix]*dy[iy][ix]+fy[iy][ix]*dx[iy][ix];
		dx[iy][ix]=buffx;
		dy[iy][ix]=buffy;
	}
	}
}

void asigna_modulo( int dimx, int dimy, double **vx, double **vy,
	double **mod)
{
	int ix,iy;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod[iy][ix]=sqrt(vx[iy][ix]*vx[iy][ix]
				+vy[iy][ix]*vy[iy][ix]);
	}
	}
}

void asigna_lista( int dimx, double *fuente, double *destino)
{
	int ix;

	for(ix=0;ix<dimx;ix++) destino[ix]=fuente[ix];
}

void asigna_resta_lista( int dimx, double *fuente, double *destino)
{
	int ix;

	for(ix=0;ix<dimx;ix++) destino[ix]-=fuente[ix];
}

void asigna_divide_lista(int dimx, double *fuente, double *destino)
{
	int ix;

	for(ix=0;ix<dimx;ix++)
	{
		if(fabs(fuente[ix])>1e-20)
			destino[ix]=destino[ix]/fuente[ix];
		else destino[ix]=0.;
	}
}

void asigna_multiplica_lista(int dimx, double *fuente, double *destino)
{
	int ix;

	for(ix=0;ix<dimx;ix++) destino[ix]=destino[ix]*fuente[ix];
}


double media( int dimx, int dimy, double **data)
{
	double suma;
	int ix,iy;

	suma=0.;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  suma+=data[iy][ix];
	}
	}
	suma=suma/((double) dimx*dimy);
	return(suma);
}

double dispersion( int dimx, int dimy, double **cont)
{
	int ix,iy;
	double sumaxx,sumax;

	sumax=0.;
	sumaxx=0.;
	for(ix=0;ix<dimx;ix++)
	{
	for(iy=0;iy<dimy;iy++)
	{
		sumax+=cont[iy][ix];
		sumaxx+=cont[iy][ix]*cont[iy][ix];
	}
	}
	sumax=sumax/((double) dimx*dimy);
	sumaxx=sumaxx/((double) dimx*dimy)-sumax*sumax;
	sumaxx=sqrt(fabs(sumaxx));

	return(sumaxx);
}

double anorma( int dimx, int dimy, double **data)
{
	double suma;
	int ix,iy;

	suma=0.;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  suma+=data[iy][ix];
	}
	}
	suma=suma/((double) dimx*dimy);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  data[iy][ix]-=suma;
	}
	}

	return(suma);
}

void denorma( int dimx, int dimy, double norma, double **data)
{
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    data[iy][ix]+=norma;
	}
	}
}

double anorma1( int dimx, int dimy, double **data)
{
	double suma;
	int ix,iy;

	suma=0.;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  suma+=fabs(data[iy][ix]);
	}
	}
	suma=suma/((double) dimx*dimy);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  data[iy][ix]/=suma;
	}
	}

	return(suma);
}

double anorma2( int dimx, int dimy, double **data)
{
	double suma;
	int ix,iy;

	suma=0.;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  suma+=data[iy][ix]*data[iy][ix];
	}
	}
	suma=sqrt(suma/((double) dimx*dimy));

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  data[iy][ix]/=suma;
	}
	}

	return(suma);
}

void extrema( int dimx, int dimy, double **data, double *extr)
{
	int ix,iy;

	extr[0]=data[0][0];
	extr[1]=data[0][0];

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  extr[0]=fMin(extr[0],data[iy][ix]);
	  extr[1]=fMax(extr[1],data[iy][ix]);
	}
	}

}

double dispersion_vec( int dimx, int dimy, double **vx, double **vy)
{
	double dx,dy,out;

	dx=dispersion(dimx,dimy,vx);
	dy=dispersion(dimx,dimy,vy);

	out=sqrt(dx*dx+dy*dy);

	return(out);
}

void normaliza_vector( int dimx, int dimy, double **vx, double **vy)
{
	double mod;
	int ix,iy;

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod=sqrt(vx[iy][ix]*vx[iy][ix]+vy[iy][ix]*vy[iy][ix]);
		if(mod>1e-30)
		{
			vx[iy][ix]=vx[iy][ix]/mod;
			vy[iy][ix]=vy[iy][ix]/mod;
		}
		else
		{
			vx[iy][ix]=0.;
			vy[iy][ix]=0.;
		}
	}
	}

}

void gira(int dimx, int dimy, double angle, double **vx, double **vy)
{
        int ix,iy;
	double normax,normay;
        double buffx,buffy;

	for(iy=0;iy<dimy;iy++)
        {
	for(ix=0;ix<dimx;ix++)
        {
		buffx=vx[iy][ix]*cos(angle)-vy[iy][ix]*sin(angle);
                buffy=vx[iy][ix]*sin(angle)+vy[iy][ix]*cos(angle);
		vx[iy][ix]=buffx;
		vy[iy][ix]=buffy;
	}
	}

}


double media_lista( int dimx, double *cont)
{
	double suma;
	int ix;

	suma=0.;
	for(ix=0;ix<dimx;ix++) suma+=cont[ix];
	suma=suma/((double) dimx);

	return(suma);
}

double anorma_lista( int dimx, double *cont)
{
	double media;
	int ix;

	media=0.;
	for(ix=0;ix<dimx;ix++)
	{
		media+=cont[ix];
	}
	media=media/((double)dimx);

	for(ix=0;ix<dimx;ix++) cont[ix]-=media;

	return(media);
}


void denorma_lista( int dimx, double norma, double *cont)
{
	int ix;

	for(ix=0;ix<dimx;ix++) cont[ix]+=norma;

}


double dispersion_lista( int dimx, double *cont)
{
	double sumaxx,sumax;
	int ix;

	sumax=0.;
	sumaxx=0.;
	for(ix=0;ix<dimx;ix++)
	{
		sumax+=cont[ix];
		sumaxx+=cont[ix]*cont[ix];
	}
	sumax=sumax/((double) dimx);
	sumaxx=sumaxx/((double) dimx)-sumax*sumax;
	sumaxx=sqrt(fabs(sumaxx));

	return(sumaxx);
}

double covarianza_lista( int dimx, double *x, double *y)
{
	double sumaxy,sumax,sumay;
	int ix;

	sumax=0.;
	sumay=0.;
	sumaxy=0.;
	for(ix=0;ix<dimx;ix++)
	{
		sumax+=x[ix];
		sumay+=y[ix];
		sumaxy+=x[ix]*y[ix];
	}
	sumax=sumax/((double) dimx);
	sumay=sumay/((double) dimx);
	sumaxy=sumaxy/((double) dimx)-sumax*sumay;

	return(sumaxy);
}


void extrema_lista( int dimx, double *data, double *extr)
{
	int ix;

	extr[0]=data[0];
	extr[1]=data[0];

	for(ix=0;ix<dimx;ix++)
	{
		extr[0]=fMin(extr[0],data[ix]);
		extr[1]=fMax(extr[1],data[ix]);
	}


}

double anorma1_lista( int dimx, double *cont)
{
        double norma;
	int ix;

	norma=0.;
	for(ix=0;ix<dimx;ix++) norma+=cont[ix];
	norma/=(double)dimx;
	for(ix=0;ix<dimx;ix++) cont[ix]/=norma;

	return(norma);
}


void coarse_resolution( int dimx, int dimy, double block, double **source,
			double **target)
{
	double bsize;
	int beff,xsize,ysize;
        int ix,iy;
        int bx,by;


	if(block<1.)
	{
		beff=(int) (1./block);
		for(ix=0;ix<dimx;ix++)
		{
		for(iy=0;iy<dimy;iy++)
		{
			for(by=0;by<beff;by++)
			{
			for(bx=0;bx<beff;bx++)
			{
				target[iy*beff+by][ix*beff+bx]=
					source[iy][ix];
			}
			}
		}
		}

	}
	else
	{
		beff=(int) block;
		bsize=(double)(beff*beff);
		xsize=dimx/beff;
		ysize=dimy/beff;
        	for(iy=0;iy<ysize;iy++)
		{
		for(ix=0;ix<xsize;ix++)
		{
                	target[iy][ix]=0.;
			for(by=0;by<beff;by++)
                	{
                	for(bx=0;bx<beff;bx++)
                	{
                	        target[iy][ix]+=
					source[beff*iy+by][beff*ix+bx];
                	}
                	}
                	target[iy][ix]=target[iy][ix]/bsize;
		}
		}
	}


}

void coarse_resolution_char( int dimx, int dimy, double block, char **source,
			     char **target)
{
    double bsize;
    double cumul;
    int beff,xsize,ysize;
    int ix,iy;
    int bx,by;


    if(block<1.)
    {
	beff=(int) (1./block);
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    for(by=0;by<beff;by++)
	    {
		for(bx=0;bx<beff;bx++)
		{
		    target[iy*beff+by][ix*beff+bx]=source[iy][ix];
		}
	    }
	}
	}
	
    }
    else
    {
	beff=(int) block;
	bsize=(double)(beff*beff);
	xsize=dimx/beff;
	ysize=dimy/beff;
	for(iy=0;iy<ysize;iy++)
	{
	for(ix=0;ix<xsize;ix++)
	{
	    cumul=0.;
	    for(by=0;by<beff;by++)
	    {
	    for(bx=0;bx<beff;bx++)
	    {
		cumul+=(double)source[beff*iy+by][beff*ix+bx]
		    +((source[beff*iy+by][beff*ix+bx]<0)?256.:0.);
	    }
	    }
	    target[iy][ix]=(char)(cumul/bsize);
	}
	}
    }


}

void coarse_resolution_mask( int dimx, int dimy, double block, int ivabs,
			     char **source, char **target)
{
    double bsize;
    double cumul;

    int *voting;
    int beff,xsize,ysize;
    int ix,iy;
    int bx,by;
    int voting0;
    int iv,iv0;


    if(block<1.)
    {
	beff=(int) (1./block);
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    for(by=0;by<beff;by++)
	    {
	    for(bx=0;bx<beff;bx++)
	    {
		target[iy*beff+by][ix*beff+bx]=source[iy][ix];
	    }
	    }
	}
	}
    }
    else
    {
	voting=(int *) calloc(256,sizeof(int));
	beff=(int) block;
	bsize=(double)(beff*beff);
	xsize=dimx/beff;
	ysize=dimy/beff;
	for(iy=0;iy<ysize;iy++)
	{
	for(ix=0;ix<xsize;ix++)
	{
	    cumul=0.;
	    for(iv=0;iv<256;iv++) voting[iv]=0;
	    for(by=0;by<beff;by++)
	    {
	    for(bx=0;bx<beff;bx++)
	    {
		iv=Mod((int)source[beff*iy+by][beff*ix+bx],256);
		voting[iv]++;
	    }
	    }

	    iv0=-1;
	    voting0=0;
	    if(ivabs>0)
	    {
		if(voting[ivabs]>0) iv0=ivabs;
	    }
	    if(iv0==-1)
	    {
		for(iv=0;iv<256;iv++)
		{
		    if(voting[iv]>voting0)
		    {
			voting0=voting[iv];
			iv0=iv;
		    }
		}
	    }
	    target[iy][ix]=(char)iv0;
	}
	}
	free(voting);
    }

}

void blur_resolution( int dimx, int dimy, int block, double **source, double
	**target)
{
        double bsize;
	int *iix,*iiy;
        int ix,iy,dx,dy;

	iix=(int *)calloc(2*block+1,sizeof(int));
	iiy=(int *)calloc(2*block+1,sizeof(int));

        bsize=(double)(2*block+1);
	bsize=bsize*bsize;
        for(iy=0;iy<dimy;iy++)
        {
		for(dy=0;dy<2*block+1;dy++) iiy[dy]=Mod(iy+dy-block,dimy);
	        for(ix=0;ix<dimx;ix++)
        	{
			for(dx=0;dx<2*block+1;dx++)
				iix[dx]=Mod(ix+dx-block,dimx);
	                target[iy][ix]=0.;
        	        for(dy=0;dy<2*block+1;dy++)
			{
			for(dx=0;dx<2*block+1;dx++)
			{
				target[iy][ix]+=source[iiy[dy]][iix[dx]];
			}
			}
			target[iy][ix]=target[iy][ix]/bsize;
        }
        }

	free(iix);
	free(iiy);

}

void coarse_resolution_lista( int dimx, double block, double *source,
	double *target)
{
	double bsize;
	int beff,size;
        int ix,bx;


	if(block<1.)
	{
	  beff=(int) (1./block);
	  for(ix=0;ix<dimx;ix++)
	  {
	    for(bx=0;bx<beff;bx++) target[ix*beff+bx]=source[ix];
	  }		
	}
	else
	{
	  beff=(int) block;
	  bsize=(double)beff;
	  size=dimx/beff;
	  for(ix=0;ix<size;ix++)
	  {
	    target[ix]=0.;
	    for(bx=0;bx<beff;bx++) target[ix]+=source[beff*ix+bx];
	    target[ix]=target[ix]/bsize;
	  }
	}


}

void de4(int dimx, int dimy, double **signal)
{
	int ix,iy;
	double **aux;
	aux=reservar_matriz(dimy/2,dimx/2);

	coarse_resolution(dimx,dimy,2.,signal,aux);
	for(iy=0;iy<dimy/2;iy++)
	{
	for(ix=0;ix<dimx/2;ix++)
	{
		signal[iy][ix]=aux[iy][ix];
		signal[iy][ix+dimx/2]=aux[iy][ix];
		signal[iy+dimy/2][ix]=aux[iy][ix];
		signal[iy+dimy/2][ix+dimx/2]=aux[iy][ix];
	}
	}

	liberar_matriz(aux,dimy/2);

}

void de4_mirror(int dimx, int dimy, double **signal)
{
	int ix,iy;
	double **aux;
	aux=reservar_matriz(dimy/2,dimx/2);

	coarse_resolution(dimx,dimy,2.,signal,aux);
	for(iy=0;iy<dimy/2;iy++)
	{
	for(ix=0;ix<dimx/2;ix++)
	{
		signal[iy][ix]=aux[iy][ix];
		signal[iy][ix+dimx/2]=aux[iy][dimx/2-1-ix];
		signal[iy+dimy/2][ix]=aux[dimy/2-1-iy][ix];
		signal[iy+dimy/2][ix+dimx/2]=aux[dimy/2-1-iy][dimx/2-1-ix];
	}
	}

	liberar_matriz(aux,dimy/2);

}



void recorta_ventana( int dimx, int dimy, int ix0, int iy0, int dx, 
	int dy, double **data, double **window)
{
	int *iix,*iiy;
	int ix,iy;

	iix=(int *) calloc(dx,sizeof(int));
	iiy=(int *) calloc(dy,sizeof(int));

	for(ix=0;ix<dx;ix++) iix[ix]=Mod(ix0+ix,dimx);
	for(iy=0;iy<dy;iy++) iiy[iy]=Mod(iy0+iy,dimy);

	for(iy=0;iy<dy;iy++)
	{
	for(ix=0;ix<dx;ix++)
	{
		window[iy][ix]=data[iiy[iy]][iix[ix]];
	}
	}

	free(iix);
	free(iiy);
}


void remete_ventana( int dimx, int dimy, int ix0, int iy0, int ix1, 
	int iy1, int tx, int ty, double **window, double **total)
{
	int *iix,*iiy;
	int ix,iy;

	iix=(int *) calloc(tx,sizeof(int));
	iiy=(int *) calloc(ty,sizeof(int));

	for(ix=0;ix<tx;ix++) iix[ix]=Mod(ix0+ix,dimx);
	for(iy=0;iy<ty;iy++) iiy[iy]=Mod(iy0+iy,dimy);

	for(iy=0;iy<ty;iy++)
	{
	for(ix=0;ix<tx;ix++)
	{
		total[iiy[iy]][iix[ix]]=window[iy1+iy][ix1+ix];
	}
	}
}

void recorta_ventana_char( int dimx, int dimy, int ix0, int iy0, int dx, 
	int dy, char **data, char **window)
{
	int *iix,*iiy;
	int ix,iy;

	iix=(int *) calloc(dx,sizeof(int));
	iiy=(int *) calloc(dy,sizeof(int));

	for(ix=0;ix<dx;ix++) iix[ix]=Mod(ix0+ix,dimx);
	for(iy=0;iy<dy;iy++) iiy[iy]=Mod(iy0+iy,dimy);

	for(iy=0;iy<dy;iy++)
	{
	for(ix=0;ix<dx;ix++)
	{
		window[iy][ix]=data[iiy[iy]][iix[ix]];
	}
	}

	free(iix);
	free(iiy);
}

int extrae_extension( char *nombre, char separador, char *ext)
{
    int longo;
    int ic,ic0;

    longo=strlen(nombre);
    for(ic=longo-1;(ic>=0)&&(nombre[ic]!=separador);ic--) ;

    ic0=ic+1;
    for(ic=ic0;ic<longo;ic++) ext[ic-ic0]=nombre[ic];
    ext[ic-ic0]='\0';
	
    return(longo-ic0);
}
