/********************************************************/
/*                                                      */
/*    utils.c - Version del 2 de Dekembriou, 2004       */
/*                                                      */
/* Useful functions for allocations, storage, display   */
/*                                                      */
/********************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

/* Personal libraries */
#include <const_wtmm.h>
#include <utils.h>

extern int VERBOSE;
extern int D_space;


/* ===================================================================== */
void registra_WTMM(  int i0, int n, char *nombre, double *x, double *y ) {
  /* ===================================================================== 
   *
   * I'm sorry, but don't expect nothing more to be stored for the moment...
   * 
   * ===================================================================== */

  FILE *canal;
  int ix;

  if( (canal=fopen(nombre,"wt")) == NULL) {
    fprintf(stderr, "\n Error while opening file %s", nombre);
    exit(-1);
  }
  
  for( ix=i0; ix<n; ix++ ) 
    fprintf(canal,"%f  %f\n",x[ix],y[ix]);
  
  
  if(canal) fclose(canal);
}

/* ===================================================================== */
int registra_WTcolor( const char *colname, double *wtcolor, 
		      int i, int nx, int ny ) {
  /* ===================================================================== 
   * 
   * Prints the color cascade in a ppm format.
   * Subroutine for constructing color coded 3D wavelet decomposition: 
   * x axis is time; y axis is wavelet  scale; z axis is wavelet 
   * 
   * coefficient (bright colors represent large values).
   * =====================================================================
   */
  
  int l, j, bin, color1, color2, color3;
  int NG1, NG2, NG3, NG4, NG5, NG6;
  double delta;
  double v, max = -1e20, min = 1e20;
  FILE *colfile;
  
  /* open the file for color representation of WT */
  if( (colfile=fopen(colname, "wt")) == NULL) {
    fprintf(stderr, "\n Error while opening color file %s", colname);
    exit(-1);
  }   
  
  /* color levels */ 
  NG1 = NG; NG2 = 2*NG; NG3 = 3*NG;
  NG4 = 4*NG; NG5 = 5*NG; NG6 = 6*NG;
  
  /* Find min and max of the representation */
  l = 0;
  while (l < i) {
    v=wtcolor[l++];
    max = MAX(v, max);
    min = MIN(v, min);
  }
  
  if(DEBUG) printf("\n\t Color: min=%g - max=%g",min,max);
  /* step for color distinction */
  delta = (max - min)/(double)(NG6-2);
  
  fprintf(colfile,"P3\n# CREATOR: multifractal\n%d %d\n255\n", nx, ny);
  
  /* Represent the different values in color */
  j = l = 0;
  while (l < i) {
    bin  =  (int)(((wtcolor[l] - min)/delta));
    
    if (bin < NG1) {
      color1 = NG;
      color2 = bin;
      color3 = 0;
    }
    else if (bin < NG2) {
      color1 = NG2-bin;
      color2 = NG;
      color3 = 0;
    }
    else if (bin < NG3) {
      color1 = 0;
      color2 = NG;
      color3 = bin-NG2;
    }
    else if (bin < NG4) {
      color1 = 0;
      color2 = NG4-bin;
      color3 = NG;
    }
    else if (bin < NG5) {
      color1 = bin-NG4;
      color2 = 0;
      color3 = NG;
    }
    else {
      color1 = NG;
      color2 = 0;
      color3 = NG6-bin;
    }
    fprintf(colfile,"%3d %3d %3d ", color1, color2, color3);
    l++;
    if (++j >= 5) {
      fprintf(colfile,"\n");
      j = 0;
    }
  }
    
  if(colfile)  fclose(colfile);
  return 0;
}


/* ===================================================================== */
int registra_PF( char *pfname, double **Z, 
		 int nq, double *qArray, int nws, double *wsArray ) {
  /* ===================================================================== 
   *   
   * Register the factors of the partition function
   * 
   * ===================================================================== */
  
  FILE *pffile;
  int ind_q, ind_ws;
  double ws;
  
  /* open the file for partition function representation */
  if( (pffile=fopen(pfname, "wt")) == NULL) {
    fprintf(stderr, "\n Error while opening partition function file %s", 
	    pfname);
    exit(-1);
  }   
  
  if(VERBOSE) 
    printf("\n Writing the partition function in %s...", pfname);
  
  /* Partition Function representation - write the moments q */
  fprintf(pffile,"# %f\n",qArray[0]);
  
  /* Print the partition function  - write the partition function for each
   * scale */
  for (ind_ws=1; ind_ws<nws; ind_ws++) {
    /* we don't dispay ws=min_ws */
    ws = wsArray[ind_ws];
    
    fprintf(pffile,"%g ", LOG(ws));
    for (ind_q = 0; ind_q < nq; ind_q++)
      fprintf(pffile,"%g ", LOG(Z[ind_ws][ind_q]));
    fprintf(pffile,"\n");
    
  }
  
  if(pffile)  fclose(pffile);
  return 0;
}

/* ===================================================================== */
int graba_PF( const char *pfname, double **Z, 
	      int nq, int nws, double *wsArray ) {
  /* ===================================================================== 
   *   
   * Register the factors of the partition function in a binary file
   * 
   * ===================================================================== */
  
  FILE *pffile;
  int  iws, iq;
  
  /* open the file for partition function representation */
  if( (pffile=fopen(pfname, "wb")) == NULL) {
    fprintf(stderr, "\n Error while opening partition function file %s", 
	    pfname);
    exit(-1);
  }   
  
  fwrite(wsArray,sizeof(double),nws,pffile);
  
  if(VERBOSE) 
    printf("\n Writing the partition function in %s...", pfname);
  
  for (iws=0; iws<nws; iws++)
    fwrite(Z[iws],sizeof(double),nq,pffile);
  
  if(pffile)  fclose(pffile);
  return 0;
}


/* ===================================================================== */
int lee_dimension_serie( char *nombre_in ) {
  /* ===================================================================== 
   *
   * Returns the number of data stored in a text file 
   *
   * ===================================================================== */

  FILE *canal;
  char nombre[90];
  double dato;
  int dimx=0;
  
  sprintf(nombre,"%s.txt",nombre_in);
  if( (canal=fopen(nombre,"rt")) == NULL) {
    fprintf(stderr,"\n Error in lee_dimension_serie while reading file %s", nombre);
    exit(-1);
  }

  while(fscanf(canal,"%lf",&dato) != EOF) dimx++;

  fclose(canal);
  
  return dimx;
}


/* ===================================================================== */
void lee_serie( int dimx, char *nombre_in, double **datos) {
  /* ===================================================================== 
   *
   * Reads data stored in a text file 
   *
   * ===================================================================== */

  FILE *canal;
  char nombre[90];
  int ix;
  
  sprintf(nombre,"%s.txt",nombre_in);
  if( (canal=fopen(nombre,"rt")) == NULL) {
    fprintf(stderr,"\n Error in lee_serie while reading file %s", nombre);
    exit(-1);
  }
  for( ix=0; (ix<dimx) && (fscanf(canal,"%lf",&(datos[0][ix])) == 1); ix++ ) ;
  fclose(canal);

}

/* ===================================================================== */
int lee_dimension_datos( char *nombre_in ) {
  /* ===================================================================== 
   *
   * Returns the number of data stored in a binary file 
   *
   * ===================================================================== */
  
  FILE *canal;
  char nombre[90],muet[90];
  int dimx=0;
  double dato;

  sprintf(nombre,"%s.hdr",nombre_in);
  if( (canal=fopen(nombre,"r")) == NULL) {
    fprintf(stderr,"\n Error in lee_dimension_datos while reading header file %s"
	    "\n Read data file to get the number of entries: may be long...", nombre);
  } else {
    fscanf(canal,"%s%d",muet,&dimx);
    return dimx;
  }

  /* if we couldn't open the header file... */
  sprintf(nombre,"%s.dat",nombre_in);
  if( (canal=fopen(nombre,"rb")) == NULL) {
    fprintf(stderr,"\n Error in lee_dimension_datos while reading file %s", nombre);
    exit(-1);
  }
  while( fread(&dato,sizeof(double),1,canal) != EOF )  dimx++;
  fclose(canal);
   
  return dimx;
}


/* ===================================================================== */
void asigna_desplaza_lista(int dimx, double shift, double *funcion) {
  /* ===================================================================== */

  int ix;

  for(ix=0;ix<dimx;ix++)    funcion[ix]+=shift;

}

#ifndef FLAG_MSM

/** Functions from the MSM library but also used in the WTMM library **/

/* ===================================================================== */
void limpia(int dimx, int dimy, double **matriz) {
/* ===================================================================== */
  int ix,iy;

  for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	  matriz[iy][ix]=0.;

}

/* ===================================================================== */
void asigna_lista( int dimx, double *fuente, double *destino) {
/* ===================================================================== */
  int ix;

  for(ix=0;ix<dimx;ix++) destino[ix]=fuente[ix];
}


/* ===================================================================== */
void lee_datos( int dimx, char *nombre_in, double **signal) {
  /* ===================================================================== */
  
  FILE *canal;
  char nombre[90];
  int iy,dimy;
  
  dimy=(D_space==1)?1:dimx;

  sprintf(nombre,"%s.dat",nombre_in);
  if( (canal=fopen(nombre,"rb")) == NULL) {
    fprintf(stderr,"\n Error in lee_datos while reading file %s", nombre);
    exit(-1);
  }
  for(iy=0;iy<dimy;iy++) fread(signal[iy],sizeof(double),dimx,canal);
  fclose(canal);
  
}


/* ===================================================================== */
void graba_datos( int dimx, char *nombre_in, double **signal) {
/* ===================================================================== */

  FILE *canal;
  char nombre[90];
  int iy,dimy;
  
  dimy=(D_space==1)?1:dimx;

  sprintf(nombre,"%s.dat",nombre_in);
  if( (canal=fopen(nombre,"wb")) == NULL) {
    fprintf(stderr,"\n Error in graba_datos while reading file %s", nombre);
    exit(-1);
  }
  for(iy=0;iy<dimy;iy++) fwrite(signal[iy],sizeof(double),dimx,canal);
  fclose(canal);

  sprintf(nombre,"%s.hdr",nombre_in);
  canal=fopen(nombre,"wt");
  fprintf(canal,"%s.dat\n%d",nombre_in,dimx);
  fclose(canal);

}


/* ===================================================================== */
void graba_serie( int dimx, char *nombre_in, double *datos) {
  /* ===================================================================== */

  FILE *canal;
  char nombre[90];
  int ix;

  sprintf(nombre,"%s.txt",nombre_in);
  if( (canal=fopen(nombre,"wt")) == NULL) {
    fprintf(stderr,"\n Error in graba_serie while reading file %s", nombre);
    exit(-1);
  }
  for(ix=0;ix<dimx;ix++)
    fprintf(canal,"%f\n",datos[ix]);
  fclose(canal);

}


/* ===================================================================== */
double **reservar_matriz( int ydim, int xdim) {
/* ===================================================================== */
  double **m;
  int j;

  m=(double **) calloc(ydim,sizeof(double *));
  for(j=0;j<ydim;++j) m[j]=(double *) calloc( xdim,sizeof(double));
	
  return m;
}

/* ===================================================================== */
void liberar_matriz( double **m, int ydim) {
/* ===================================================================== */
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);

}

/* ===================================================================== */
float **reservar_matriz_float( int ydim, int xdim)
{
/* ===================================================================== */
  float **m;
  int j;

  m=(float **) calloc(ydim,sizeof(float *));
  for(j=0;j<ydim;++j) m[j]=(float *) calloc( xdim,sizeof(float));
	
  return m;
}

/* ===================================================================== */
void liberar_matriz_float( float **m, int ydim)
{
/* ===================================================================== */
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);
}

/* ===================================================================== */
int **reservar_matriz_int( int ydim, int xdim)
{
/* ===================================================================== */
  int **m;
  int j;

  m=(int **) calloc(ydim,sizeof(int *));
  for(j=0;j<ydim;++j) m[j]=(int *) calloc( xdim,sizeof(int));
	
  return m;
}

/* ===================================================================== */
void liberar_matriz_int( int **m, int ydim)
{
/* ===================================================================== */
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);
}

/* ===================================================================== */
char **reservar_matriz_char( int ydim, int xdim)
{
/* ===================================================================== */
  char **m;
  int j;

  m=(char **) calloc(ydim,sizeof(char *));
  for(j=0;j<ydim;++j) m[j]=(char *) calloc( xdim,sizeof(char));
	
  return m;
}

/* ===================================================================== */
void liberar_matriz_char( char **m, int ydim)
{
/* ===================================================================== */
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);
}


#endif
