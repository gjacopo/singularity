/*    derivacion.c.  Version del 24 de Agosto, 2004         */

#define DERIVACION_C

#ifndef FFT_C
#include <FFT.c>
#endif

/*     Run-time variables     */

int DERIVA_MODO=0; // Derivation mode.
                   // 0: Half-pixel by Fourier transform. Generates aliasing 
                   // in wide areas of low gradient
                   // 1: Forward 1-pixel increment. Simple and less sensitive 
                   // to artifacts, but reconstruction from it is of less 
                   // quality (unknown reason)




/*     Function prototypes   */

int parsing_derivacion( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type);

void gradiente_1D( int dimx, double *gx);
void reconstruye_1D( int dimx, double *gx);
void gradiente_naif_1D( int dimx, double *gx);
void reconstruye_naif_1D( int dimx, double *gx);
void gradiente_FFT_1D( int dimx, double *gx);
void reconstruye_FFT_1D( int dimx, double *gx);
void gradiente_naif_inter_1D( int dimx, double *gx);
void filtro( int dimx, double norma, double expon, double *funcion);


void gradiente_2D( int dimx, int dimy,  double **gx, double **gy);
void reconstruye_2D( int dimx, int dimy, double **gx, double **gy);
void gradiente_naif_2D_bak( int dimx, int dimy,  double **gx, double **gy);
void gradiente_naif_2D( int dimx, int dimy,  double **gx, double **gy);
void reconstruye_naif_2D( int dimx, int dimy, double **gx, double **gy);
void gradiente_FFT_2D( int dimx, int dimy,  double **gx, double **gy);
void reconstruye_FFT_2D( int dimx, int dimy, double **gx, double **gy);
void gradiente_naif_inter_2D( int dimx, int dimy,  double **gx, double **gy);
void gradiente_complex( int dimx, int dimy,  double **gx, double **gy);
void reconstruye_complex( int dimx, int dimy, double **gx, double **gy);
void deriva_complex( int dimx, int dimy, int mode, double **gx, double **gy);
void deriva_complex_naif( int dimx, int dimy, int mode, double **gx, 
			  double **gy);
void deriva_complex_FFT( int dimx, int dimy, int mode, double **gx, 
			 double **gy);
void filtro_2D( int dimx, int dimy, double expon, double **data);


/*     Function declarations       */

int parsing_derivacion( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type)
{
	int in;

	in=in0;

//    Argument DERIVA_MODO. Type 1: Integer

	sprintf(olarg[in],"%s","-dermode");
	strcpy(olval[in],"mode");
	sprintf(olexp[in],"%s\n %s : %s %d\n",
		"\nDERIVATIVE PARAMETERS\n=================",olarg[in],
		"Method for calculating derivatives\n    0: Half-pixel FFT kernel derivative\n    1: Direct 1-pixel increment\n    2: Derivative interpolation\nDefault:",DERIVA_MODO);
	if(siflag)
	{
	  type[in]=4;
	  ptrflag[in]=deflag;
	}
	else type[in]=1;
	ptrvar_i[in]=&DERIVA_MODO;
	ptrval_i[in][0]=0;
	ptrval_i[in][1]=2;
	in++;

	return(in);
}

void gradiente_1D( int dimx,  double *gx)
{
    switch(DERIVA_MODO)
    {
        case 0:
	  gradiente_FFT_1D(dimx,gx);
	  break;
        case 1:
	  gradiente_naif_1D(dimx,gx);
	  break;
        case 2:
	  gradiente_naif_inter_1D(dimx,gx);
	  break;
        default:
	  printf("No definition associated to the derivative mode\n");
	  exit(-1);
    }
}

void reconstruye_1D( int dimx, double *gx)
{
    switch(DERIVA_MODO)
    {
        case 0:
	  reconstruye_FFT_1D(dimx,gx);
	  break;
        case 1:
	  reconstruye_naif_1D(dimx,gx);
	  break;
        case 2:
	  reconstruye_naif_1D(dimx,gx);
	  break;
        default:
	  printf("No definition associated to the derivative mode\n");
	  exit(-1);
    }
}

void gradiente_naif_1D( int dimx, double *gx)
{
	double cierre;
	int ix;

	cierre=gx[0]-gx[dimx-1];
	for(ix=0;ix<dimx-1;ix++) gx[ix]=gx[ix+1]-gx[ix];

	gx[dimx-1]=cierre;
}

void reconstruye_naif_1D( int dimx, double *gx)
{
	double *aux;
	int ix;

	aux=(double *) calloc(dimx,sizeof(double));
	for(ix=0;ix<dimx;ix++) aux[ix]=gx[ix];
	gx[0]=0.;
	for(ix=1;ix<dimx;ix++)
		gx[ix]=gx[ix-1]+aux[ix-1];

	free(aux);
}

void gradiente_FFT_1D( int dimx, double *gx)
{
         double *aux;
	 double dR,x;
	 int ix;

	 aux=(double *) calloc(dimx,sizeof(double));
	 Fourier1D(dimx,gx,aux,-1);
	 for(ix=0;ix<dimx;ix++)
	 {

	   x=((double)ix)/((double)dimx);
	   if(ix>=dimx/2) x-=1.;
	   x=sin(PI*x);
	   dR=-x*aux[ix];
	   aux[ix]=x*gx[ix];
	   gx[ix]=dR;
	 }
	 Fourier1D(dimx,gx,aux,1);
	 free(aux);
}

void reconstruye_FFT_1D( int dimx, double *gx)
{
         double *aux;
	 double dR,x;
	 int ix;

	 aux=(double *) calloc(dimx,sizeof(double));
	 Fourier1D(dimx,gx,aux,-1);
	 for(ix=0;ix<dimx;ix++)
	 {

	   x=((double)ix)/((double)dimx);
	   if(ix>=dimx/2) x-=1.;
	   if(x>1e-30) x=1./sin(PI*x);
	   else x=0.;
	   dR=x*aux[ix];
	   aux[ix]=-x*gx[ix];
	   gx[ix]=dR;
	 }
	 Fourier1D(dimx,gx,aux,1);
	 free(aux);
}

void gradiente_naif_inter_1D( int dimx, double *gx)
{
         double *aux;
	 double buff;
	 int ix,ic;

	 aux=(double *) calloc(dimx,sizeof(double));
	 for(ix=0;ix<dimx;ix++) aux[ix]=gx[ix];
	 gradiente_naif_1D(dimx,aux);
	 for(ix=0;ix<dimx;ix++)
	 {
	   for(ic=1,buff=aux[ix];(ic<10)&&(fabs(buff)<1e-3);ic++)
	     buff+=aux[Mod(ix-ic,dimx)]+aux[Mod(ix+ic,dimx)];
	   buff/=(double)(2*ic-1);
	   gx[ix]=buff;
	 }


	 free(aux);
}


void filtro_1D( int dimx,double expon, double *funcion)
{
	double *funcR,*funcI;
	double x,f;
	double dx,dy;
	int ix;
	int xeff;

	xeff=mdimensiona(dimx);
	funcR=(double *) calloc(xeff,sizeof(double));
	funcI=(double *) calloc(xeff,sizeof(double));

	asigna_lista(dimx,funcion,funcR);
	Fourier1D(xeff,funcR,funcI,-1);
	for(ix=0;ix<xeff;ix++)
	{
	  x=((double)ix)/((double)xeff);
	  if(ix>=xeff/2) x-=1.;

	  dx=2.*sin(PI*x);

	  f=fabs(dx);
	  if(f>1e-30) f=pow(f,expon);
	  else f=0.;
	  
	  funcR[ix]=f*funcR[ix];
	  funcI[ix]=f*funcI[ix];
	}
	Fourier1D(xeff,funcR,funcI,1);
	asigna_lista(dimx,funcR,funcion);

	free(funcR);

}


void gradiente_2D( int dimx, int dimy,  double **gx, double **gy)
{
    switch(DERIVA_MODO)
    {
        case 0:
	  gradiente_FFT_2D(dimx,dimy,gx,gy);
	  break;
        case 1:
	  gradiente_naif_2D(dimx,dimy,gx,gy);
	  break;
        case 2:
	  gradiente_naif_inter_2D(dimx,dimy,gx,gy);
	  break;
        default:
	  printf("Unrecognized derivative mode\n");
	  exit(-1);
    }
}

void reconstruye_2D( int dimx, int dimy, double **gx, double **gy)
{
    switch(DERIVA_MODO)
    {
        case 0:
	  reconstruye_FFT_2D(dimx,dimy,gx,gy);
	  break;
        case 1:
	  reconstruye_naif_2D(dimx,dimy,gx,gy);
	  break;
        case 2:
	  reconstruye_naif_2D(dimx,dimy,gx,gy);
	  break;
        default:
	  printf("Unrecognized derivative mode\n");
	  exit(-1);
    }
}

void gradiente_naif_2D_bak( int dimx, int dimy, double **gx, double **gy)
{
         double *columna_derecha;
	 int ix,iy;

	 columna_derecha=(double *) calloc(dimy,sizeof(double));

	 for(iy=0;iy<dimy;iy++) columna_derecha[iy]=gx[iy][0]-gx[iy][dimx-1];

	 for(iy=0;iy<dimy-1;iy++)
	 {
	   for(ix=0;ix<dimx;ix++) gy[iy][ix]=gx[iy+1][ix]-gx[iy][ix];
	 }
	 for(ix=0;ix<dimx;ix++) gy[dimy-1][ix]=gx[0][ix]-gx[dimy-1][ix];

	 for(ix=0;ix<dimx-1;ix++)
	 {
	   for(iy=0;iy<dimy;iy++) gx[iy][ix]=gx[iy][ix+1]-gx[iy][ix];
	 }
	 for(iy=0;iy<dimy;iy++) gx[iy][dimx-1]=columna_derecha[iy];

	 free(columna_derecha);

}

void gradiente_naif_2D( int dimx, int dimy, double **gx, double **gy)
{
         double **gxI,**gyI;
	 double dxR,dxI,dyR,dyI;
	 double aux,x,y;
	 int ix,iy;

	 gxI=reservar_matriz(dimy,dimx);
	 gyI=reservar_matriz(dimy,dimx);

	 Fourier2D(dimx,dimy,gx,gxI,-1);

	 for(iy=0;iy<dimy;iy++)
	 {
	   y=((double)iy)/((double)dimy);
	   dyR=cos(2.*PI*y)-1.;
	   dyI=sin(2.*PI*y);
	   for(ix=0;ix<dimx;ix++)
	   {
	     x=((double)ix)/((double)dimx);
	     dxR=cos(2.*PI*x)-1.;
	     dxI=sin(2.*PI*x);

	     gy[iy][ix]=dyR*gx[iy][ix]-dyI*gxI[iy][ix];
	     gyI[iy][ix]=dyR*gxI[iy][ix]+dyI*gx[iy][ix];
	     aux=dxR*gx[iy][ix]-dxI*gxI[iy][ix];
	     gxI[iy][ix]=dxR*gxI[iy][ix]+dxI*gx[iy][ix];
	     gx[iy][ix]=aux;
	   }

	 }
	 Fourier2D(dimx,dimy,gx,gxI,1);
	 Fourier2D(dimx,dimy,gy,gyI,1);


	 liberar_matriz(gxI,dimy);
	 liberar_matriz(gyI,dimy);

}

void reconstruye_naif_2D( int dimx, int dimy, double **gx, double **gy)
{
         double **gxI,**gyI;
	 double dxR,dxI,dyR,dyI,modx,mody;
	 double auxR,auxI,x,y;
	 int ix,iy;

	 gxI=reservar_matriz(dimy,dimx);
	 gyI=reservar_matriz(dimy,dimx);

	 Fourier2D(dimx,dimy,gx,gxI,-1);
	 Fourier2D(dimx,dimy,gy,gyI,-1);

	 for(iy=0;iy<dimy;iy++)
	 {
	   y=((double)iy)/((double)dimy);
	   dyR=cos(2.*PI*y)-1.;
	   dyI=-sin(2.*PI*y);
	   mody=dyR*dyR+dyI*dyI;
	   for(ix=0;ix<dimx;ix++)
	   {
	     x=((double)ix)/((double)dimx);
	     dxR=cos(2.*PI*x)-1.;
	     dxI=-sin(2.*PI*x);
	     modx=dxR*dxR+dxI*dxI;

	     if(modx+mody>1e-30)
	     {
	       auxR=(dxR*gx[iy][ix]-dxI*gxI[iy][ix]
		 +dyR*gy[iy][ix]-dyI*gyI[iy][ix])/(modx+mody);
	       auxI=(dxR*gxI[iy][ix]+dxI*gx[iy][ix]
		 +dyR*gyI[iy][ix]+dyI*gy[iy][ix])/(modx+mody);
	     }
	     else
	     {
	       auxR=0.;
	       auxI=0.;
	     }
	     gx[iy][ix]=auxR;
	     gxI[iy][ix]=auxI;
	   }

	 }
	 Fourier2D(dimx,dimy,gx,gxI,1);


	 liberar_matriz(gxI,dimy);
	 liberar_matriz(gyI,dimy);

}

void gradiente_FFT_2D( int dimx, int dimy, double **gx, double **gy)
{
         double **gxI,**gyI;
	 double aux,x,y;
	 int ix,iy;

	 gxI=reservar_matriz(dimy,dimx);
	 gyI=reservar_matriz(dimy,dimx);

	 Fourier2D(dimx,dimy,gx,gxI,-1);

	 for(iy=0;iy<dimy;iy++)
	 {
	   y=((double)iy)/((double)dimy);
	   if((iy>0)&&(iy>=dimy/2)) y-=1.;
	   y=sin(PI*y);
	   for(ix=0;ix<dimx;ix++)
	   {
	     x=((double)ix)/((double)dimx);
	     if((ix>0)&&(ix>=dimx/2)) x-=1.;
	     x=sin(PI*x);

	     gy[iy][ix]=-y*gxI[iy][ix];
	     gyI[iy][ix]=y*gx[iy][ix];
	     aux=-x*gxI[iy][ix];
	     gxI[iy][ix]=x*gx[iy][ix];
	     gx[iy][ix]=aux;
	   }

	 }

	 Fourier2D(dimx,dimy,gx,gxI,1);
	 Fourier2D(dimx,dimy,gy,gyI,1);


	 liberar_matriz(gxI,dimy);
	 liberar_matriz(gyI,dimy);

}

void reconstruye_FFT_2D( int dimx, int dimy, double **gx, double **gy)
{
         double **gxI,**gyI;
	 double aux,x,y,f;
	 int ix,iy;

	 gxI=reservar_matriz(dimy,dimx);
	 gyI=reservar_matriz(dimy,dimx);

	 Fourier2D(dimx,dimy,gx,gxI,-1);
	 Fourier2D(dimx,dimy,gy,gyI,-1);

	 for(iy=0;iy<dimy;iy++)
	 {
	   y=((double)iy)/((double)dimy);
	   if((iy>0)&&(iy>=dimy/2)) y-=1.;
	   y=sin(PI*y);
	   for(ix=0;ix<dimx;ix++)
	   {
	     x=((double)ix)/((double)dimx);
	     if((ix>0)&&(ix>=dimx/2)) x-=1.;
	     x=sin(PI*x);

	     f=x*x+y*y;
	     if(f>1e-30)
	     {
	       aux=(x*gxI[iy][ix]+y*gyI[iy][ix])/f;
	       gxI[iy][ix]=-(x*gx[iy][ix]+y*gy[iy][ix])/f;
	       gx[iy][ix]=aux;
	     }
	     else
	     {
	       gx[iy][ix]=0.;
	       gxI[iy][ix]=0.;
	     }
	     gy[iy][ix]=0.;
	     gyI[iy][ix]=0.;
	   }

	 }

	 Fourier2D(dimx,dimy,gx,gxI,1);


	 liberar_matriz(gxI,dimy);
	 liberar_matriz(gyI,dimy);

}

void gradiente_naif_inter_2D( int dimx, int dimy, double **gx, double **gy)
{
	 double x,y;
	 int ix,iy;


	 for(iy=0;iy<dimy;iy++)
	 {
	     y=(double)iy;
	     if(iy>dimy/2) y-=(double)dimy;
	     for(ix=0;ix<dimx;ix++)
	     {
		 x=(double)ix;
		 if(ix>dimx/2) x-=(double)dimx;
		 gy[iy][ix]=1./(.00001+x*x+y*y);
	     }
	 }
	 convuelto_2D(dimx,dimy,gy,gx);
	 gradiente_naif_2D(dimx,dimy,gx,gy);


}


void gradiente_complex( int dimx, int dimy, double **gx, double **gy)
{
  limpia(dimx,dimy,gy);
  deriva_complex(dimx,dimy,0,gx,gy);
}

void reconstruye_complex( int dimx, int dimy, double **gx, double **gy)
{
  deriva_complex(dimx,dimy,1,gx,gy);
  limpia(dimx,dimy,gy);
}

void deriva_complex( int dimx, int dimy, int mode, double **gx, double **gy)
{
    switch(DERIVA_MODO)
    {
        case 0:
	  deriva_complex_FFT(dimx,dimy,mode,gx,gy);
	  break;
        case 1:
	  deriva_complex_naif(dimx,dimy,mode,gx,gy);
	  break;
         default:
	  printf("Unrecognized derivative mode\n");
	  exit(-1);
    }
}

void deriva_complex_naif( int dimx, int dimy, int mode, double **gx, 
			  double **gy)
{
	 double dxR,dxI,dyR,dyI;
	 double dR,dI,modd;
	 double aux,x,y;
	 int ix,iy;

	 Fourier2D(dimx,dimy,gx,gy,-1);

	 for(iy=0;iy<dimy;iy++)
	 {
	   y=((double)iy)/((double)dimy);
	   dyR=-(cos(2.*PI*y)-1.);
	   dyI=-sin(2.*PI*y);

	   for(ix=0;ix<dimx;ix++)
	   {
	     x=((double)ix)/((double)dimx);
	     dxR=cos(2.*PI*x)-1.;
	     dxI=sin(2.*PI*x);



	     dR=dxR+dyI;
	     dI=dxI-dyR;
	     modd=dR*dR+dI*dI;
	     if(mode)
	     {
	       if(modd>1e-30)
	       {
		 dR=dR/modd;
		 dI=-dI/modd;
	       }
	       else
	       {
		 dR=0.;
		 dI=0;
	       }
	     }

	     aux=dR*gx[iy][ix]-dI*gy[iy][ix];
	     gy[iy][ix]=dR*gy[iy][ix]+dI*gx[iy][ix];
	     gx[iy][ix]=aux;
	   }
	 }

	 Fourier2D(dimx,dimy,gx,gy,1);

}

void deriva_complex_FFT( int dimx, int dimy, int mode, double **gx, 
			 double **gy)
{
	 double dxR,dxI,dyR,dyI;
	 double dR,dI,modd;
	 double aux,x,y;
	 int ix,iy;

	 Fourier2D(dimx,dimy,gx,gy,-1);

	 for(iy=0;iy<dimy;iy++)
	 {
	   y=((double)iy)/((double)dimy);
	   if(iy>dimy/2) y-=1.;
	   dyR=0.;
	   dyI=-sin(PI*y);

	   for(ix=0;ix<dimx;ix++)
	   {
	     x=((double)ix)/((double)dimx);
	     if(ix>dimx/2) x-=1.;
	     dxR=0.;
	     dxI=sin(PI*x);



	     dR=dxR+dyI;
	     dI=dxI-dyR;
	     modd=dR*dR+dI*dI;
	     if(mode)
	     {
	       if(modd>1e-30)
	       {
		 dR=dR/modd;
		 dI=-dI/modd;
	       }
	       else
	       {
		 dR=0.;
		 dI=0;
	       }
	     }

	     aux=dR*gx[iy][ix]-dI*gy[iy][ix];
	     gy[iy][ix]=dR*gy[iy][ix]+dI*gx[iy][ix];
	     gx[iy][ix]=aux;
	   }
	 }

	 Fourier2D(dimx,dimy,gx,gy,1);

}

void filtro_2D( int dimx, int dimy, double expon, double **data)
{
        double **auxR,**auxI;
	double x,y,f;
	int xeff,yeff;
	int ix,iy;

	xeff=mdimensiona(dimx);
	yeff=mdimensiona(dimy);


	auxR=reservar_matriz(yeff,xeff);
	auxI=reservar_matriz(yeff,xeff);

	asigna(dimx,dimy,data,auxR);

	Fourier2D(xeff,yeff,auxR,auxI,-1);

	for(iy=0;iy<yeff;iy++)
	{
	   y=((double)iy)/((double)yeff);
	   if(iy>=yeff/2) y-=1.;
	   y=sin(PI*y);
	   for(ix=0;ix<xeff;ix++)
	   {
	     x=((double)ix)/((double)xeff);
	     if(ix>=xeff/2) x-=1.;
	     x=sin(PI*x);
	     f=sqrt(x*x+y*y);
	     if(f>1e-30) f=pow(f,expon);
	     else f=0.;
	     auxR[iy][ix]*=f;
	     auxI[iy][ix]*=f;
	  }
	}
	Fourier2D(xeff,yeff,auxR,auxI,1);

	asigna(dimx,dimy,auxR,data);


	liberar_matriz(auxR,yeff);
	liberar_matriz(auxI,yeff);

}

