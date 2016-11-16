/*       FFT.c. Version del 21 de Septiembre, 2005    */

#define FFT_C

#ifndef TENSOR_C
#include <tensor.c>
#endif

#ifndef OPERACIONES_C
#include <operaciones.c>
#endif

#ifndef PI
#define PI   3.1415926535997932
#endif

/*        Run-time variables    */

int MEMORY=0;  // Defaults to the use of FFT (consumes more memory, it 
               // requires the use of dimensions being power of 2 but it
               // is faster) with respecto to FFFT




/*       Function prototypes      */

int parsing_memory( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type);
int mdimensiona( int dim);
void Fourier1D(int dimx, double *funcionR, double *funcionI, int signo);
void FFT1D( int dim, double *funcionR, double *funcionI, int signo);
void FFFT1D( int dim, double *funcionR, double *funcionI, int signo);
void PFFT1D( int dima, int dimb, int dimc, double *uinR, double *uinI, 
	   double *uoutR, double *uoutI, int sign);
void convuelto_1D( int dimx, double *f1, double *f2);
void convuelve_1D( int dimx, double *f1, double *f2, double *salida);

void Fourier2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	   int signo);
void convuelve_2D( int dimx, int dimy, double **f1, double **f2, 
		   double **salida);
void convuelto_2D( int dimx, int dimy, double **f1, double **f2);
void convuelto_vec_2D( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y);
void deconvuelve_2D( int dimx, int dimy, double **f1, double **f2, 
	double **salida);
void deconvuelto_2D( int dimx, int dimy, double **f1, double **f2);
void deconvuelto_vec_2D( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y);



void FFT2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	   int signo);
void FFThorizontal(int dimx, int dimy, double **funcionhR, double **funcionhI, 
	int signo);
void FFTvertical(int dimx, int dimy, double **funcionvR, double **funcionvI, 
	int signo);
void FFFT2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	int signo);
void FFFThorizontal(int dimx, int dimy, double **funcionhR, double **funcionhI,
		    int signo);
void PFFThorizontal( int dima, int dimb, int dimc, int iy, 
		     double **uinR, double **uinI, 
		     double **uoutR, double **uoutI, int sign);
void FFFTvertical(int dimx, int dimy, double **funcionvR, double **funcionvI, 
	int signo);
void PFFTvertical( int dima, int dimb, int dimc, int ix, 
		   double **uinR, double **uinI, 
		   double **uoutR, double **uoutI, int sign);




void convuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida);
void convuelto_FFT( int dimx, int dimy, double **f1, double **f2);
void convuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y);
void deconvuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida);
void deconvuelto_FFT( int dimx, int dimy, double **f1, double **f2);
void deconvuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y);

void band_pass( int dimx, int dimy, double fmin, double fmax, 
		double **func);



void convuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida);
void convuelto_FFFT( int dimx, int dimy, double **f1, double **f2);
void convuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y);
void deconvuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida);
void deconvuelto_FFFT( int dimx, int dimy, double **f1, double **f2);
void deconvuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y);


/*      Function declarations    */


int parsing_memory( int in0, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    float **ptrvar_f, float **ptrval_f, 
		    int **ptrvar_i, int **ptrval_i,
		    int **ptrflag, int *type)
{
	int in;

	in=in0;

//    Argument MEMORY. Type 0: flag

	sprintf(olarg[in],"%s","-memory");
	sprintf(olexp[in],"%s\n %s : %s\n",
		"\nMEMORY PARAMETERS\n=================",olarg[in],
		"Flag. If enabled, the program makes a better use of memory (at the\ncost of longer processing times). Default: DISABLED");
	if(siflag)
	{
	  type[in]=3;
	  ptrvar_i[in]=&MEMORY;
	  ptrflag[in]=deflag;
	}
	else
	{
	  type[in]=0;
	  ptrflag[in]=&MEMORY;
	}
	in++;

	return(in);
}

int mdimensiona( int dim)
{
  if(MEMORY) return(dim);
  else return(dimensiona(dim));
}


void Fourier1D(int dimx, double *funcionR, double *funcionI, int signo)
{
  if(MEMORY) FFFT1D(dimx,funcionR,funcionI,signo);
  else FFT1D(dimx,funcionR,funcionI,signo);
}

void FFT1D(int dim, double *funcionR, double *funcionI, int signo)
{

	double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
	int ix,je,mm,mmax,istep;



	je=1;
	for(ix=0;ix<dim;ix++)
	{
		if(je>ix+1)
		{
			tempR=funcionR[je-1];
			tempI=funcionI[je-1];
			funcionR[je-1]=funcionR[ix];
			funcionI[je-1]=funcionI[ix];
			funcionR[ix]=tempR;
			funcionI[ix]=tempI;
		}
		mm=dim/2;

		while((mm>1)&&(je>mm))
		{
			je=je-mm;
			mm=mm/2;
		}
		je=je+mm;

	}	


	mmax=1;
	while(dim>mmax)
	{
		istep=2*mmax;
		wpasoR=cos(PI/((double) mmax));
		wpasoI=signo*sin(PI/((double) mmax));
		wwR=1.;
		wwI=0.;

		for(mm=1;mm<=mmax;mm++)
		{
		for(ix=mm-1;ix<dim;ix+=istep)
		{
			je=ix+mmax;
			C_mult(wwR,wwI,funcionR[je],funcionI[je],
				&tempR,&tempI);
			funcionR[je]=funcionR[ix]-tempR;
			funcionI[je]=funcionI[ix]-tempI;
			funcionR[ix]=funcionR[ix]+tempR;
			funcionI[ix]=funcionI[ix]+tempI;
		}
		C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
		}
		mmax=istep;

	}

	for(ix=0;ix<dim;ix++)
	{
		funcionR[ix]=funcionR[ix]/sqrt((double)dim);
 		funcionI[ix]=funcionI[ix]/sqrt((double)dim);
	}



}

void FFFT1D( int dim, double *funcionR, double *funcionI, int signo)
{
	double *workR,*workI;

	int inu;
	int dima,dimb,dimc;
	int i,i1,i2,i3;


	workR=(double *) calloc(dim,sizeof(double));
	workI=(double *) calloc(dim,sizeof(double));

	dima=1;
	dimb=dim;
	dimc=1;

	inu=1;

	while(dimb>1)
	{
	        dima=dimc*dima;
		dimc=2;
		while(dimb%dimc!=0) dimc++;

		dimb=dimb/dimc;

		if(inu==1) PFFT1D(dima,dimb,dimc,funcionR,funcionI,
			     workR,workI,signo);
		else PFFT1D(dima,dimb,dimc,workR,workI,funcionR,funcionI,
			  signo);
		inu=1-inu;
	}
	if(inu==0)
	{
	        for(i=0;i<dim;i++)
		{
		       funcionR[i]=workR[i];
		       funcionI[i]=workI[i];
		}
	}


	for(i=0;i<dim;i++)
	{
	       funcionR[i]=funcionR[i]/sqrt((double)dim);
	       funcionI[i]=funcionI[i]/sqrt((double)dim);
	}

	free(workR);
	free(workI);

}

void PFFT1D( int dima, int dimb, int dimc, double *uinR, double *uinI, 
	   double *uoutR, double *uoutI, int sign)
{
        double angle;
	double deltaR,deltaI;
	double omegaR,omegaI;
	double sumR,sumI;

	int ia,ib,ic,i,jcr,jc;

	angle=2.*PI/((double)(dima*dimc));
	omegaR=1.;
	omegaI=0.;

	deltaR=cos(angle);
	deltaI=((double)sign)*sin(angle);


	for(ic=0;ic<dimc;ic++)
	{
	for(ia=0;ia<dima;ia++)
	{
	       for(ib=0;ib<dimb;ib++)
	       {
		      i=ib+dimb*(dimc-1+ia*dimc);
		      sumR=uinR[i];
		      sumI=uinI[i];
		      for(jcr=1;jcr<dimc;jcr++)
		      {
			jc=dimc-1-jcr;
			i=ib+dimb*(jc+ia*dimc);
			C_mult(omegaR,omegaI,sumR,sumI,&sumR,&sumI);
			sumR+=uinR[i];
			sumI+=uinI[i];
		      }
		      i=ib+dimb*(ia+ic*dima);
		      uoutR[i]=sumR;
		      uoutI[i]=sumI;
	       }
	       C_mult(omegaR,omegaI,deltaR,deltaI,&omegaR,&omegaI);
	}
	}


}

void convuelto_1D( int dimx, double *f1, double *f2)
{
    double *Rf1,*Rf2,*If1,*If2;
    double buffR,buffI;
    int xeff;
    int ix;
    
    xeff=mdimensiona(dimx);
    Rf1=(double *)calloc(xeff,sizeof(double));
    Rf2=(double *)calloc(xeff,sizeof(double));
    If1=(double *)calloc(xeff,sizeof(double));
    If2=(double *)calloc(xeff,sizeof(double));
    
    asigna_lista(dimx,f1,Rf1);
    asigna_lista(dimx,f2,Rf2);

    Fourier1D(xeff,Rf2,If2,-1);
    Fourier1D(xeff,Rf1,If1,-1);
    for(ix=0;ix<xeff;ix++)
    {
	buffR=Rf1[ix]*Rf2[ix]-If1[ix]*If2[ix];
	buffI=Rf1[ix]*If2[ix]+If1[ix]*Rf2[ix];
	Rf2[ix]=buffR*sqrt(dimx);
	If2[ix]=buffI*sqrt(dimx);
    }
    Fourier1D(xeff,Rf1,If1,1);
    Fourier1D(xeff,Rf2,If2,1);

    asigna_lista(dimx,Rf2,f2);

    free(Rf1);
    free(Rf2);
    free(If1);
    free(If2);
}

void convuelve_1D( int dimx, double *f1, double *f2, double *salida)
{
    double *Rf1,*Rf2,*If1,*If2;
    double buffR,buffI;
    int xeff;
    int ix;
    
    xeff=mdimensiona(dimx);
    Rf1=(double *)calloc(xeff,sizeof(double));
    Rf2=(double *)calloc(xeff,sizeof(double));
    If1=(double *)calloc(xeff,sizeof(double));
    If2=(double *)calloc(xeff,sizeof(double));
    
    asigna_lista(dimx,f1,Rf1);
    asigna_lista(dimx,f2,Rf2);

    Fourier1D(xeff,Rf2,If2,-1);
    Fourier1D(xeff,Rf1,If1,-1);
    for(ix=0;ix<xeff;ix++)
    {
	buffR=Rf1[ix]*Rf2[ix]-If1[ix]*If2[ix];
	buffI=Rf1[ix]*If2[ix]+If1[ix]*Rf2[ix];
	Rf2[ix]=buffR*sqrt(dimx);
	If2[ix]=buffI*sqrt(dimx);
    }
    Fourier1D(xeff,Rf1,If1,1);
    Fourier1D(xeff,Rf2,If2,1);

    asigna_lista(dimx,Rf2,salida);

    free(Rf1);
    free(Rf2);
    free(If1);
    free(If2);
}


void Fourier2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	       int signo)
{
  if(MEMORY) FFFT2D(dimx,dimy,funcionR,funcionI,signo);
  else FFT2D(dimx,dimy,funcionR,funcionI,signo);
}

void convuelve_2D( int dimx, int dimy, double **f1, double **f2, 
		   double **salida)
{
  if(MEMORY) convuelve_FFFT(dimx,dimy,f1,f2,salida);
  else  convuelve_FFT(dimx,dimy,f1,f2,salida);
}

void convuelto_2D( int dimx, int dimy, double **f1, double **f2)
{
  if(MEMORY) convuelto_FFFT(dimx,dimy,f1,f2);
  else  convuelto_FFT(dimx,dimy,f1,f2);
}

void convuelto_vec_2D( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y)
{
  if(MEMORY) convuelto_vec_FFFT(dimx,dimy,f1x,f1y,f2x,f2y);
  else  convuelto_vec_FFT(dimx,dimy,f1x,f1y,f2x,f2y);
}

void deconvuelve_2D( int dimx, int dimy, double **f1, double **f2, 
	double **salida)
{
  if(MEMORY) deconvuelve_FFFT(dimx,dimy,f1,f2,salida);
  else  deconvuelve_FFT(dimx,dimy,f1,f2,salida);
}

void deconvuelto_2D( int dimx, int dimy, double **f1, double **f2)
{
  if(MEMORY) deconvuelto_FFFT(dimx,dimy,f1,f2);
  else  deconvuelto_FFT(dimx,dimy,f1,f2);
}

void deconvuelto_vec_2D( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y)
{
  if(MEMORY) deconvuelto_vec_FFFT(dimx,dimy,f1x,f1y,f2x,f2y);
  else  deconvuelto_vec_FFT(dimx,dimy,f1x,f1y,f2x,f2y);
}



void FFT2D(int dimx, int dimy, double **funcionR, double **funcionI, int signo)
{
	FFThorizontal(dimx,dimy,funcionR,funcionI,signo);
	FFTvertical(dimx,dimy,funcionR,funcionI,signo);
}

void FFThorizontal(int dimx, int dimy, double **funcionhR, double **funcionhI,
		   int signo)
{

	double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
	int iy;
	int ix,je,mm,mmax,istep;

	for(iy=0;iy<dimy;iy++)
	{

/*	INICIO DE LINEA HORIZONTAL		*/


	je=1;
	for(ix=0;ix<dimx;ix++)
	{
		if(je>ix+1)
		{
			tempR=funcionhR[iy][je-1];
			tempI=funcionhI[iy][je-1];
			funcionhR[iy][je-1]=funcionhR[iy][ix];
			funcionhI[iy][je-1]=funcionhI[iy][ix];
			funcionhR[iy][ix]=tempR;
			funcionhI[iy][ix]=tempI;
		}
		mm=dimx/2;

		while((mm>1)&&(je>mm))
		{
			je=je-mm;
			mm=mm/2;
		}
		je=je+mm;

	}	


	mmax=1;
	while(dimx>mmax)
	{
		istep=2*mmax;
		wpasoR=cos(PI/((double) mmax));
		wpasoI=signo*sin(PI/((double) mmax));
		wwR=1.;
		wwI=0.;

		for(mm=1;mm<=mmax;mm++)
		{
		for(ix=mm-1;ix<dimx;ix+=istep)
		{
			je=ix+mmax;
			C_mult(wwR,wwI,funcionhR[iy][je],funcionhI[iy][je],
				&tempR,&tempI);
			funcionhR[iy][je]=funcionhR[iy][ix]-tempR;
			funcionhI[iy][je]=funcionhI[iy][ix]-tempI;
			funcionhR[iy][ix]=funcionhR[iy][ix]+tempR;
			funcionhI[iy][ix]=funcionhI[iy][ix]+tempI;
		}
		C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
		}
		mmax=istep;

	}

	for(ix=0;ix<dimx;ix++)
	{
		funcionhR[iy][ix]=funcionhR[iy][ix]/sqrt((double)dimx);
		funcionhI[iy][ix]=funcionhI[iy][ix]/sqrt((double)dimx);
	}


/*	FIN DE LINEA HORIZONTAL			*/

	}
	


}

void FFTvertical(int dimx, int dimy, double **funcionvR, double **funcionvI, int signo)
{

	double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
	int ix;
	int iy,je,mm,mmax,istep;

	for(ix=0;ix<dimx;ix++)
	{

/*	INICIO DE LINEA VERTICAL		*/


	je=1;
	for(iy=0;iy<dimy;iy++)
	{
		if(je>iy+1)
		{
			tempR=funcionvR[je-1][ix];
			tempI=funcionvI[je-1][ix];
			funcionvR[je-1][ix]=funcionvR[iy][ix];
			funcionvI[je-1][ix]=funcionvI[iy][ix];
			funcionvR[iy][ix]=tempR;
			funcionvI[iy][ix]=tempI;
		}
		mm=dimy/2;

		while((mm>1)&&(je>mm))
		{
			je=je-mm;
			mm=mm/2;
		}
		je=je+mm;

	}	


	mmax=1;
	while(dimy>mmax)
	{
		istep=2*mmax;
		wpasoR=cos(PI/((double) mmax));
		wpasoI=signo*sin(PI/((double) mmax));
		wwR=1.;
		wwI=0.;

		for(mm=1;mm<=mmax;mm++)
		{
		for(iy=mm-1;iy<dimy;iy+=istep)
		{
			je=iy+mmax;
			C_mult(wwR,wwI,funcionvR[je][ix],funcionvI[je][ix],
				&tempR,&tempI);
			funcionvR[je][ix]=funcionvR[iy][ix]-tempR;
			funcionvI[je][ix]=funcionvI[iy][ix]-tempI;
			funcionvR[iy][ix]=funcionvR[iy][ix]+tempR;
			funcionvI[iy][ix]=funcionvI[iy][ix]+tempI;
		}
		C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
		}
		mmax=istep;

	}

	for(iy=0;iy<dimy;iy++)
	{
		funcionvR[iy][ix]=funcionvR[iy][ix]/sqrt((double)dimy);
		funcionvI[iy][ix]=funcionvI[iy][ix]/sqrt((double)dimy);
	}


/*	FIN DE LINEA VERTICAL			*/

	}
	


}

void FFT(int dim, double *funcionR, double *funcionI, int signo)
{

	double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
	int ix,je,mm,mmax,istep;



	je=1;
	for(ix=0;ix<dim;ix++)
	{
		if(je>ix+1)
		{
			tempR=funcionR[je-1];
			tempI=funcionI[je-1];
			funcionR[je-1]=funcionR[ix];
			funcionI[je-1]=funcionI[ix];
			funcionR[ix]=tempR;
			funcionI[ix]=tempI;
		}
		mm=dim/2;

		while((mm>1)&&(je>mm))
		{
			je=je-mm;
			mm=mm/2;
		}
		je=je+mm;

	}	


	mmax=1;
	while(dim>mmax)
	{
		istep=2*mmax;
		wpasoR=cos(PI/((double) mmax));
		wpasoI=signo*sin(PI/((double) mmax));
		wwR=1.;
		wwI=0.;

		for(mm=1;mm<=mmax;mm++)
		{
		for(ix=mm-1;ix<dim;ix+=istep)
		{
			je=ix+mmax;
			C_mult(wwR,wwI,funcionR[je],funcionI[je],
				&tempR,&tempI);
			funcionR[je]=funcionR[ix]-tempR;
			funcionI[je]=funcionI[ix]-tempI;
			funcionR[ix]=funcionR[ix]+tempR;
			funcionI[ix]=funcionI[ix]+tempI;
		}
		C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
		}
		mmax=istep;

	}

	for(ix=0;ix<dim;ix++)
	{
		funcionR[ix]=funcionR[ix]/sqrt((double)dim);
		funcionI[ix]=funcionI[ix]/sqrt((double)dim);
	}



}



void FFFT2D( int dimx, int dimy, double **funcionR, double **funcionI, int signo)
{
  	FFFThorizontal(dimx,dimy,funcionR,funcionI,signo);
	FFFTvertical(dimx,dimy,funcionR,funcionI,signo);
}

void FFFThorizontal( int dimx, int dimy,  double **funcionR, double **funcionI, 
		     int signo)
{
	double **workR,**workI;

	int inu;
	int dima,dimb,dimc;
	int ix,iy,i1,i2,i3;


	workR=reservar_matriz(dimy,dimx);
	workI=reservar_matriz(dimy,dimx);

	for(iy=0;iy<dimy;iy++)
	{

	  dima=1;
	  dimb=dimx;
	  dimc=1;

	  inu=1;

	  while(dimb>1)
	  {
	    dima=dimc*dima;
	    dimc=2;
	    while(dimb%dimc!=0) dimc++;

	    dimb=dimb/dimc;

	    if(inu==1) PFFThorizontal(dima,dimb,dimc,iy,funcionR,funcionI,
			     workR,workI,signo);
	    else PFFThorizontal(dima,dimb,dimc,iy,workR,workI,funcionR,funcionI,
				signo);
	    inu=1-inu;
	  }
	  if(inu==0)
	  {
	    for(ix=0;ix<dimx;ix++)
	    {
	      funcionR[iy][ix]=workR[iy][ix];
	      funcionI[iy][ix]=workI[iy][ix];
	    }
	  }


	  for(ix=0;ix<dimx;ix++)
	  {
	    funcionR[iy][ix]=funcionR[iy][ix]/sqrt((double)dimx);
	    funcionI[iy][ix]=funcionI[iy][ix]/sqrt((double)dimx);
	  }

	}

	liberar_matriz(workR,dimy);
	liberar_matriz(workI,dimy);

}

void PFFThorizontal( int dima, int dimb, int dimc, int iy, 
		     double **uinR, double **uinI, 
		     double **uoutR, double **uoutI, int sign)
{
        double angle;
	double deltaR,deltaI;
	double omegaR,omegaI;
	double sumR,sumI;

	int ia,ib,ic,ix,jcr,jc;

	angle=2.*PI/((double)(dima*dimc));
	omegaR=1.;
	omegaI=0.;

	deltaR=cos(angle);
	deltaI=((double)sign)*sin(angle);


	for(ic=0;ic<dimc;ic++)
	{
	for(ia=0;ia<dima;ia++)
	{
	       for(ib=0;ib<dimb;ib++)
	       {
		      sumR=0.;
		      sumI=0.;
		      for(jcr=0;jcr<dimc;jcr++)
		      {
			jc=dimc-1-jcr;
			ix=ib+dimb*(jc+ia*dimc);
			C_mult(omegaR,omegaI,sumR,sumI,&sumR,&sumI);
			sumR+=uinR[iy][ix];
			sumI+=uinI[iy][ix];
		      }
		      ix=ib+dimb*(ia+ic*dima);
		      uoutR[iy][ix]=sumR;
		      uoutI[iy][ix]=sumI;
	       }
	       C_mult(omegaR,omegaI,deltaR,deltaI,&omegaR,&omegaI);
	}
	}

}

void FFFTvertical( int dimx, int dimy,  double **funcionR, double **funcionI, 
		   int signo)
{
	double **workR,**workI;

	int inu;
	int dima,dimb,dimc;
	int ix,iy,i1,i2,i3;


	workR=reservar_matriz(dimy,dimx);
	workI=reservar_matriz(dimy,dimx);

	for(ix=0;ix<dimx;ix++)
	{
	  dima=1;
	  dimb=dimy;
	  dimc=1;

	  inu=1;

	  while(dimb>1)
	  {
	    dima=dimc*dima;
	    dimc=2;
	    while(dimb%dimc!=0) dimc++;

	    dimb=dimb/dimc;

	    if(inu==1) PFFTvertical(dima,dimb,dimc,ix,funcionR,funcionI,
			    workR,workI,signo);
	    else PFFTvertical(dima,dimb,dimc,ix,workR,workI,funcionR,funcionI,
		      signo);
	    inu=1-inu;
	  }
	  if(inu==0)
	  {
	    for(iy=0;iy<dimy;iy++)
	    {
	      funcionR[iy][ix]=workR[iy][ix];
	      funcionI[iy][ix]=workI[iy][ix];
	    }
	  }


	  for(iy=0;iy<dimy;iy++)
	  {
	    funcionR[iy][ix]=funcionR[iy][ix]/sqrt((double)dimy);
	    funcionI[iy][ix]=funcionI[iy][ix]/sqrt((double)dimy);
	  }

	}

	liberar_matriz(workR,dimy);
	liberar_matriz(workI,dimy);

}

void PFFTvertical( int dima, int dimb, int dimc, int ix,
		   double **uinR, double **uinI, 
		   double **uoutR, double **uoutI, int sign)
{
        double angle;
	double deltaR,deltaI;
	double omegaR,omegaI;
	double sumR,sumI;

	int ia,ib,ic,iy,jcr,jc;

	angle=2.*PI/((double)(dima*dimc));
	omegaR=1.;
	omegaI=0.;

	deltaR=cos(angle);
	deltaI=((double)sign)*sin(angle);


	for(ic=0;ic<dimc;ic++)
	{
	for(ia=0;ia<dima;ia++)
	{
	       for(ib=0;ib<dimb;ib++)
	       {
		      sumR=0.;
		      sumI=0.;
		      for(jcr=0;jcr<dimc;jcr++)
		      {
			jc=dimc-1-jcr;
			iy=ib+dimb*(jc+ia*dimc);
			C_mult(omegaR,omegaI,sumR,sumI,&sumR,&sumI);
			sumR+=uinR[iy][ix];
			sumI+=uinI[iy][ix];
		      }
		      iy=ib+dimb*(ia+ic*dima);
		      uoutR[iy][ix]=sumR;
		      uoutI[iy][ix]=sumI;
	       }
	       C_mult(omegaR,omegaI,deltaR,deltaI,&omegaR,&omegaI);
	}
	}

}





void convuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
		    double **salida)
{
	double **If1,**If2,**Isalida;
	double norma;
	int ix,iy;

	norma=sqrt((double)(dimx*dimy));
	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	Isalida=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);
	limpia(dimx,dimy,Isalida);

	FFT2D(dimx,dimy,f1,If1,-1);
	FFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		salida[iy][ix]=norma*(f1[iy][ix]*f2[iy][ix]
			-If1[iy][ix]*If2[iy][ix]);
		Isalida[iy][ix]=norma*(f1[iy][ix]*If2[iy][ix]
			+If1[iy][ix]*f2[iy][ix]);
	}
	}

	FFT2D(dimx,dimy,salida,Isalida,1);
	FFT2D(dimx,dimy,f1,If1,1);
	FFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);
	liberar_matriz(Isalida,dimy);

}

void convuelto_FFT( int dimx, int dimy, double **f1, double **f2)
{
	int ix,iy;
	double **If1,**If2;
	double buffR,buffI;
	double norma;

	norma=sqrt((double)(dimx*dimy));
	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);

	FFT2D(dimx,dimy,f1,If1,-1);
	FFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		buffR=f1[iy][ix]*f2[iy][ix]
			-If1[iy][ix]*If2[iy][ix];
		buffI=f1[iy][ix]*If2[iy][ix]
			+If1[iy][ix]*f2[iy][ix];
		f2[iy][ix]=norma*buffR;
		If2[iy][ix]=norma*buffI;
	}
	}

	FFT2D(dimx,dimy,f1,If1,1);
	FFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);

}

void convuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y)
{
	int ix,iy;
	double **If1x,**If1y,**If2x,**If2y;
	double buffRx,buffIx;
	double buffRy,buffIy;
	double norma;

	norma=sqrt((double)(dimx*dimy));

	If1x=reservar_matriz(dimy,dimx);
	If1y=reservar_matriz(dimy,dimx);
	If2x=reservar_matriz(dimy,dimx);
	If2y=reservar_matriz(dimy,dimx);

	limpia(dimx,dimy,If1x);
	limpia(dimx,dimy,If1y);
	limpia(dimx,dimy,If2x);
	limpia(dimx,dimy,If2y);

	FFT2D(dimx,dimy,f1x,If1x,-1);
	FFT2D(dimx,dimy,f1y,If1y,-1);
	FFT2D(dimx,dimy,f2x,If2y,-1);
	FFT2D(dimx,dimy,f2y,If2y,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		buffRx=f1x[iy][ix]*f2x[iy][ix]-If1x[iy][ix]*If2x[iy][ix]
		    -(f1y[iy][ix]*f2y[iy][ix]-If1y[iy][ix]*If2y[iy][ix]);
		buffIx=f1x[iy][ix]*If2x[iy][ix]+If1x[iy][ix]*f2x[iy][ix]
		    -(f1y[iy][ix]*If2y[iy][ix]+If1y[iy][ix]*f2y[iy][ix]);

		buffRy=f1x[iy][ix]*f2y[iy][ix]-If1x[iy][ix]*If2y[iy][ix]
		    +f1y[iy][ix]*f2x[iy][ix]-If1y[iy][ix]*If2x[iy][ix];
		buffIy=f1x[iy][ix]*If2y[iy][ix]+If1x[iy][ix]*f2y[iy][ix]
		    +f1y[iy][ix]*If2x[iy][ix]+If1y[iy][ix]*f2x[iy][ix];


		f2x[iy][ix]=norma*buffRx;
		f2y[iy][ix]=norma*buffRy;
		If2x[iy][ix]=norma*buffIx;
		If2y[iy][ix]=norma*buffIy;
	}
	}

	FFT2D(dimx,dimy,f1x,If1y,1);
	FFT2D(dimx,dimy,f1y,If1y,1);
	FFT2D(dimx,dimy,f2x,If2x,1);
	FFT2D(dimx,dimy,f2y,If2y,1);

	liberar_matriz(If1x,dimy);
	liberar_matriz(If1y,dimy);
	liberar_matriz(If2x,dimy);
	liberar_matriz(If2y,dimy);

}

void deconvuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida)
{
	int ix,iy;
	double **If1,**If2,**Isalida;
	double mod;

	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	Isalida=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);

	FFT2D(dimx,dimy,f1,If1,-1);
	FFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod=f2[iy][ix]*f2[iy][ix]+If2[iy][ix]*If2[iy][ix];
		if(mod>1e-30)
		{
			salida[iy][ix]=(f1[iy][ix]*f2[iy][ix]
				+If1[iy][ix]*If2[iy][ix])/mod;
			Isalida[iy][ix]=(-f1[iy][ix]*If2[iy][ix]
				+If1[iy][ix]*f2[iy][ix])/mod;
		}
		else
		{
			salida[iy][ix]=0.;
			Isalida[iy][ix]=0.;
		}
	}
	}

	FFT2D(dimx,dimy,salida,Isalida,1);
	FFT2D(dimx,dimy,f1,If1,1);
	FFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);
	liberar_matriz(Isalida,dimy);

}

void deconvuelto_FFT( int dimx, int dimy, double **f1, double **f2)
{
	int ix,iy;
	double **If1,**If2;
	double buffR,buffI,mod;

	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);

	FFT2D(dimx,dimy,f1,If1,-1);
	FFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod=f1[iy][ix]*f1[iy][ix]+If1[iy][ix]*If1[iy][ix];
		if(mod>1e-30)
		{
			buffR=(f1[iy][ix]*f2[iy][ix]
				+If1[iy][ix]*If2[iy][ix])/mod;
			buffI=(f1[iy][ix]*If2[iy][ix]
				-If1[iy][ix]*f2[iy][ix])/mod;
			f2[iy][ix]=buffR;
			If2[iy][ix]=buffI;
		}
		else
		{
			f2[iy][ix]=0.;
			If2[iy][ix]=0.;
		}
	}
	}

	FFT2D(dimx,dimy,f1,If1,1);
	FFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);

}

void deconvuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y)
{
	int ix,iy;
	double **If1x,**If1y,**If2x,**If2y;
	double buffRx,buffIx;
	double buffRy,buffIy;
	double modR,modI,modd;
	double norma;

	norma=sqrt((double)(dimx*dimy));

	If1x=reservar_matriz(dimy,dimx);
	If1y=reservar_matriz(dimy,dimx);
	If2x=reservar_matriz(dimy,dimx);
	If2y=reservar_matriz(dimy,dimx);

	limpia(dimx,dimy,If1x);
	limpia(dimx,dimy,If1y);
	limpia(dimx,dimy,If2x);
	limpia(dimx,dimy,If2y);

	FFT2D(dimx,dimy,f1x,If1x,-1);
	FFT2D(dimx,dimy,f1y,If1y,-1);
	FFT2D(dimx,dimy,f2x,If2y,-1);
	FFT2D(dimx,dimy,f2y,If2y,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{

		modR=(f1x[iy][ix]*f1x[iy][ix]
			-If1x[iy][ix]*If1x[iy][ix]
			+f1y[iy][ix]*f1y[iy][ix]
			-If1y[iy][ix]*If1y[iy][ix])*norma*norma;
		modI=(2.*f1x[iy][ix]*If1x[iy][ix]
			+2.*f1y[iy][ix]*If1y[iy][ix])*norma*norma;

		modd=modR*modR+modI*modI;

		buffRx=(f1x[iy][ix]*f2x[iy][ix]
			-If1x[iy][ix]*If2x[iy][ix]
			+f1y[iy][ix]*f2y[iy][ix]
			-If1y[iy][ix]*If2y[iy][ix])*norma;
		buffIx=(f1x[iy][ix]*If2x[iy][ix]
			+If1x[iy][ix]*f2x[iy][ix]
			+f1y[iy][ix]*If2y[iy][ix]
			+If1y[iy][ix]*f2y[iy][ix])*norma;

		buffRy=(f1x[iy][ix]*f2y[iy][ix]
			-If1x[iy][ix]*If2y[iy][ix]
		    	-f1y[iy][ix]*f2x[iy][ix]
			+If1y[iy][ix]*If2x[iy][ix])*norma;
		buffIy=(f1x[iy][ix]*If2y[iy][ix]
			+If1x[iy][ix]*f2y[iy][ix]
			-f1y[iy][ix]*If2x[iy][ix]
			-If1y[iy][ix]*f2x[iy][ix])*norma;

		if(modd>1e-30)
		{
			f2x[iy][ix]=(buffRx*modR+buffIx*modI)/modd;
			If2x[iy][ix]=(buffIx*modR-buffRx*modI)/modd;

			f2y[iy][ix]=(buffRy*modR+buffIy*modI)/modd;
			If2y[iy][ix]=(buffIy*modR-buffRy*modI)/modd;

		}
		else
		{
			f2x[iy][ix]=0.;
			If2x[iy][ix]=0.;
			f2y[iy][ix]=0.;
			If2y[iy][ix]=0.;
		}

	}
	}

	FFT2D(dimx,dimy,f1x,If1y,1);
	FFT2D(dimx,dimy,f1y,If1y,1);
	FFT2D(dimx,dimy,f2x,If2x,1);
	FFT2D(dimx,dimy,f2y,If2y,1);

	liberar_matriz(If1x,dimy);
	liberar_matriz(If1y,dimy);
	liberar_matriz(If2x,dimy);
	liberar_matriz(If2y,dimy);

}


void band_pass( int dimx, int dimy, double fmin, double fmax, 
	double **func)
{
	double **fR,**fI;
	double x,y,f;
	int xeff,yeff;
	int ix,iy;

	xeff=dimensiona(dimx);
	yeff=dimensiona(dimy);

	fR=reservar_matriz(yeff,xeff);
	fI=reservar_matriz(yeff,xeff);

	limpia(dimx,dimy,fR);
	limpia(dimx,dimy,fI);

	asigna(dimx,dimy,func,fR);
	FFT2D(xeff,yeff,fR,fI,-1);
	for(iy=0;iy<yeff;iy++)
	{
		y=((double)iy)/((double)yeff);
		if(iy>=yeff/2) y-=1.;
		for(ix=0;ix<xeff;ix++)
		{
			x=((double)ix)/((double)xeff);
			if(ix>=xeff/2) x-=1.;
			f=sqrt(x*x+y*y);
			if((f<fmin)||(f>fmax))
			{
				fR[iy][ix]=0.;
				fI[iy][ix]=0.;
			}

		}
	}
	FFT2D(xeff,yeff,fR,fI,1);
	asigna(dimx,dimy,fR,func);


	liberar_matriz(fR,yeff);
	liberar_matriz(fI,yeff);
}


void convuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida)
{
	double **If1,**If2,**Isalida;
	double norma;
	int ix,iy;

	norma=sqrt((double)(dimx*dimy));
	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	Isalida=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);
	limpia(dimx,dimy,Isalida);

	FFFT2D(dimx,dimy,f1,If1,-1);
	FFFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		salida[iy][ix]=norma*(f1[iy][ix]*f2[iy][ix]
			-If1[iy][ix]*If2[iy][ix]);
		Isalida[iy][ix]=norma*(f1[iy][ix]*If2[iy][ix]
			+If1[iy][ix]*f2[iy][ix]);
	}
	}

	FFFT2D(dimx,dimy,salida,Isalida,1);
	FFFT2D(dimx,dimy,f1,If1,1);
	FFFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);
	liberar_matriz(Isalida,dimy);

}

void convuelto_FFFT( int dimx, int dimy, double **f1, double **f2)
{
	int ix,iy;
	double **If1,**If2;
	double buffR,buffI;
	double norma;

	norma=sqrt((double)(dimx*dimy));
	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);

	FFFT2D(dimx,dimy,f1,If1,-1);
	FFFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		buffR=f1[iy][ix]*f2[iy][ix]
			-If1[iy][ix]*If2[iy][ix];
		buffI=f1[iy][ix]*If2[iy][ix]
			+If1[iy][ix]*f2[iy][ix];
		f2[iy][ix]=norma*buffR;
		If2[iy][ix]=norma*buffI;
	}
	}

	FFFT2D(dimx,dimy,f1,If1,1);
	FFFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);

}

void convuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y)
{
	int ix,iy;
	double **If1x,**If1y,**If2x,**If2y;
	double buffRx,buffIx;
	double buffRy,buffIy;
	double norma;

	norma=sqrt((double)(dimx*dimy));

	If1x=reservar_matriz(dimy,dimx);
	If1y=reservar_matriz(dimy,dimx);
	If2x=reservar_matriz(dimy,dimx);
	If2y=reservar_matriz(dimy,dimx);

	limpia(dimx,dimy,If1x);
	limpia(dimx,dimy,If1y);
	limpia(dimx,dimy,If2x);
	limpia(dimx,dimy,If2y);

	FFFT2D(dimx,dimy,f1x,If1x,-1);
	FFFT2D(dimx,dimy,f1y,If1y,-1);
	FFFT2D(dimx,dimy,f2x,If2y,-1);
	FFFT2D(dimx,dimy,f2y,If2y,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		buffRx=f1x[iy][ix]*f2x[iy][ix]-If1x[iy][ix]*If2x[iy][ix]
		    -(f1y[iy][ix]*f2y[iy][ix]-If1y[iy][ix]*If2y[iy][ix]);
		buffIx=f1x[iy][ix]*If2x[iy][ix]+If1x[iy][ix]*f2x[iy][ix]
		    -(f1y[iy][ix]*If2y[iy][ix]+If1y[iy][ix]*f2y[iy][ix]);

		buffRy=f1x[iy][ix]*f2y[iy][ix]-If1x[iy][ix]*If2y[iy][ix]
		    +f1y[iy][ix]*f2x[iy][ix]-If1y[iy][ix]*If2x[iy][ix];
		buffIy=f1x[iy][ix]*If2y[iy][ix]+If1x[iy][ix]*f2y[iy][ix]
		    +f1y[iy][ix]*If2x[iy][ix]+If1y[iy][ix]*f2x[iy][ix];


		f2x[iy][ix]=norma*buffRx;
		f2y[iy][ix]=norma*buffRy;
		If2x[iy][ix]=norma*buffIx;
		If2y[iy][ix]=norma*buffIy;
	}
	}

	FFFT2D(dimx,dimy,f1x,If1y,1);
	FFFT2D(dimx,dimy,f1y,If1y,1);
	FFFT2D(dimx,dimy,f2x,If2x,1);
	FFFT2D(dimx,dimy,f2y,If2y,1);

	liberar_matriz(If1x,dimy);
	liberar_matriz(If1y,dimy);
	liberar_matriz(If2x,dimy);
	liberar_matriz(If2y,dimy);

}

void deconvuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
	double **salida)
{
	int ix,iy;
	double **If1,**If2,**Isalida;
	double mod;

	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	Isalida=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);

	FFFT2D(dimx,dimy,f1,If1,-1);
	FFFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod=f2[iy][ix]*f2[iy][ix]+If2[iy][ix]*If2[iy][ix];
		if(mod>1e-30)
		{
			salida[iy][ix]=(f1[iy][ix]*f2[iy][ix]
				+If1[iy][ix]*If2[iy][ix])/mod;
			Isalida[iy][ix]=(-f1[iy][ix]*If2[iy][ix]
				+If1[iy][ix]*f2[iy][ix])/mod;
		}
		else
		{
			salida[iy][ix]=0.;
			Isalida[iy][ix]=0.;
		}
	}
	}

	FFFT2D(dimx,dimy,salida,Isalida,1);
	FFFT2D(dimx,dimy,f1,If1,1);
	FFFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);
	liberar_matriz(Isalida,dimy);

}

void deconvuelto_FFFT( int dimx, int dimy, double **f1, double **f2)
{
	int ix,iy;
	double **If1,**If2;
	double buffR,buffI,mod;

	If1=reservar_matriz(dimy,dimx);
	If2=reservar_matriz(dimy,dimx);
	limpia(dimx,dimy,If1);
	limpia(dimx,dimy,If2);

	FFFT2D(dimx,dimy,f1,If1,-1);
	FFFT2D(dimx,dimy,f2,If2,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
		mod=f1[iy][ix]*f1[iy][ix]+If1[iy][ix]*If1[iy][ix];
		if(mod>1e-30)
		{
			buffR=(f1[iy][ix]*f2[iy][ix]
				+If1[iy][ix]*If2[iy][ix])/mod;
			buffI=(f1[iy][ix]*If2[iy][ix]
				-If1[iy][ix]*f2[iy][ix])/mod;
			f2[iy][ix]=buffR;
			If2[iy][ix]=buffI;
		}
		else
		{
			f2[iy][ix]=0.;
			If2[iy][ix]=0.;
		}
	}
	}

	FFFT2D(dimx,dimy,f1,If1,1);
	FFFT2D(dimx,dimy,f2,If2,1);

	liberar_matriz(If1,dimy);
	liberar_matriz(If2,dimy);

}

void deconvuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
	double **f2x, double **f2y)
{
	int ix,iy;
	double **If1x,**If1y,**If2x,**If2y;
	double buffRx,buffIx;
	double buffRy,buffIy;
	double modR,modI,modd;
	double norma;

	norma=sqrt((double)(dimx*dimy));

	If1x=reservar_matriz(dimy,dimx);
	If1y=reservar_matriz(dimy,dimx);
	If2x=reservar_matriz(dimy,dimx);
	If2y=reservar_matriz(dimy,dimx);

	limpia(dimx,dimy,If1x);
	limpia(dimx,dimy,If1y);
	limpia(dimx,dimy,If2x);
	limpia(dimx,dimy,If2y);

	FFFT2D(dimx,dimy,f1x,If1x,-1);
	FFFT2D(dimx,dimy,f1y,If1y,-1);
	FFFT2D(dimx,dimy,f2x,If2y,-1);
	FFFT2D(dimx,dimy,f2y,If2y,-1);

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{

		modR=(f1x[iy][ix]*f1x[iy][ix]
			-If1x[iy][ix]*If1x[iy][ix]
			+f1y[iy][ix]*f1y[iy][ix]
			-If1y[iy][ix]*If1y[iy][ix])*norma*norma;
		modI=(2.*f1x[iy][ix]*If1x[iy][ix]
			+2.*f1y[iy][ix]*If1y[iy][ix])*norma*norma;

		modd=modR*modR+modI*modI;

		buffRx=(f1x[iy][ix]*f2x[iy][ix]
			-If1x[iy][ix]*If2x[iy][ix]
			+f1y[iy][ix]*f2y[iy][ix]
			-If1y[iy][ix]*If2y[iy][ix])*norma;
		buffIx=(f1x[iy][ix]*If2x[iy][ix]
			+If1x[iy][ix]*f2x[iy][ix]
			+f1y[iy][ix]*If2y[iy][ix]
			+If1y[iy][ix]*f2y[iy][ix])*norma;

		buffRy=(f1x[iy][ix]*f2y[iy][ix]
			-If1x[iy][ix]*If2y[iy][ix]
		    	-f1y[iy][ix]*f2x[iy][ix]
			+If1y[iy][ix]*If2x[iy][ix])*norma;
		buffIy=(f1x[iy][ix]*If2y[iy][ix]
			+If1x[iy][ix]*f2y[iy][ix]
			-f1y[iy][ix]*If2x[iy][ix]
			-If1y[iy][ix]*f2x[iy][ix])*norma;

		if(modd>1e-30)
		{
			f2x[iy][ix]=(buffRx*modR+buffIx*modI)/modd;
			If2x[iy][ix]=(buffIx*modR-buffRx*modI)/modd;

			f2y[iy][ix]=(buffRy*modR+buffIy*modI)/modd;
			If2y[iy][ix]=(buffIy*modR-buffRy*modI)/modd;

		}
		else
		{
			f2x[iy][ix]=0.;
			If2x[iy][ix]=0.;
			f2y[iy][ix]=0.;
			If2y[iy][ix]=0.;
		}

	}
	}

	FFFT2D(dimx,dimy,f1x,If1y,1);
	FFFT2D(dimx,dimy,f1y,If1y,1);
	FFFT2D(dimx,dimy,f2x,If2x,1);
	FFFT2D(dimx,dimy,f2y,If2y,1);

	liberar_matriz(If1x,dimy);
	liberar_matriz(If1y,dimy);
	liberar_matriz(If2x,dimy);
	liberar_matriz(If2y,dimy);

}



