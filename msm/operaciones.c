/*	operaciones.c. Version del 24 de Agosto, 2004		*/

#define OPERACIONES_C

#ifndef PI
#define PI   3.1415926535997932
#endif

/*     Function prototypes    */

double fit(double *x, double *y, int n, double *a, double *b, double *corr);
double best_fit(double *x, double *y, int n, int hmin, 
		double **repx, double **repy, int *N);

int dimensiona( int dim);
int adimensiona( int size);

double fMax( double a, double b);
double fMin( double a, double b);
int Max( int a, int b);
int Min( int a, int b);
int Mod( int a, int b);
int Round( double a);
double TanH( double x);
void C_mult( double a1, double b1, double a2, 
	double b2, double *a, double *b );
void C_sqrt( double a0, double b0, double *a, double *b);

double cuantil( int dimx, double *data, double prob0);
void quicksort( int low, int high, double *data );
int partition( int low, int high, double *data);
void quicksort_ref( int low, int high, int *ref, double *data );
int partition_ref( int low, int high, int *ref, double *data);
double cuantil_2D( int dimx, int dimy, double **data, double prob0);
double cuantil_2D_sort( int dimx, int dimy, double **data, double prob0);
double cuantil_2D_histo( int dimx, int dimy, double **data, double prob0);


/*   Function declarations   */


double fit(double *x, double *y, int n, double *a, double *b, double *corr)
{
	double sumx,sumy,sumxx,sumyy,sumxy;
	double sx,sy,sxy;
	int i;


	sumx=0.;
	sumy=0.;
	sumxx=0.;
	sumxy=0.;
	sumyy=0.;
	for(i=0;i<n;i++)
	{
		sumx+=x[i];
		sumy+=y[i];
		sumxx+=x[i]*x[i];
		sumxy+=x[i]*y[i];
		sumyy+=y[i]*y[i];
	}

	sxy=sumxy/((double)n)-sumx*sumy/((double)(n*n));
	sx=sqrt(fabs(sumxx/((double)n)-sumx*sumx/((double)(n*n))));
	sy=sqrt(fabs(sumyy/((double)n)-sumy*sumy/((double)(n*n))));
	if(sx>1e-30)
	{
	  *a=sxy/(sx*sx);
	  *b=(sumy-*a*sumx)/((double)n);
	  if(sy<1e-30) *corr=1.;
	  else *corr=sxy/(sx*sy);
	}
	else
	{
	  *a=0.;
	  *b=sumy/((double)n);
	  *corr=1.;
	}

	return(sy*sy-a[0]*a[0]*sx*sx);
}


double best_fit(double *x, double *y, int n, int hmin, 
		double **repx, double **repy, int *N)
{
        double *histoc,*rep1,*repe2;
	double mm1[2];
	double norm;
	double m2,me2,var2,vare2,regr;
	double histac,rep1ac,repe2ac;
       	int ix,ip,ip0;


	*N=n;
	histoc=(double *) calloc(n,sizeof(double));
	rep1=(double *) calloc(n,sizeof(double));
	repe2=(double *) calloc(n,sizeof(double));


	mm1[0]=mm1[1]=x[0];
	for(ix=0;ix<n;ix++)
	{
	  mm1[0]=fMin(mm1[0],x[ix]);
	  mm1[1]=fMax(mm1[1],x[ix]);
	}
	mm1[1]-=mm1[0];

	m2=0.;
	var2=0.;
	for(ix=0;ix<n;ix++)
	{
	  ip=(int) (((double)*N)*(x[ix]-mm1[0])/mm1[1]);
	  if(ip>=*N) ip=*N-1;
	  histoc[ip]+=1.;
	  rep1[ip]+=x[ix];
	  repe2[ip]+=y[ix];
	  m2+=y[ix];
	  var2+=y[ix]*y[ix];
	}
	m2/=(double)n;
	var2=var2/((double)n)-m2*m2;

	histac=0.;
	rep1ac=0.;
	repe2ac=0.;
	ip0=0;
	for(ip=0;ip<*N;ip++)
	{
	  if(histac<hmin)
	  {
	    histac+=histoc[ip];
	    rep1ac+=rep1[ip];
	    repe2ac+=repe2[ip];
	  }
	  else
	  {
	    histoc[ip0]=histac;
	    rep1[ip0]=rep1ac;
	    repe2[ip0]=repe2ac;
	    histac=histoc[ip];
	    rep1ac=rep1[ip];
	    repe2ac=repe2[ip];
	    ip0++;
	  }
	}
	histoc[ip0]=histac;
	rep1[ip0]=rep1ac;
	repe2[ip0]=repe2ac;
	*N=ip0+1;
	
	for(ip=0;ip<*N;ip++)
	{
	    rep1[ip]/=histoc[ip];
	    repe2[ip]/=histoc[ip];
	}

	me2=vare2=norm=0.;
	for(ip=0;ip<*N;ip++)
	{
	  norm+=histoc[ip];
	  me2+=repe2[ip]*histoc[ip];
	  vare2+=repe2[ip]*repe2[ip]*histoc[ip];
	}
	me2/=norm;
	vare2=vare2/norm-me2*me2;

	if(var2>1e-30) regr=vare2/var2;
	else regr=1.;

	

	*repx=(double *) calloc(*N,sizeof(double));
	*repy=(double *) calloc(*N,sizeof(double));

	for(ip=0;ip<*N;ip++)
	{
	  repx[0][ip]=rep1[ip];
	  repy[0][ip]=repe2[ip];
	}


	free(histoc);
	free(rep1);
	free(repe2);

	return(regr);
}


int dimensiona( int dim)
{
	int out;

	out=1;
	while(out<dim) out=2*out;

	return(out);
}

int adimensiona( int size)
{
	int out,s0;

	for(out=0,s0=size;s0>1;s0=s0/2) out++;
	

	return(out);
}

double fMax(double a,double b)
{
	double salida;

	if(a>b) salida=a;
	else salida=b;

	return(salida);
}

double fMin(double a,double b)
{
	double salida;

	if(a<b) salida=a;
	else salida=b;

	return(salida);
}

int Max(int a,int b)
{
	int salida;

	if(a>b) salida=a;
	else salida=b;

	return(salida);
}

int Min(int a,int b)
{
	int salida;

	if(a<b) salida=a;
	else salida=b;

	return(salida);
}

int Mod(int a, int b)
{
	int output;

	output=a/b;
	output=a-output*b;
	if(output<0) output+=b;

	return(output);
}

int Round( double a)
{
	int out;

	out=(int)(a+0.5);
	if(a<-0.5) out--;

	return(out);
}


void C_mult( double a1, double b1, double a2, double b2, double *a, double *b )
{
	*a=a1*a2-b1*b2;
	*b=a1*b2+a2*b1;
}

void C_sqrt( double a0, double b0, double *a, double *b)
{
	double mod;

	mod=sqrt(a0*a0+b0*b0);

	*a=sqrt(0.5*(mod+a0));
	if(b0<0) *b=-sqrt(0.5*(mod-a0));
	else *b=sqrt(0.5*(mod-a0));

}

double cuantil( int dimx, double *data, double prob0)
{
	double *sort;
	double quant;
	int ix;

	sort=(double *) calloc(dimx,sizeof(double));

	for(ix=0;ix<dimx;ix++) sort[ix]=data[ix];
	quicksort(0,dimx-1,sort);

	ix=(int)(prob0*(double)(dimx-1));
	quant=sort[ix];

	free(sort);
	return(quant);
}


void quicksort( int low, int high, double *data )
{
        int pivot;

	if(low<high)
	{   
	      pivot=partition(low,high,data);
	      quicksort(low,pivot-1,data);
	      quicksort(pivot+1,high,data);
	}
}

int partition( int low, int high, double *data)
{
        double pivot_item;
	double buff;
        int left,right;

	pivot_item=data[low];
	left=low;
	right=high;
	while(left<right)
	{
	       while((data[left]<=pivot_item)&&(left<=right)) left++;
	       while((data[right]>=pivot_item)&&(left<=right)) right--;
	       if(left<right)
	       {
		      buff=data[left];
		      data[left]=data[right];
		      data[right]=buff;
	       }
	}
	data[low]=data[right];
	data[right]=pivot_item;

	return(right);

}

void quicksort_ref( int low, int high, int *ref, double *data )
{
        int pivot;

	if(low<high)
	{   
	      pivot=partition_ref(low,high,ref,data);
	      quicksort_ref(low,pivot-1,ref,data);
	      quicksort_ref(pivot+1,high,ref,data);
	}
}

int partition_ref( int low, int high, int *ref, double *data)
{
        double pivot_item;
	double buff;
	int pivot_item_i,buffi;
        int left,right;

	pivot_item=data[low];
	pivot_item_i=ref[low];
	left=low;
	right=high;
	while(left<right)
	{
	       while((data[left]<=pivot_item)&&(left<=right)) left++;
	       while((data[right]>=pivot_item)&&(left<=right)) right--;
	       if(left<right)
	       {
		      buff=data[left];
		      data[left]=data[right];
		      data[right]=buff;
		      buffi=ref[left];
		      ref[left]=ref[right];
		      ref[right]=buffi;
	       }
	}
	data[low]=data[right];
	ref[low]=ref[right];
	data[right]=pivot_item;
	ref[right]=pivot_item_i;

	return(right);

}

double cuantil_2D( int dimx, int dimy, double **data, double prob0)
{
        if(dimx*dimy<10000) return(cuantil_2D_sort(dimx,dimy,data,prob0));
	else return(cuantil_2D_histo(dimx,dimy,data,prob0));
}




double cuantil_2D_sort( int dimx, int dimy, double **data, double prob0)
{
	double *sort;
	double quant;
	int ix,iy;


	sort=(double *) calloc(dimx*dimy,sizeof(double));

	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  sort[ix+dimx*iy]=data[iy][ix];
	}
	}
	quicksort(0,dimx*dimy-1,sort);

	ix=(int)(prob0*(double)(dimx*dimy-1));
	quant=sort[ix];


	free(sort);
	return(quant);
}


double cuantil_2D_histo( int dimx, int dimy, double **data, double prob0)
{
	double *histo;
	double maxd,mind,norma;
	double quant;
	const int Nhisto=10000;
	int ix,iy,ip;
	
	histo=(double *)calloc(Nhisto+1,sizeof(double));
	norma=1./((double)(dimx*dimy));

	mind=maxd=data[0][0];
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	  mind=fMin(mind,data[iy][ix]);
	  maxd=fMax(maxd,data[iy][ix]);
	}
	}
	maxd-=mind;

	if(maxd>1e-30)
	{
	  for(iy=0;iy<dimy;iy++)
	  {
	  for(ix=0;ix<dimx;ix++)
	  {
	    ip=(int)(Nhisto*(data[iy][ix]-mind)/maxd);
	    histo[ip]+=norma;
	  }
	  }
	  norma=0.;
	  for(ip=0;(ip<=Nhisto)&&(norma<prob0);ip++) norma+=histo[ip];
	  quant=mind+maxd*((double)ip)/((double)Nhisto);
	}
	else quant=mind;


	free(histo);
	return(quant);
}
