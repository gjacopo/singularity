/*	graficos.c.  Version del 14 de Septiembre de 2004		*/


/*   
      Comment: Starting in September, 14th, matrices will be assumed to be
      referred to physical axis, that is, the (0,0) corresponds to the 
      left-bottom corner and the positive sense of y axis is upwards. Images
      are stored and read from left-top corner, with positive y axis downwards;
      for that reason since now reading and storing data is done reverting 
      y axis (the boxes [iy] are changed to [dimy-1-iy]). This change would 
      cause no real effect in the programs using these routines neither 
      backwards compatibility should be compromised. 

      I adopted this new convention to avoid difficulties when retrieving 2D
      fields of physically formatted files.
*/


#define GRAFICOS_C


#ifndef STRUCT_DEF
#include <struct_def.c>
#endif

#ifndef TENSOR_C
#include <tensor.c>
#endif

/*	Dimension parameters (INRIA and vH)	*/

#define XMAXVH 1536
#define YMAXVH 1024
#define XMAXINRIA 1350
#define YMAXINRIA 460


/*                GIF   Parameters                              */

#define DELAY 50 // centesimas de segundo de retardo entre frames
#define WR 0.3
#define WG 0.59 // pesos de las componentes cromaticas en las calibraciones
#define WB 0.11




/*	        Function prototypes		*/


int ext_grafico( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		 char *nombre, Read2D *p_lee);
int ext_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	     char *nombre, Read2D *p_lee);
int ext_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	     char *nombre, Read2D *p_lee);
int ext_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		     char *nombre, Read2D *p_lee);
		 
int check_grafico( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		    char *nombre, Read2D *p_lee);
int check_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre, Read2D *p_lee);
int check_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre, Read2D *p_lee);
int check_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		       char *nombre, Read2D *p_lee);



int gif( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
	 char *nombre, double ***datos);
int ppm( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
	 char *nombre, double ***datos);
int unformatted( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
		 char *nombre, double ***datos);

int lee_en_ascii( FILE* canal);
int lee_double_int( int xmax, int ymax, int ix0, int iy0, int dimx, int dimy,
		    int littleindian, char* nombre, double **data);


int lee_mask_block( int dimx, int dimy, double block, char *nombre,
		    char **mask);
int lee( int dimx, int dimy, char* nombre_out, double **data);



int lee_gif( int dimx, int dimy, int nimag, char *nombre, char **Red,
	     char **Green, char **Blue);
int interpreta_extension_gif( FILE *canal, int verbose);
int lee_graphic_ce_gif( FILE *canal, int verbose);
int lee_comment_gif( FILE *canal, int verbose);
int lee_aplicacion_gif(FILE *canal, int verbose);
int interpreta_stream_gif( FILE *canal, unsigned char **encoded);
void graba_stream_gif( int l_code, unsigned char *encoded, FILE *canal);
int annade_decoded_gif(int w0, int letra, int *tab0, int *tab1, int id0, 
	int *decoded);
int primer_car_gif( int w0, int letra, int *tab0);


int graba_gif( int dimx, int dimy, char *nombre, char **Red, char **Green,
	       char **Blue);
int graba_gif_animado( int dimx, int dimy, int modo, int colort, 
		       char *nombre, char **Red, char **Green, char **Blue);
int adimensiona_gif( int size);
int crea_tabla_color_gif( int dimx, int dimy, char **Red, char **Green, 
			  char **Blue, int *dimI, int *dimC, int *dimY,
			  int **tabla_ICY);
void RGB_ICY( int dimI, int dimC, int dimY, int mode, char *Red, char *Green,
	      char *Blue, int *iI, int *iC, int *iY);
void codifica_imagen_con_tabla_ICY( int dimx, int dimy, int dimI, int dimC,
				    int dimY, char **Red, char **Green, 
				    char **Blue, int color_size, 
				    int *tabla_ICY, int *decoded);


void graba_foto(int dimx, int dimy, char* nombre, double **datos);
void graba_foto_block(int dimx, int dimy, double block, char* nombre,
	double **datos);
void graba_foto_block_limites(int dimx, int dimy, double block, char* nombre,
	double min_c, double max_c, double **datos);
void graba_foto_vec_block(int dimx, int dimy, double block, char* nombre,
	double **vx, double **vy);
void graba_foto_4(int dimx, int dimy, char* nombre, double **datos);
void graba_mask_block(int dimx, int dimy, double block, char* nombre,
	char **mask);
void graba_RGB_block(int dimx, int dimy, double block, char *nombre,
	char **Red, char **Green, char **Blue);
void graba_foto_color_block( int dimx, int dimy, int n_cr, double block,
	char *nombre, double ***datos);
void graba_video_block( int dimx, int dimy, int n_cr, double block, 
	int modo, char *nombre, double ***datos);
void graba_video_RGB_block( int dimx, int dimy, double block, 
	int modo, char *nombre, char **Red, char **Green, char **Blue);



void prepara_foto_block( int dimx, int dimy, double block, double **datos,
	char **foto);
void prepara_foto_fijo_block( int dimx, int dimy, double block, int levels,
	double **datos, char **foto);
void prepara_foto_block_limites( int dimx, int dimy, double min_c, 
	double max_c, double block, double **datos, char **foto);
void prepara_foto_log( int dimx, int dimy, double **datos, char **foto);
void prepara_foto_4( int dimx, int dimy, double **datos, char **foto);
void prepara_foto_log_4( int dimx, int dimy, double **datos, char **foto);
void prepara_char_block(int dimx, int dimy, double block, char **grayin,
	char **grayout);
int graba_pgm( int dimx, int dimy, int bin, char* nombre_out, char **cont);
int graba_ppm( int dimx, int dimy, int bin, char *nombre_out, char **Red, 
	char **Green, char **Blue);
int graba( int dimx, int dimy, char* nombre_out, double **data);


/*     Function declarations   */

int ext_grafico( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		 char *nombre, Read2D *p_lee)
{
    int error=-1;

    if(error) error=ext_gif(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
    if(error) error=ext_ppm(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
    if(error) error=ext_unformatted(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
  
    return(error);

}

int ext_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	     char *nombre, Read2D *p_lee)
{
    char ext[90];
    int error=0;
    int lext;

    lext=extrae_extension(nombre,'.',ext);


    if(lext!=3) error=-1;
    else if(
	((ext[0]!='g')&&(ext[0]!='G'))||
	((ext[1]!='i')&&(ext[1]!='I'))||
	((ext[2]!='f')&&(ext[2]!='F'))
	) error=-1;

    if(!error) error=check_gif(dimx,dimy,dimv,dimz,bd,nombre,p_lee);

    return(error);
}

int ext_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	     char *nombre, Read2D *p_lee)
{
    char ext[90];
    int error=0;
    int lext;

    lext=extrae_extension(nombre,'.',ext);

    if(lext!=3) error=-1;
    else if(
	((ext[0]!='p')&&(ext[0]!='P'))||
	((ext[1]!='p')&&(ext[1]!='P')&&(ext[1]!='g')&&(ext[1]!='G'))||
	((ext[2]!='m')&&(ext[2]!='M'))
	) error=-1;
 
    if(!error) error=check_ppm(dimx,dimy,dimv,dimz,bd,nombre,p_lee);

    return(error);
}

int ext_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		     char *nombre, Read2D *p_lee)
{
    char ext[90];
    int error=0;
    int lext;

    lext=extrae_extension(nombre,'.',ext);

    if(lext!=3) error=-1;
    else if(
	((ext[0]=='i')||(ext[0]=='I'))&&
	((ext[1]=='m')||(ext[1]=='M'))&&
	((ext[2]=='c')||(ext[2]=='C'))
	) error=0;
    else if(
	((ext[0]!='r')&&(ext[0]!='R'))||
	((ext[1]!='a')&&(ext[1]!='A'))||
	((ext[2]!='w')&&(ext[2]!='W'))
	) error=-1;

    if(!error) error=check_unformatted(dimx,dimy,dimv,dimz,bd,nombre,p_lee);

    return(error);
}


int check_grafico( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		    char *nombre, Read2D *p_lee)
{
    int error=-1;

    if(error) error=check_gif(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
    if(error) error=check_ppm(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
    if(error) error=check_unformatted(dimx,dimy,dimv,dimz,bd,nombre,p_lee);
    
    return(error);
}


int check_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre,  Read2D *p_lee)
{
    FILE *canal;
    char R,G,B;
    unsigned char buffer[255];
    unsigned char quest;
    unsigned char b0;

    int error;
    int gdimv;
    int xmax,ymax;
    int x0,y0;
    int ic,ie,id,id0,it,in;
    int lt,last_cd;
    int ix,iy;
    unsigned int color_res,color_size;
    int bits,bit0,w0,w1,l_code,letra;
    int verbose,global_ct,local_ct,fin,apren;


    verbose=0; // For debugging purposes
	
    canal=fopen(nombre,"rb");
    if(canal==NULL) return(-1);
 
    fread(buffer,sizeof(char),6,canal);
    if(feof(canal)) return(-1);
    buffer[6]='\0';
    if(verbose) printf("Version: %s\n",buffer);
    fread(buffer,sizeof(char),7,canal);
    if(feof(canal)) return(-1);
   

    *dimx=256*(int)buffer[1]+buffer[0];
    *dimy=256*(int)buffer[3]+buffer[2];
    gdimv=*dimv=1; // unless we see color we asume that images are monochrome

    if(verbose) printf("Dimensiones: %d x %d\n",*dimx,*dimy);
    if((buffer[4]&0x80)&&verbose) 
	printf("Color de fondo: %d\n",(int)buffer[5]);
    color_res=1+(int)(buffer[4]&0x70)/16;
    if(verbose) printf("Profundidad de color: %d bits\n",color_res);
    if(buffer[4]&0x80)
    {
	global_ct=1;
	color_size=1+(int)(buffer[4]&0x07);
	color_size=(int)pow(2.,(double)color_size);
	if((buffer[4]&0x2)&&verbose)
	    printf("Tabla de color ordenada\n");
	if(verbose)
	{
	    printf("Tamanno de la tabla de color: %d\n",
		   color_size);
	    printf("Razon de aspecto del pixel: %f\n",
		   ((double)buffer[6]+15.)/64);
	}
    }
    else
    {
	global_ct=0;
	if(verbose) printf("Imagen Gif sin tabla de color global\n");
    }
    if(global_ct==1)
    {
	for(ic=0;ic<color_size;ic++)
	{
	    R=getc(canal);
	    G=getc(canal);
	    B=getc(canal);
	    if((R!=G)||(R!=B)||(G!=B)) gdimv=3;
	}
    }

    for(quest=getc(canal),in=0;quest!=0x3B;quest=getc(canal))
    {
	if(quest==0x21) interpreta_extension_gif(canal,verbose);
	else if(quest!=0x2c)
	{
	    if(verbose) printf("Identificador desconocido!! Codigo: %d\n",
			       (int)quest);
	    return(-1);
	}
	else
	{
	    fread(buffer,sizeof(char),9,canal);
	    if(feof(canal)) return(-1);
	    x0=256*(int)buffer[1]+(int)buffer[0];
	    y0=256*(int)buffer[3]+(int)buffer[2];
	    xmax=256*(int)buffer[5]+(int)buffer[4];
	    ymax=256*(int)buffer[7]+(int)buffer[6];
	    local_ct=(int)(buffer[8]&0x80);
	    if(local_ct)
	    {
		color_size=1+(int)(buffer[8]&0x07);
		color_size=(int)pow(2.,(double)color_size);
	    }
	    else *dimv=Max(*dimv,gdimv);
	    if(verbose)
	    {
		printf("Imagen numero: %d\n",in);
		printf("Posicion de la ventana: (%d,%d)\n",
		       x0,y0);
		printf("Dimensiones: %d x %d\n",xmax,ymax);
		if(local_ct)
		{
		    printf("Tamanno de la tabla de color local: %d\n",
			   color_size);
		    if(buffer[8]&0x20) printf("Tabla ordenada\n");
		}
		if(buffer[8]&0x40) printf("Imagen entrelazada\n");
	    }
	    if(local_ct)
	    {
		for(ic=0;ic<color_size;ic++)
		{
		    R=getc(canal);
		    G=getc(canal);
		    B=getc(canal);
		    if((R!=G)||(R!=B)||(G!=B)) *dimv=3;
		}
	    }
	    
	    if((local_ct==0)&&(global_ct==0))
	    {
		if(verbose)
		{
		    printf("Imagen sin tabla de colores\n");
		    printf("Decodificacion imposible\n");
		}
		return(-1);
	    }
	    
/*		lectura y decodificacion de la imagen		*/

	    bit0=(int)getc(canal);
	    error=lee_comment_gif(canal,0);
	    in++;
	}
    }
    fclose(canal);

    if(!error)
    {
	if(verbose) printf("Numero de imagenes: %d\n",in);
	
	*dimz=in;
	if(*dimv==3) *bd=3;
	else *bd=1;

	*p_lee=&gif;
    }

    return(error);
}

int check_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre, Read2D *p_lee)
{
    FILE* canal;
    char iden[3],buffer[90];

    double value;

    int error=0;
    int levels;
    int ix,iy;

    *dimx=*dimy=*dimz=*dimv=*bd=0;

    canal=fopen(nombre,"rb");
    if(canal==NULL) return(-1);
    fread(iden,sizeof(char),3,canal);
    if(feof(canal)) return(-1);

    if((iden[0]!='P')||
       ((iden[1]!='2')&&(iden[1]!='3')&&(iden[1]!='5')&&(iden[1]!='6'))
	) return(-1);
    
    *dimx=lee_en_ascii(canal);  
    *dimy=lee_en_ascii(canal);
    if((iden[1]=='3')||(iden[1]=='6')) *dimv=3;
    else *dimv=1;
    *dimz=1;
    levels=1+lee_en_ascii(canal);
    *bd=(int)(log((double)levels)/log(256.));

    error=0;
    switch(iden[1])
    {

	case '2': // ASCII graylevel
	    for(iy=0;(iy<*dimy)&&(!error);iy++)
	    {
	    for(ix=0;(ix<*dimx)&&(!error);ix++)
	    {
		lee_en_ascii(canal);
		error=feof(canal);
	    }
	    }
	    break;
	case '3': // ASCII color
	    for(iy=0;(iy<*dimy)&&(!error);iy++)
	    {
	    for(ix=0;(ix<*dimx)&&(!error);ix++)
	    {
		lee_en_ascii(canal);
		lee_en_ascii(canal);
		lee_en_ascii(canal);
		error=feof(canal);
	    }
	    }
	    break;
	case '5': // Binary graylevel
	    for(iy=0;(iy<*dimy)&&(!error);iy++)
	    {
	    for(ix=0;(ix<*dimx)&&(!error);ix++)
	    {
		error=(fscanf(canal,"%c",&buffer)==EOF)?-1:0;
	    }
	    }
	    break;
	case '6': // Binary color
 	    for(iy=0;(iy<*dimy)&&(!error);iy++)
	    {
	    for(ix=0;(ix<*dimx)&&(!error);ix++)
	    {
		error=(fscanf(canal,"%c",&buffer)==EOF)?-1:0;
		error=(fscanf(canal,"%c",&buffer)==EOF)?-1:0;
		error=(fscanf(canal,"%c",&buffer)==EOF)?-1:0;
	    }
	    }
	    break;
    }
 
    if((ix==*dimx)&&(iy==*dimy))
    {	
 	if(fscanf(canal,"%c",buffer)==EOF) error=0; // File is over, ok!
	else error=-1;
    }
    else error=-1;
    fclose(canal);

    if(!error) *p_lee=&ppm;

    return(error);
}

int check_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		       char *nombre, Read2D *p_lee)
{
    FILE *canal;
    long int length;
    int error=0;

    canal=fopen(nombre,"rb");
    if(canal==NULL) return(-1);
    fseek(canal,0,SEEK_SET);
    length=-ftell(canal);
    fseek(canal,0,SEEK_END);
    length+=ftell(canal);
    fclose(canal);
    if(length==3145728) 
    {
	*dimx=XMAXVH;
	*dimy=YMAXVH;
	*dimv=1;
	*dimz=1;
	*bd=2;
    }
    else if(length==1242000)
    {
	*dimx=XMAXINRIA;
	*dimy=YMAXINRIA;
	*dimv=1;
	*dimz=1;
	*bd=1;
    }
    else error=-1;
    if(!error) *p_lee=&unformatted;


    return(error);
}


int gif( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
	 char *nombre, double ***datos)
{
    FILE *canal;
    double value;
    unsigned char buffer[255];
    unsigned char quest;
    unsigned char *Rg,*Gg,*Bg;
    unsigned char *Rl,*Gl,*Bl;
    unsigned char *encoded;
    unsigned char b0;
    int *decoded,*tab0,*tab1;
    int xmax,ymax;
    int x0,y0;
    int ibe,base;
    int ic,ie,id,id0,it,in;
    int lt,last_cd;
    int ix,iy;
    unsigned int color_res,color_size;
    int bits,bit0,w0,w1,l_code,letra;
    int verbose,global_ct,local_ct,fin,apren;


    verbose=0; // used for debugging purposes

    canal=fopen(nombre,"rb");
    fread(buffer,sizeof(char),6,canal);
    buffer[6]='\0';
    if(verbose) printf("Version: %s\n",buffer);
    fread(buffer,sizeof(char),7,canal);
    xmax=256*(int)buffer[1]+buffer[0];
    ymax=256*(int)buffer[3]+buffer[2];
    if(verbose) printf("Dimensiones: %d x %d\n",xmax,ymax);
    if((xmax>dimx)||(ymax>dimy))
    {
	printf("Imagen demasiado grande!!!\n");
	return(-1);
    }
    if((buffer[4]&0x80)&&verbose) printf("Color de fondo: %d\n",
					 (int)buffer[5]);
    color_res=1+(int)(buffer[4]&0x70)/16;
    if(verbose) printf("Profundidad de color: %d bits\n",color_res);
    if(buffer[4]&0x80)
    {
	global_ct=1;
	color_size=1+(int)(buffer[4]&0x07);
	color_size=(int)pow(2.,(double)color_size);
	if((buffer[4]&0x2)&&verbose)
	    printf("Tabla de color ordenada\n");
	if(verbose)
	{
	    printf("Tamanno de la tabla de color: %d\n",
		   color_size);
	    printf("Razon de aspecto del pixel: %f\n",
		   ((double)buffer[6]+15.)/64);
	}
    }
    else
    {
	global_ct=0;
	if(verbose) printf("Imagen Gif sin tabla de color global\n");
    }
    if(global_ct==1)
    {
	Rg=(char *) calloc(color_size,sizeof(unsigned char));
	Gg=(char *) calloc(color_size,sizeof(unsigned char));
	Bg=(char *) calloc(color_size,sizeof(unsigned char));

	for(ic=0;ic<color_size;ic++)
	{
	    Rg[ic]=getc(canal);
	    Gg[ic]=getc(canal);
	    Bg[ic]=getc(canal);
	}
    }

    for(quest=getc(canal),in=0;quest!=0x3B;quest=getc(canal))
    {
	if(quest==0x21) interpreta_extension_gif(canal,verbose);
	else if(quest!=0x2c)
	{
	    printf("Identificador desconocido!! Codigo: %d\n",
		   (int)quest);
	    return(-1);
	}
	else
	{
	    fread(buffer,sizeof(char),9,canal);
	    x0=256*(int)buffer[1]+(int)buffer[0];
	    y0=256*(int)buffer[3]+(int)buffer[2];
	    xmax=256*(int)buffer[5]+(int)buffer[4];
	    ymax=256*(int)buffer[7]+(int)buffer[6];
	    local_ct=(int)(buffer[8]&0x80);
	    if(local_ct)
	    {
		color_size=1+(int)(buffer[8]&0x07);
		color_size=(int)pow(2.,(double)color_size);
	    }
	    if(verbose)
	    {
		printf("Imagen numero: %d\n",in);
		printf("Posicion de la ventana: (%d,%d)\n",x0,y0);
		printf("Dimensiones: %d x %d\n",xmax,ymax);
		if(local_ct)
		{
		    printf("Tamanno de la tabla de color local: %d\n",
			   color_size);
		    if(buffer[8]&0x20) printf("Tabla ordenada\n");
		}
		if(buffer[8]&0x40) printf("Imagen entrelazada\n");
	    }
	    if(local_ct)
	    {
		Rl=(char *) calloc(color_size,sizeof(unsigned char));
		Gl=(char *) calloc(color_size,sizeof(unsigned char));
		Bl=(char *) calloc(color_size,sizeof(unsigned char));

		for(ic=0;ic<color_size;ic++)
		{
		    Rl[ic]=getc(canal);
		    Gl[ic]=getc(canal);
		    Bl[ic]=getc(canal);
		}
	    }

	    if((local_ct==0)&&(global_ct==0)&&(in==iz))
	    {
		printf("Imagen sin tabla de colores\n");
		printf("Decodificacion imposible\n");
		return(-1);
	    }

/*		lectura y decodificacion de la imagen		*/

	    bit0=(int)getc(canal);
	    if(in==iz)
	    {
		if(verbose)
		{
		    printf("\nEmpieza la decodificacion...\n");
		    printf("Tamanno inicial de palabra: %d bits\n",
			   bit0);
		}
		l_code=interpreta_stream_gif(canal,&encoded);
		if(verbose) printf("Tamanno del codigo: %d bytes\n",
				   l_code);
		decoded=(int *) calloc(xmax*ymax,sizeof(int));
		tab0=(int *) calloc(4096,sizeof(int));
		tab1=(int *) calloc(4096,sizeof(int));

	
		w0=(int) pow(2.,(double)bit0);
		w1=2*w0;
		bits=bit0+1;


		lt=0;
		last_cd=w0;
		id0=0;
		for(id=0,ie=0,b0=0x01,fin=1;fin==1;)
		{
		    if(id>=xmax*ymax) fin=3;
/*			Interpretando el caracter vigente		*/


		    base=1;
		    letra=0;
		    for(ibe=0;(ibe<bits)&&(fin!=2);ibe++)
		    {
			if(encoded[ie]&b0) letra+=base;
			base=2*base;
			b0=2*b0;
			if(b0==0x00)
			{
			    b0=0x01;
			    ie++;
			    if(ie>=l_code) fin=2;
			}
		    }

/*			Produciendo la salida asociada			*/

		    if(letra==w0)
		    {
			lt=0;
			bits=bit0+1;
			w1=2*w0;
		    }
		    else if(letra==w0+1) fin=0;
		    else if(fin==1)
		    {
			if(letra>lt+w0+2)
			{
			    fin=4;
			}
			else if(letra==lt+w0+2)
			{
/*		Aprendizaje de un codigo no visto	*/

			    apren=1;
			    tab0[lt]=last_cd;
			    tab1[lt]=primer_car_gif(w0,last_cd,tab0);
			    lt++;
			}
/*		Reproduccion del codigo			*/

			id=annade_decoded_gif(w0,letra,tab0,tab1,id,decoded);

		    }
/*		Aprendizaje de codigos vistos		*/


		    if((letra==w0)||(letra==w0+1)||(fin!=1)) apren=1;
		    if((apren!=1)&&(last_cd!=w0))
		    {

			tab0[lt]=last_cd;
			tab1[lt]=decoded[id0];

			for(ic=0;(ic<lt)&&(id0<id);ic++)
			{

			    if((tab0[ic]==tab0[lt])&&(tab1[ic]==tab1[lt]))
			    {
				tab0[lt]=ic+w0+2;
				id0++;
				if(id0<id) tab1[lt]=decoded[id0];
			    }
			}
			if(id0<id) lt++;
		    }
		    

/*		Actualizacion de variables		*/

		    if((lt+w0+2==w1)&&(letra!=w0))
		    {
			w1=2*w1;
			bits++;
			if(bits>12) bits=12;
		    }

		    last_cd=letra;
		    apren=0;
		    id0=id;
		    
		}
		if((fin==0)&&verbose) printf("Lectura exitosa!\n");
		if(fin==2) 
		    printf("Error: sennal fin de imagen no encontrada\n");
		if(fin==3) printf("Error: codigo desborda imagen\n");
		if(fin==4) printf("Error de lectura\n");
				

/*		Fin de la decodificacion		*/

		free(encoded);
		free(tab0);
		free(tab1);


		for(id=0,ix=0,iy=0;id<xmax*ymax;id++,ix++)
		{
		    if(ix>=xmax)
		    {
			ix-=xmax;
			iy++;
		    }
		    if(dimv==3)
		    {
			if(local_ct)
			{
			    datos[0][dimy-1-iy][ix]=(double)Rl[decoded[id]];
			    datos[1][dimy-1-iy][ix]=(double)Gl[decoded[id]];
			    datos[2][dimy-1-iy][ix]=(double)Bl[decoded[id]];
			}
			else
			{
			    datos[0][dimy-1-iy][ix]=(double)Rg[decoded[id]];
			    datos[1][dimy-1-iy][ix]=(double)Gg[decoded[id]];
			    datos[2][dimy-1-iy][ix]=(double)Bg[decoded[id]];
			}
		    }
		    else
		    {
			if(local_ct) 
			    datos[0][dimy-1-iy][ix]=
				WR*(double)Rl[decoded[id]]+
				WG*(double)Gl[decoded[id]]+
				WB*(double)Bl[decoded[id]];
			else
			    datos[0][dimy-1-iy][ix]=
				WR*(double)Rg[decoded[id]]+
				WG*(double)Gg[decoded[id]]+
				WB*(double)Bg[decoded[id]];
		    }
		}
		free(decoded);
	    }
	    else lee_comment_gif(canal,0);
	    in++;
	    
	    if(local_ct)
	    {
		free(Rl);
		free(Gl);
		free(Bl);
	    }
	    
	}
    }
    fclose(canal);
    if(global_ct)
    {
	free(Rg);
	free(Gg);
	free(Bg);
    }
    
    if(verbose) printf("Numero de imagenes: %d\n",in);
    
    if(in<iz)
    {
	printf("Imagen %d de %d no encontrada\n",iz,in);
	return(-1);
    }
    else return(color_size);

}

int ppm( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
	 char *nombre, double ***datos)
{
    FILE* canal;
    char iden[3];

    double value;

    int xmax,ymax,levels;
    int ix,iy;

    canal=fopen(nombre,"rb");
    fread(iden,sizeof(char),3,canal);

    xmax=lee_en_ascii(canal);  
    ymax=lee_en_ascii(canal);
    levels=1+lee_en_ascii(canal);

    switch((int)iden[1])
    {

	case '2': // ASCII graylevel
	    for(iy=0;iy<dimy;iy++)
	    {
	    for(ix=0;ix<dimx;ix++)
	    {
		datos[0][dimy-1-iy][ix]=(double) lee_en_ascii(canal);
	    }
	    }
	    break;
	case '3': // ASCII color
	    for(iy=0;iy<ymax;iy++)
	    {
            for(ix=0;ix<xmax;ix++)
	    {
		datos[0][dimy-1-iy][ix]=(double) lee_en_ascii(canal);
		datos[1][dimy-1-iy][ix]=(double) lee_en_ascii(canal);
		datos[2][dimy-1-iy][ix]=(double) lee_en_ascii(canal);
	    }
	    }
	    break;
	case '5': // Binary graylevel
	    for(iy=0;iy<dimy;iy++)
	    {
	    for(ix=0;ix<dimx;ix++)
	    {
		value=(double)fgetc(canal);
		if(value<0) value+=256.;
		datos[0][dimy-1-iy][ix]=value;
	    }
	    }
	    break;
	case '6': // Binary color
 	    for(iy=0;iy<dimy;iy++)
	    {
            for(ix=0;ix<dimx;ix++)
	    {
		value=(double)fgetc(canal);
		if(value<0) value+=256.;
		datos[0][dimy-1-iy][ix]=value;
		value=(double)fgetc(canal);
		if(value<0) value+=256.;
		datos[1][dimy-1-iy][ix]=value;
		value=(double)fgetc(canal);
		if(value<0) value+=256.;
		datos[2][dimy-1-iy][ix]=value;
		
	    }
	    }
	    break;
    }
	
    fclose(canal);
    return(levels);

}


int unformatted( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
		 char *nombre, double ***datos)
{
    int levels;

    if((dimx==XMAXVH)&&(dimy==YMAXVH))
	levels=lee_double_int(XMAXVH,YMAXVH,0,0,dimx,dimy,0,nombre,datos[0]);
    else levels=lee_double_int(XMAXINRIA,YMAXINRIA,0,0,dimx,dimy,1,
			       nombre,datos[0]);
    
}

int lee_en_ascii(FILE* canal)
{
	int comienzo,fin,comentario;
	int leo;
	int salida;

	comienzo=0;
	comentario=0;
	fin=0;
	salida=0;
	while(!fin)
	{
		leo=getc(canal);
		if(leo=='#') comentario=1;
		if(leo=='\n') comentario=0;
		if(!comentario)
		{
			if((leo>='0')&&(leo<='9')) 
			{
				salida=salida*10+(leo-'0');
				comienzo=1;
			}
			else if(comienzo) fin=1;
		}
	}

	return(salida);
}


int lee_double_int( int xmax, int ymax, int ix0, int iy0, int dimx, int dimy,
		    int littleindian, char* nombre, double **data)
{
    FILE* canal;
    double value;
    int dat1,dat2,ix,iy;

    canal=fopen(nombre,"rb");
    for(iy=0;iy<ymax;iy++)
    {
    for(ix=0;ix<xmax;ix++)
    {
	dat1=(int) getc(canal);
	if(dat1<0) dat1=dat1+256;
	dat2=(int) getc(canal);
	if(dat2<0) dat2=dat2+256;

	if(littleindian) value=(double)(dat1+256*dat2);
	else value=(double)(dat2+256*dat1);
	
	if((ix>=ix0)&&(ix<ix0+dimx)&&(iy>=iy0)&&(iy<iy0+dimy))
	    data[dimy-1-(iy-iy0)][ix-ix0]=value;
    }
    }
    fclose(canal);
    
    return(32768);
}


int lee( int dimx, int dimy, char *nombre_out, double **data)
{
	FILE* canal;
	int ix,iy;

	canal=fopen(nombre_out,"rb");
	for(iy=0;iy<dimy;iy++) fread(data[iy],sizeof(double),dimx,canal);
	fclose(canal);
	return(0);
}



int lee_gif( int dimx, int dimy, int nimag, char *nombre, char **Red, 
	     char **Green, char **Blue)
{
    FILE *canal;
    unsigned char buffer[255];
    unsigned char quest;
    char *Rg,*Gg,*Bg;
    char *Rl,*Gl,*Bl;
    unsigned char *encoded;
    unsigned char b0;
    int *decoded,*tab0,*tab1;
    int xmax,ymax;
    int x0,y0;
    int ibe,base;
    int ic,ie,id,id0,it,in;
    int lt,last_cd;
    int ix,iy;
    unsigned int color_res,color_size;
    int bits,bit0,w0,w1,l_code,letra;
    int verbose,global_ct,local_ct,fin,apren;


    verbose=0; // used for debugging purposes

    canal=fopen(nombre,"rb");
    fread(buffer,sizeof(char),6,canal);
    buffer[6]='\0';
    if(verbose) printf("Version: %s\n",buffer);
    fread(buffer,sizeof(char),7,canal);
    xmax=256*(int)buffer[1]+buffer[0];
    ymax=256*(int)buffer[3]+buffer[2];
    if(verbose) printf("Dimensiones: %d x %d\n",xmax,ymax);
    if((xmax>dimx)||(ymax>dimy))
    {
	printf("Imagen demasiado grande!!!\n");
	return(-1);
    }
    if((buffer[4]&0x80)&&verbose) printf("Color de fondo: %d\n",
					 (int)buffer[5]);
    color_res=1+(int)(buffer[4]&0x70)/16;
    if(verbose) printf("Profundidad de color: %d bits\n",color_res);
    if(buffer[4]&0x80)
    {
	global_ct=1;
	color_size=1+(int)(buffer[4]&0x07);
	color_size=(int)pow(2.,(double)color_size);
	if((buffer[4]&0x2)&&verbose)
	    printf("Tabla de color ordenada\n");
	if(verbose)
	{
	    printf("Tamanno de la tabla de color: %d\n",
		   color_size);
	    printf("Razon de aspecto del pixel: %f\n",
		   ((double)buffer[6]+15.)/64);
	}
    }
    else
    {
	global_ct=0;
	if(verbose) printf("Imagen Gif sin tabla de color global\n");
    }
    if(global_ct==1)
    {
	Rg=(char *) calloc(color_size,sizeof(char));
	Gg=(char *) calloc(color_size,sizeof(char));
	Bg=(char *) calloc(color_size,sizeof(char));

	for(ic=0;ic<color_size;ic++)
	{
	    Rg[ic]=getc(canal);
	    Gg[ic]=getc(canal);
	    Bg[ic]=getc(canal);
	}
    }

    for(quest=getc(canal),in=0;quest!=0x3B;quest=getc(canal))
    {
	if(quest==0x21) interpreta_extension_gif(canal,verbose);
	else if(quest!=0x2c)
	{
	    printf("Identificador desconocido!! Codigo: %d\n",
		   (int)quest);
	    return(-1);
	}
	else
	{
	    fread(buffer,sizeof(char),9,canal);
	    x0=256*(int)buffer[1]+(int)buffer[0];
	    y0=256*(int)buffer[3]+(int)buffer[2];
	    xmax=256*(int)buffer[5]+(int)buffer[4];
	    ymax=256*(int)buffer[7]+(int)buffer[6];
	    local_ct=(int)(buffer[8]&0x80);
	    if(local_ct)
	    {
		color_size=1+(int)(buffer[8]&0x07);
		color_size=(int)pow(2.,(double)color_size);
	    }
	    if(verbose)
	    {
		printf("Imagen numero: %d\n",in);
		printf("Posicion de la ventana: (%d,%d)\n",x0,y0);
		printf("Dimensiones: %d x %d\n",xmax,ymax);
		if(local_ct)
		{
		    printf("Tamanno de la tabla de color local: %d\n",
			   color_size);
		    if(buffer[8]&0x20) printf("Tabla ordenada\n");
		}
		if(buffer[8]&0x40) printf("Imagen entrelazada\n");
	    }
	    if(local_ct)
	    {
		Rl=(char *) calloc(color_size,sizeof(char));
		Gl=(char *) calloc(color_size,sizeof(char));
		Bl=(char *) calloc(color_size,sizeof(char));

		for(ic=0;ic<color_size;ic++)
		{
		    Rl[ic]=getc(canal);
		    Gl[ic]=getc(canal);
		    Bl[ic]=getc(canal);
		}
	    }

	    if((local_ct==0)&&(global_ct==0)&&(in==nimag))
	    {
		printf("Imagen sin tabla de colores\n");
		printf("Decodificacion imposible\n");
		return(-1);
	    }

/*		lectura y decodificacion de la imagen		*/

	    bit0=(int)getc(canal);
	    if(in==nimag)
	    {
		if(verbose)
		{
		    printf("\nEmpieza la decodificacion...\n");
		    printf("Tamanno inicial de palabra: %d bits\n",
			   bit0);
		}
		l_code=interpreta_stream_gif(canal,&encoded);
		if(verbose) printf("Tamanno del codigo: %d bytes\n",
				   l_code);
		decoded=(int *) calloc(xmax*ymax,sizeof(int));
		tab0=(int *) calloc(4096,sizeof(int));
		tab1=(int *) calloc(4096,sizeof(int));

	
		w0=(int) pow(2.,(double)bit0);
		w1=2*w0;
		bits=bit0+1;


		lt=0;
		last_cd=w0;
		id0=0;
		for(id=0,ie=0,b0=0x01,fin=1;fin==1;)
		{
		    if(id>=xmax*ymax) fin=3;
/*			Interpretando el caracter vigente		*/


		    base=1;
		    letra=0;
		    for(ibe=0;(ibe<bits)&&(fin!=2);ibe++)
		    {
			if(encoded[ie]&b0) letra+=base;
			base=2*base;
			b0=2*b0;
			if(b0==0x00)
			{
			    b0=0x01;
			    ie++;
			    if(ie>=l_code) fin=2;
			}
		    }

/*			Produciendo la salida asociada			*/

		    if(letra==w0)
		    {
			lt=0;
			bits=bit0+1;
			w1=2*w0;
		    }
		    else if(letra==w0+1) fin=0;
		    else if(fin==1)
		    {
			if(letra>lt+w0+2)
			{
			    fin=4;
			}
			else if(letra==lt+w0+2)
			{
/*		Aprendizaje de un codigo no visto	*/

			    apren=1;
			    tab0[lt]=last_cd;
			    tab1[lt]=primer_car_gif(w0,last_cd,tab0);
			    lt++;
			}
/*		Reproduccion del codigo			*/

			id=annade_decoded_gif(w0,letra,tab0,tab1,id,decoded);

		    }
/*		Aprendizaje de codigos vistos		*/


		    if((letra==w0)||(letra==w0+1)||(fin!=1)) apren=1;
		    if((apren!=1)&&(last_cd!=w0))
		    {

			tab0[lt]=last_cd;
			tab1[lt]=decoded[id0];

			for(ic=0;(ic<lt)&&(id0<id);ic++)
			{

			    if((tab0[ic]==tab0[lt])&&(tab1[ic]==tab1[lt]))
			    {
				tab0[lt]=ic+w0+2;
				id0++;
				if(id0<id) tab1[lt]=decoded[id0];
			    }
			}
			if(id0<id) lt++;
		    }
		    

/*		Actualizacion de variables		*/

		    if((lt+w0+2==w1)&&(letra!=w0))
		    {
			w1=2*w1;
			bits++;
			if(bits>12) bits=12;
		    }

		    last_cd=letra;
		    apren=0;
		    id0=id;
		    
		}
		if((fin==0)&&verbose) printf("Lectura exitosa!\n");
		if(fin==2) 
		    printf("Error: sennal fin de imagen no encontrada\n");
		if(fin==3) printf("Error: codigo desborda imagen\n");
		if(fin==4) printf("Error de lectura\n");
				

/*		Fin de la decodificacion		*/

		free(encoded);
		free(tab0);
		free(tab1);


		for(id=0,ix=0,iy=0;id<xmax*ymax;id++,ix++)
		{
		    if(ix>=xmax)
		    {
			ix-=xmax;
			iy++;
		    }
		    if(local_ct)
		    {
			Red[dimy-1-iy][ix]=Rl[decoded[id]];
			Green[dimy-1-iy][ix]=Gl[decoded[id]];
			Blue[dimy-1-iy][ix]=Bl[decoded[id]];
		    }
		    else
		    {
			Red[dimy-1-iy][ix]=Rg[decoded[id]];
			Green[dimy-1-iy][ix]=Gg[decoded[id]];
			Blue[dimy-1-iy][ix]=Bg[decoded[id]];
		    }
		}
		free(decoded);
	    }
	    else lee_comment_gif(canal,0);
	    in++;
	    
	    if(local_ct)
	    {
		free(Rl);
		free(Gl);
		free(Bl);
	    }
	    
	}
    }
    fclose(canal);
    if(global_ct)
    {
	free(Rg);
	free(Gg);
	free(Bg);
    }
    
    if(verbose) printf("Numero de imagenes: %d\n",in);
    
    if(in<nimag)
    {
	printf("Imagen no encontrada\n");
	return(-1);
    }
    else return(color_size);
    
}


int interpreta_extension_gif( FILE *canal, int verbose)
{
    unsigned char decide;
    int error;

    decide=getc(canal);

    switch(decide)
    {
	case 0xF9:
	    error=lee_graphic_ce_gif(canal,verbose);
	    break;
	case 0xFE:
	    error=lee_comment_gif(canal,verbose);
	    break;
	case 0xFF:
	    error=lee_aplicacion_gif(canal,verbose);
	    break;
	default:
	    error=lee_comment_gif(canal,0);
    }
    return(error);
}

int lee_graphic_ce_gif( FILE *canal, int verbose)
{
    unsigned char buffer[6];
    unsigned int disposal,delay;

    fread(buffer,sizeof(char),6,canal);
    if(feof(canal)) return(-1);

    if(verbose)
    {
	printf("\nGraphical control extension parameters:\n");
	printf("=======================================\n");
	disposal=(int) (buffer[1]&0x70)/16;
	switch(disposal)
	{
	    case 0:
		printf("Disposicion: ninguna\n");
		break;
	    case 1:
		printf("Disposicion: ultimo\n");
		break;
	    case 2:
		printf("Disposicion: fondo\n");
		break;
	    case 3:
		printf("Disposicion: previo\n");
		break;
	    default:
		printf("Disposicion: indefinida\n");
	}
	if(buffer[1]&0x02) printf("Interaccion con el usuario\n");
	if(buffer[1]&0x01)
	    printf("Color transparente: %d\n",(int) buffer[4]);
	delay=10*(256*(int) buffer[3]+(int)buffer[2]);
	printf("Tiempo de retardo: %d ms\n",delay);
    }
    return(0);
}

int lee_comment_gif( FILE *canal, int verbose)
{
    unsigned char *buff;
    unsigned char quest;
    unsigned int size;
    int error;

    if(verbose) printf("Comentario:\n");
	
    error=0;
    for(quest=getc(canal);(quest!=0x00)&&(!error);quest=getc(canal))
    {
	size=(int)quest;
	buff=(unsigned char *)calloc(size,sizeof(unsigned char));
	fread(buff,sizeof(unsigned char),size,canal);
	if(feof(canal)) error=-1;
	if(verbose) printf("%s",buff);
	free(buff);
    }
    if(verbose) printf("\n");
    return(error);
}

int lee_aplicacion_gif(FILE *canal, int verbose)
{
    unsigned char *buff;
    unsigned char quest;
    int size,ic;
    int error;

    buff=(unsigned char *)calloc(12,sizeof(unsigned char));
    fread(buff,sizeof(unsigned char),1,canal);
    if(feof(canal)) return(-1);
    fread(buff,sizeof(unsigned char),8,canal);
    if(feof(canal)) return(-1);
    buff[8]='\0';
    fread(buff+9,sizeof(unsigned char),3,canal);
    if(feof(canal)) return(-1);
    if(verbose)
    {
	printf("Aplicacion: %s\n",buff);
	printf("Bytes de autentificacion: %d %d %d\n",
	       (int)buff[9],(int)buff[10],(int)buff[11]);
    }
    free(buff);

    error=0;
    for(quest=getc(canal),ic=0;(quest!=0x00)&&(!error);quest=getc(canal))
    {
	size=(int)quest;
	ic+=size;
	buff=(unsigned char *)calloc(size,sizeof(unsigned char));
	fread(buff,sizeof(unsigned char),size,canal);
	if(feof(canal)) error=-1;
	free(buff);
    }
    if(verbose) printf("Tamanno de aplicacion: %d\n",ic);

    return(error);
}


int interpreta_stream_gif( FILE *canal, unsigned char **encoded)
{
    char debuff[255];
    unsigned char quest;
    long int pos;
    int l_code,ib;

    pos=ftell(canal);
    l_code=0;
    for(quest=getc(canal);quest!=0x00;quest=getc(canal))
    {
	l_code+=(int) quest;
	fseek(canal,(int)quest,SEEK_CUR);
    }
    encoded[0]=(unsigned char *)calloc(l_code,sizeof(unsigned char));

    fseek(canal,pos,SEEK_SET);
    for(quest=getc(canal),ib=0;quest!=0x00;quest=getc(canal))
    {
	fread(encoded[0]+ib,sizeof(unsigned char),(int)quest,canal);
	ib+=(int)quest;
    }

    return(l_code);
}

void graba_stream_gif( int l_code, unsigned char *encoded, FILE *canal)
{
    unsigned char bsize;
    int ib;

    bsize=0xFF;
    for(ib=0;ib<l_code-(int)bsize;ib+=(int)bsize)
    {
	fwrite(&bsize,sizeof(char),1,canal);
	fwrite(&(encoded[ib]),sizeof(char),(int)bsize,canal);
    }
    bsize=(char)(l_code-ib);
    fwrite(&bsize,sizeof(char),1,canal);
    fwrite(&(encoded[ib]),sizeof(unsigned char),(int)bsize,canal);
    bsize=0x00;
    fwrite(&bsize,sizeof(char),1,canal);

}


int annade_decoded_gif(int w0, int letra, int *tab0, int *tab1, int id0, 
	int *decoded)
{
    int letra0,letra1;
    int id;

    id=id0;

    if(letra<w0)
 	decoded[id++]=letra;
    else
    {
	letra0=tab0[letra-w0-2];
	letra1=tab1[letra-w0-2];
	id=annade_decoded_gif(w0,letra0,tab0,tab1,id,decoded);
	id=annade_decoded_gif(w0,letra1,tab0,tab1,id,decoded);
    }
    
    return(id);

}

int primer_car_gif( int w0, int letra, int *tab0)
{
    int out;

    if(letra<w0) out=letra;
    else out=primer_car_gif(w0,tab0[letra-w0-2],tab0);

    return(out);

}



int graba_gif( int dimx, int dimy, char *nombre, char **Red, char **Green,
	       char **Blue)
{
    FILE *canal;
    unsigned char buffer[255];
    unsigned char desize;
    char auxR,auxG,auxB;

    int *tabla_ICY;
    int color_size;
    int aux;

    unsigned char *encoded;
    unsigned char b0;
    int *decoded,*tab0,*tab1;
    int dimI,dimC,dimY;
    int x0,y0;
    int ibe;
    int ic,ie,id;
    int lt,last_cd;
    int ix,iy;
    int iI,iC,iY;
    int bits,bit0,w0,w1,letra;
    int verbose,local_ct;


    verbose=0;

/*	Creacion de la tabla de color asociada		*/

    color_size=crea_tabla_color_gif(dimx,dimy,Red,Green,Blue,
				    &dimI,&dimC,&dimY,&tabla_ICY);

    if(color_size>256) bit0=8;
    else if(color_size==1) bit0=1;
    else bit0=adimensiona_gif(color_size-1);
    desize=(char) (bit0-1);




    canal=fopen(nombre,"wb");
    sprintf(buffer,"GIF87a");
    fwrite(buffer,sizeof(char),6,canal);

    buffer[1]=(char)(dimx/256);
    buffer[0]=(char)(dimx-256*(int)buffer[1]);
    buffer[3]=(char)(dimy/256);
    buffer[2]=(char)(dimy-256*(int)buffer[3]);
    buffer[4]=0x70|desize; 
// lo cual significa: no hay tabla de color global, profundidad de color
// maxima, tabla no ordenada; se introduce, no obstante, el tamanno de
// la tabla de color para que las imagenes sean reconocibles al xv.

    buffer[5]=0x00; // color de fondo global, no dado
    buffer[6]=0x00; //pixel aspect ratio
	
    fwrite(buffer,sizeof(char),7,canal);


//	Comienza la grabacion de la imagen



    buffer[0]=0x2C; //identificador de imagen
    buffer[1]=0x00;
    buffer[2]=0x00; //Left position
    buffer[3]=0x00;
    buffer[4]=0x00; //Top position
    buffer[6]=(char)(dimx/256);
    buffer[5]=(char)(dimx-256*(int)buffer[6]);
    buffer[8]=(char)(dimy/256);
    buffer[7]=(char)(dimy-256*(int)buffer[8]);
    buffer[9]=0xA0|desize; // que significa: tabla local si,
// entrelazado no, tabla ordenada y desize codifica el tamanno de la tabla


    fwrite(buffer,sizeof(char),10,canal);

    for(ic=0;ic<color_size;ic++)
    {
	iI=tabla_ICY[ic]%dimI;
	iC=((tabla_ICY[ic]-iI)/dimI)%dimC;
	iY=(tabla_ICY[ic]-iI-dimI*iC)/(dimI*dimC);
	RGB_ICY(dimI,dimC,dimY,1,&auxR,&auxG,&auxB,&iI,&iC,&iY);

	putc(auxR,canal);
	putc(auxG,canal);
	putc(auxB,canal);
    }
    if((color_size>1)&&(color_size<256)) aux=dimensiona(color_size);
    else aux=2;
    if(color_size<aux)
    {
	for(ic=color_size;ic<aux;ic++)
	{
	    putc(0,canal);
	    putc(0,canal);
	    putc(0,canal);
	}
    }

	

/*		Codificacion de acuerdo a la tabla		*/

    decoded=(int *) calloc(dimx*dimy,sizeof(int));
    codifica_imagen_con_tabla_ICY(dimx,dimy,dimI,dimC,dimY,Red,Green,Blue,
				  color_size,tabla_ICY,decoded);

/*		Codificacion de la imagen		*/


    encoded=(unsigned char *) calloc((3*dimx*dimy)/2,sizeof(char));
    tab0=(int *) calloc(4096,sizeof(int));
    tab1=(int *) calloc(4096,sizeof(int));


/*		Registro del tamanno inicial de palabra	*/

    if(bit0<2) bit0=2;
    buffer[0]=(char) bit0;
    fwrite(buffer,sizeof(unsigned char),1,canal);


/*		Inicializacion			*/

    w0=(int) pow(2.,(double)bit0);
    encoded[0]=0;
    last_cd=w0;
    w1=2*w0;
    bits=bit0+1;
 
/*		Graba un reset como primer codigo	*/
    ie=0;
    b0=0x01;
    letra=w0;
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)
    {
	if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
	b0=2*b0;
	if(b0==0x00)
	{
	    b0=0x01;
	    ie++;
	    encoded[ie]=0;
	}
    }

/*		Bucle codificador		*/

    for(id=0;id<dimx*dimy;)
    {

/*		Examinando la cadena		*/

	if(bits==13)
	{
	    bits=12;
	    letra=w0;
	}
	else
	{
	    if(last_cd==w0)
	    {
		lt=0;
		bits=bit0+1;
		w1=2*w0;
	    }

	    letra=decoded[id++];
	    for(ic=0;(ic<lt)&&(id<dimx*dimy);ic++)
	    {
		if((tab0[ic]==letra)&&(tab1[ic]==decoded[id]))
		{
		    letra=w0+2+ic;
		    id++;
		}
	    }

//	Las cadenas aprendidas corresponden siempre a la ultima cadena 
//		registrada mas el caracter siguiente.

	    if(id<dimx*dimy)
	    {
		tab0[lt]=letra;
		tab1[lt++]=decoded[id];
	    }

	}

/*		Actualizacion de variables		*/


	last_cd=letra;


/*		Registrando el caracter vigente		*/

	for(ibe=0;ibe<bits;ibe++,letra=letra/2)
	{
	    if(Mod(letra,2)) encoded[ie]=encoded[ie]|b0;
	    b0=2*b0;
	    if(b0==0x00)
	    {
		b0=0x01;
		ie++;
		encoded[ie]=0;
	    }
	}

	if(lt+w0+1==w1)
	{
	    w1=2*w1;
	    bits++;
	}

    }
    if(bits==13) bits=12;
    letra=w0+1;
    for(ibe=0;ibe<bits;ibe++,letra=letra/2)
    {
	if(letra%2) encoded[ie]=encoded[ie]|b0;
	b0=2*b0;
	if(b0==0x00)
	{
	    b0=0x01;
	    ie++;
	    encoded[ie]=0;
	}
    }
    ie++; //ie representa ahora la longitud del codificado


/*		Fin de la codificacion		*/

    free(decoded);
    free(tab0);
    free(tab1);


/*	    Grabacion del codigo resultante	*/

    graba_stream_gif(ie,encoded,canal);
    free(encoded);


/*		Trailer y cierre		*/

    buffer[0]=0x3B;
    fwrite(buffer,sizeof(char),1,canal);
    fclose(canal);

    return(0);

}


int graba_gif_animado( int dimx, int dimy, int modo, int colort, 
	char *nombre, char **Red, char **Green, char **Blue)
{
    FILE *canal;
    char cabecera[11]="NETSCAPE2.0";
    unsigned char buffer[255];
    unsigned char desize;
    char auxR,auxG,auxB;

    int *tabla_ICY;
    int aux,color_size;

    unsigned char *encoded;
    unsigned char b0;
    unsigned char bsize;
    int *decoded,*tab0,*tab1;
    int dimI,dimC,dimY;
    int x0,y0;
    int ib;
    int ibe;
    int ic,ie,id;
    int lt,last_cd;
    int ix,iy;
    int iI,iC,iY;
    int bits,bit0,w0,w1,letra;
    int verbose;


    verbose=0;


    switch(modo)
    {
	case 0:

//	MODO 0: Grabacion de la cabecera


	    canal=fopen(nombre,"wb");
	    sprintf(buffer,"GIF87a");
	    fwrite(buffer,sizeof(char),6,canal);

	    buffer[1]=(char)(dimx/256);
	    buffer[0]=(char)(dimx-256*(int)buffer[1]);
	    buffer[3]=(char)(dimy/256);
	    buffer[2]=(char)(dimy-256*(int)buffer[3]);
	    
	    buffer[4]=0xFF; 
// lo cual significa: tabla de color global, profundidad de color
// maxima, tabla ordenada y de tamanno maximo


	    buffer[5]=0x00; // color de fondo global: el negro
	    buffer[6]=0x00; //pixel aspect ratio
	    fwrite(buffer,sizeof(char),7,canal);

		
	    for(ic=0;ic<256;ic++)
	    {
		putc((char)ic,canal);
		putc((char)ic,canal);
		putc((char)ic,canal);
	    }
		


//	Insercion de aplicacion asociada a la animacion

	    buffer[0]=0x21; //extension
	    buffer[1]=0xFF; //de aplicacion
	    buffer[2]=0x0B; //que mide 11 bytes
	    buffer[3]='\0';
	    strcat(buffer,cabecera); //identificada como NETSCAPE
	    buffer[14]=0x03;
	    buffer[15]=0x01;
	    buffer[16]=0x00;
	    buffer[17]=0x00;
	    buffer[18]=0x00;
	    fwrite(buffer,sizeof(unsigned char),19,canal);



	    fclose(canal);
	    break;
	case 1:

//	MOD0 1: Insercion de un frame

/*	Creacion de la tabla de color asociada		*/

	    if(colort)
	    {
		color_size=crea_tabla_color_gif(dimx,dimy,Red,Green,Blue,
						&dimI,&dimC,&dimY,
						&tabla_ICY);

		if(color_size>256) bit0=8;
		else if(color_size==1) bit0=1;
		else bit0=adimensiona_gif(color_size-1);
		desize=(char) (bit0-1);
	    }
	    else bit0=8;
	    
	    canal=fopen(nombre,"ab");
	    fseek(canal,0,SEEK_END);

//	Grabacion de la extension grafica de control

	    buffer[0]=0x21; //extension
	    buffer[1]=0xF9; //grafica de control
	    buffer[2]=0x04; //que mide 4 bytes
	    buffer[3]=0x00; // ni disposicion ni usuario ni transparencia
	    buffer[4]=Mod(DELAY,256);
	    buffer[5]=DELAY/256;  // tiempo de retardo entre frames
	    buffer[6]=0x00; // color transparente
	    buffer[7]=0x00; //terminador
	    fwrite(buffer,sizeof(unsigned char),8,canal);

//	Comienza la grabacion de la imagen


	    buffer[0]=0x2C; //identificador de imagen
	    buffer[1]=0x00;
	    buffer[2]=0x00; //Left position
	    buffer[3]=0x00;
	    buffer[4]=0x00; //Top position
	    buffer[6]=(char)(dimx/256);
	    buffer[5]=(char)(dimx-256*(int)buffer[6]);
	    buffer[8]=(char)(dimy/256);
	    buffer[7]=(char)(dimy-256*(int)buffer[8]);
		

		
	    if(colort) buffer[9]=0xAF; // que significa:
// tabla local si,  entrelazado no, tabla ordenada de tamanno maximo, forzado
	    else buffer[9]=0x0;


	    fwrite(buffer,sizeof(char),10,canal);

	    if(colort)
	    {
		for(ic=0;ic<color_size;ic++)
		{
		    iI=tabla_ICY[ic]%dimI;
		    iC=((tabla_ICY[ic]-iI)/dimI)%dimC;
		    iY=(tabla_ICY[ic]-iI-dimI*iC)/(dimI*dimC);
		    RGB_ICY(dimI,dimC,dimY,1,&auxR,&auxG,&auxB,&iI,&iC,&iY);

		    putc(auxR,canal);
		    putc(auxG,canal);
		    putc(auxB,canal);
 		}
		if((color_size>1)&&(color_size<256)) 
		    aux=dimensiona(color_size);
		else aux=2;

		if(color_size<256)
		{
		    for(ic=color_size;ic<256;ic++)
		    {
			putc(0,canal);
			putc(0,canal);
			putc(0,canal);
		    }
		}
	    }
	    color_size=Max(color_size,256);
	    bit0=8; // La tabla pasa a tener al menos 256 colores, 
	    // forzado por compatibilidad con xv.



/*		Codificacion de acuerdo a la tabla		*/

	    decoded=(int *) calloc(dimx*dimy,sizeof(int));
		

	    codifica_imagen_con_tabla_ICY(dimx,dimy,dimI,dimC,dimY,
					  Red,Green,Blue,
					  color_size,tabla_ICY,
					  decoded);

/*		Codificacion de la imagen		*/


	    encoded=(unsigned char *) calloc((3*dimx*dimy)/2,sizeof(char));
	    tab0=(int *) calloc(4096,sizeof(int));
	    tab1=(int *) calloc(4096,sizeof(int));
	

/*		Registro del tamanno inicial de palabra	*/

	    if(bit0<2) bit0=2;
	    buffer[0]=(char) bit0;
	    fwrite(buffer,sizeof(unsigned char),1,canal);


/*		Inicializacion			*/

	    w0=(int) pow(2.,(double)bit0);
	    encoded[0]=0;
	    last_cd=w0;
	    w1=2*w0;
	    bits=bit0+1;
 
/*		Graba un reset como primer codigo	*/

	    ie=0;
	    b0=0x01;
	    letra=w0;
	    for(ibe=0;ibe<bits;ibe++,letra=letra/2)
	    {
		if(letra%2) encoded[ie]=encoded[ie]|b0;
		b0=2*b0;
		if(b0==0x00)
		{
		    b0=0x01;
		    ie++;
		    encoded[ie]=0;
		}
	    }

/*		Bucle codificador		*/

	    for(id=0;id<dimx*dimy;)
	    {
		
/*		Examinando la cadena		*/

		if(bits==13)
		{
		    bits=12;
		    letra=w0;
		}
		else
		{
		    if(last_cd==w0)
		    {
			lt=0;
			bits=bit0+1;
			w1=2*w0;
		    }
		    
		    letra=decoded[id++];
		    for(ic=0;(ic<lt)&&(id<dimx*dimy);ic++)
		    {
			if((tab0[ic]==letra)&&(tab1[ic]==decoded[id]))
			{
			    letra=w0+2+ic;
			    id++;
			}
		    }
		    
//	Las cadenas aprendidas corresponden siempre a la ultima cadena 
//		registrada mas el caracter siguiente.

		    if(id<dimx*dimy)
		    {
			tab0[lt]=letra;
			tab1[lt++]=decoded[id];
		    }
		    
		}
		
/*		Actualizacion de variables		*/
		

		last_cd=letra;
		
		
/*		Registrando el caracter vigente		*/

		for(ibe=0;ibe<bits;ibe++,letra=letra/2)
		{
		    if(letra%2) encoded[ie]=encoded[ie]|b0;
		    b0=2*b0;
		    if(b0==0x00)
		    {
			b0=0x01;
			ie++;
			encoded[ie]=0;
		    }
		}

		if(lt+w0+1==w1)
		{
		    w1=2*w1;
		    bits++;
		}
		
	    }
	    if(bits==13) bits=12;
	    letra=w0+1;
	    for(ibe=0;ibe<bits;ibe++,letra=letra/2)
	    {
		if(letra%2) encoded[ie]=encoded[ie]|b0;
		b0=2*b0;
		if(b0==0x00)
		{
		    b0=0x01;
		    ie++;
		    encoded[ie]=0;
		}
	    }
	    ie++; //ie representa ahora la longitud del codificado
	    
	    
/*		Fin de la codificacion		*/

	    free(decoded);
	    free(tab0);
	    free(tab1);
	

/*	    Grabacion del codigo resultante	*/


	    graba_stream_gif(ie,encoded,canal);
	    free(encoded);
	    fclose(canal);
	    break;
	case 2:
	default:

//	Modo por defecto: inserta el trailer y libera la memoria
//	asociada a la tabla de color


/*		Trailer y cierre		*/
	
	    canal=fopen(nombre,"ab");
	    fseek(canal,0,SEEK_END);
	    buffer[0]=0x3B;
	    fwrite(buffer,sizeof(char),1,canal);
	    fclose(canal);
	    break;
    }

    return(0);

}

int adimensiona_gif( int size)
{
	int out,s0;


	for(out=0,s0=size;s0>0;s0=s0/2)	out++;

	return(out);
}

int crea_tabla_color_gif( int dimx, int dimy, char **Red, char **Green, 
			  char **Blue, int *dimI, int *dimC, int *dimY,
			  int **tabla_ICY)
{
    char auxR,auxG,auxB;
    double *freq;
    int *histo;
    int color_size;
    int ix,iy,iI,iC,iY,it;
    


/*   Initializing Intensity-Chrominance-Yellow space   */

    *dimI=256; // 2**8
    *dimC=64;  // 2**6
    *dimY=64;  // 2**6. With these choices, the ICY space is described with
               // 2**20 values, 16 times less than the maximum possible.
               // This represents 4 Mb in memory, reasonable nowadays.


    histo=(int *) calloc(*dimI**dimC**dimY,sizeof(int));


/*	Obtaining the color table already reduced in ICY space	*/

    color_size=0;
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	RGB_ICY(*dimI,*dimC,*dimY,0,&Red[iy][ix],&Green[iy][ix],&Blue[iy][ix],
		&iI,&iC,&iY);

	if(histo[iI+*dimI*(iC+*dimC*iY)]==0) color_size++;
	histo[iI+*dimI*(iC+*dimC*iY)]+=1;
    }
    }


    for(;color_size>256;) // performing chromatical reduction
    {
	color_size=0;
	*dimI/=2;
	*dimC/=2;
	*dimY/=2;
	for(iY=0;iY<*dimY;iY++)
	{
	for(iC=0;iC<*dimC;iC++)
	{
	for(iI=0;iI<*dimI;iI++)
	{
	    histo[iI+*dimI*(iC+*dimC*iY)]=
		histo[2*iI+2**dimI*(2*iC+2**dimC*2*iY)]+
		histo[2*iI+1+2**dimI*(2*iC+2**dimC*2*iY)]+
		histo[2*iI+2**dimI*(2*iC+1+2**dimC*2*iY)]+
		histo[2*iI+1+2**dimI*(2*iC+1+2**dimC*2*iY)]+
		histo[2*iI+2**dimI*(2*iC+2**dimC*(2*iY+1))]+
		histo[2*iI+1+2**dimI*(2*iC+2**dimC*(2*iY+1))]+
		histo[2*iI+2**dimI*(2*iC+1+2**dimC*(2*iY+1))]+
		histo[2*iI+1+2**dimI*(2*iC+1+2**dimC*(2*iY+1))];

	    if(histo[iI+*dimI*(iC+*dimC*iY)]>0) color_size++;
	}
	}
	}
    }

/*   Creating the color table and its frequency   */


    *tabla_ICY=(int *) calloc(*dimI**dimC**dimY,sizeof(int));
    freq=(double *) calloc(color_size,sizeof(double));

    it=0;
    for(iY=0;iY<*dimY;iY++)
    {
    for(iC=0;iC<*dimC;iC++)
    {
    for(iI=0;iI<*dimI;iI++)
    {
	if(histo[iI+*dimI*(iC+*dimC*iY)]>0)
	{
	    tabla_ICY[0][it]=iI+*dimI*(iC+*dimC*iY);
	    freq[it++]=(double)histo[iI+*dimI*(iC+*dimC*iY)];
	}
    }   
    }
    }
   

/*   Ordering table according to the frequency value   */

    quicksort_ref(0,color_size-1,tabla_ICY[0],freq);

/*        Freeing memory before finishing              */

    free(histo);
    free(freq);


    return(color_size);
}

void RGB_ICY( int dimI, int dimC, int dimY, int mode, char *Red, char *Green,
	      char *Blue, int *iI, int *iC, int *iY)
{
    double dR,dG,dB;
    double dI,dC,dY;
    double norm;

    if(mode==0)
    {
	dR=((double)*Red)/256.;
	if(dR<0.) dR+=1.;
	dG=((double)*Green)/256.;
	if(dG<0.) dG+=1.;
	dB=((double)*Blue)/256.;
	if(dB<0.) dB+=1.;

	dI=WR*dR+WG*dG+WB*dB;
	dC=0.5*(dG-dR+1.);
	dY=0.25*(dG+dR-2.*dB+2);

	*iI=(int)(dI*(double)dimI);
	if(*iI>dimI) *iI=dimI-1;
	*iC=(int)(dC*(double)dimC);
	if(*iC>dimC) *iC=dimC-1;
	*iY=(int)(dY*(double)dimY);
	if(*iY>dimI) *iY=dimY-1;

    }
    else
    {
	dI=((double)*iI)/((double)dimI);
	dC=2.*((double)*iC)/((double)dimC)-1.;
	dY=4.*((double)*iY)/((double)dimY)-2.;

	norm=2.*(WR+WG+WB);
	dR=(2.*dI+WB*dY-(2.*WG+WB)*dC)/norm;
	dG=(2.*dI+WB*dY+(2.*WR+WB)*dC)/norm;
	dB=(2.*dI-(WR+WG)*dY+(WR-WG)*dC)/norm;

	if(dR<0) dR=0.;
	if(dG<0) dG=0.;
	if(dB<0) dB=0.;

	dR=256*dR;
	dG=256*dG;
	dB=256*dB;
	*Red=(((int)dR)>255)?(char)255:(char)(int)dR;
	*Green=(((int)dG)>255)?(char)255:(char)(int)dG;
	*Blue=(((int)dB)>255)?(char)255:(char)(int)dB;
	
    }


}


void codifica_imagen_con_tabla_ICY( int dimx, int dimy, int dimI, int dimC,
				    int dimY, char **Red, char **Green, 
				    char **Blue, int color_size, 
				    int *tabla_ICY, int *decoded)
{
    int apren;
    int ix,iy,ic;
    int iI,iC,iY,iI0,iC0,iY0;


    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	RGB_ICY(dimI,dimC,dimY,0,&Red[iy][ix],&Green[iy][ix],&Blue[iy][ix],
		&iI0,&iC0,&iY0);

	for(ic=0,apren=0;(ic<color_size)&&(apren==0);ic++)
	{
	    iI=tabla_ICY[ic]%dimI;
	    iC=((tabla_ICY[ic]-iI)/dimI)%dimC;
	    iY=(tabla_ICY[ic]-iI-dimI*iC)/(dimI*dimC);

	    
	    if((iI==iI0)&&(iC==iC0)&&(iY==iY0))
	    {
		apren=1;
		decoded[ix+dimx*(dimy-1-iy)]=ic;
	    }
	    
	}

	if(apren==0)
	    printf("Error: Color no encontrado en la tabla: %d %d\n",ix,iy);
    }
    }
    
    free(tabla_ICY);
}



void graba_foto(int dimx, int dimy, char* nombre, double **datos)
{
    graba_foto_block(dimx,dimy,1.,nombre,datos);
}

void graba_foto_block(int dimx, int dimy, double block, char* nombre,
	double **datos)
{
    char **foto;
    int sizex,sizey;
    int ix,iy,kx,ky;

    sizex=dimx*block;
    sizey=dimy*block;

    foto=reservar_matriz_char(sizey,sizex);
    prepara_foto_block(dimx,dimy,block,datos,foto);
    graba_gif(sizex,sizey,nombre,foto,foto,foto);
    liberar_matriz_char(foto,sizey);

}

void graba_foto_block_limites(int dimx, int dimy, double block, 
	char* nombre, double min_c, double max_c, double **datos)
{
    char **foto;
    int sizex,sizey;
    int ix,iy,kx,ky;

    sizex=dimx*block;
    sizey=dimy*block;
    
    foto=reservar_matriz_char(sizey,sizex);
    prepara_foto_block_limites(dimx,dimy,min_c,max_c,block,datos,foto);
    graba_gif(sizex,sizey,nombre,foto,foto,foto);
    liberar_matriz_char(foto,sizey);

}


void graba_foto_vec_block(int dimx, int dimy, double block, char* nombre,
			  double **vx, double **vy)
{
    double **datos;
    int ix,iy;

    datos=reservar_matriz(dimy,dimx);
    
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	datos[iy][ix]=sqrt(vx[iy][ix]*vx[iy][ix]+
			   vy[iy][ix]*vy[iy][ix]);
	if(datos[iy][ix]>1e-30) datos[iy][ix]=log(datos[iy][ix]);
	else datos[iy][ix]=-30.*log(10.);
    }
    }
    graba_foto_block(dimx,dimy,block,nombre,datos);

    liberar_matriz(datos,dimy);
}

void graba_foto_4(int dimx, int dimy, char* nombre, double **datos)
{
    char **foto;

    foto=reservar_matriz_char(dimy,dimx);
    prepara_foto_4(dimx,dimy,datos,foto);
    graba_gif(dimx,dimy,nombre,foto,foto,foto);
    liberar_matriz_char(foto,dimy);

}

void graba_mask_block(int dimx, int dimy, double block, char* nombre,
	char **mask)
{
    char **foto;
    int sizex,sizey;
    int ix,iy,ibx,iby;
    int beff,prov,buff;


    sizex=dimx*block;
    sizey=dimy*block;

    foto=reservar_matriz_char(sizey,sizex);
    if(block>=1)
    {
     	beff=(int)block;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		foto[iby+beff*iy][ibx+beff*ix]=mask[iy][ix];
	    }
	    }
	}
	}
    }
    else
    {
	beff=(int)(1./block);
	for(iy=0;iy<sizey;iy++)
	{
	for(ix=0;ix<sizex;ix++)
	{
	    buff=0;
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		prov=(int)mask[iy*beff+iby][ix*beff+ibx];
		if(prov<0) prov+=256;
		buff+=prov;
	    }
	    }
	    buff=buff/(beff*beff);
	    foto[iy][ix]=(char)buff;
	}
	}
	}
    graba_gif(sizex,sizey,nombre,foto,foto,foto);
    liberar_matriz_char(foto,sizey);

}

void graba_RGB_block(int dimx, int dimy, double block, char* nombre,
		     char **Red, char **Green, char **Blue)
{
    char **fR,**fG,**fB;
    int sizex,sizey;
    int ix,iy;

    sizex=dimx*block;
    sizey=dimy*block;
    
    fR=reservar_matriz_char(sizey,sizex);
    fG=reservar_matriz_char(sizey,sizex);
    fB=reservar_matriz_char(sizey,sizex);

    prepara_char_block(dimx,dimy,block,Red,fR);
    prepara_char_block(dimx,dimy,block,Green,fG);
    prepara_char_block(dimx,dimy,block,Blue,fB);

    graba_gif(sizex,sizey,nombre,fR,fG,fB);
    liberar_matriz_char(fR,sizey);
    liberar_matriz_char(fG,sizey);
    liberar_matriz_char(fB,sizey);

}

void graba_foto_color_block( int dimx, int dimy, int n_cr, double block,
			     char *nombre, double ***datos)
{

    char ***foto;
    int dimxf,dimyf;
    int icr;

    dimxf=dimx*block;
    dimyf=dimy*block;

    foto=reservar_tritensor_char(n_cr,dimyf,dimxf);
    
    for(icr=0;icr<n_cr;icr++) 
	prepara_foto_block(dimx,dimy,block,datos[icr],foto[icr]);

    if(n_cr==1) graba_gif(dimxf,dimyf,nombre,foto[0],foto[0],foto[0]);
    else graba_ppm(dimxf,dimyf,1,nombre,foto[0],foto[1],foto[2]);


    liberar_tritensor_char(foto,n_cr,dimyf);

}

void graba_video_block( int dimx, int dimy, int n_cr, double block, 
			int modo, char *nombre, double ***datos)
{

    char ***video;
    int dimxf,dimyf;
    int icr;

    dimxf=dimx*block;
    dimyf=dimy*block;

    video=reservar_tritensor_char(n_cr,dimyf,dimxf);
    if(modo==1)
    {

	for(icr=0;icr<n_cr;icr++)
	    prepara_foto_block(dimx,dimy,block,datos[icr],video[icr]);
	
	if(n_cr==1) graba_gif_animado(dimxf,dimyf,modo,1,nombre,
				      video[0],video[0],video[0]);
	else graba_gif_animado(dimxf,dimyf,modo,1,nombre,
			       video[0],video[1],video[2]);
		
    }
    else graba_gif_animado(dimxf,dimyf,modo,1,nombre,
			   video[0],video[0],video[0]);
    

    liberar_tritensor_char(video,n_cr,dimyf);

}

void graba_video_RGB_block( int dimx, int dimy, double block, int modo, 
			    char *nombre, char **Red, char **Green, 
			    char **Blue)
{
    char ***video;
    int dimxf,dimyf;

    dimxf=dimx*block;
    dimyf=dimy*block;


    video=reservar_tritensor_char(3,dimyf,dimxf);
    if(modo==1)
    {
	prepara_char_block(dimx,dimy,block,Red,video[0]);
	prepara_char_block(dimx,dimy,block,Green,video[1]);
	prepara_char_block(dimx,dimy,block,Blue,video[2]);
	
	graba_gif_animado(dimxf,dimyf,modo,1,nombre,
			  video[0],video[1],video[2]);
    }
    else if(modo==0) graba_gif_animado(dimxf,dimyf,modo,0,nombre,
				       video[0],video[1],video[2]);
    else graba_gif_animado(dimxf,dimyf,modo,1,nombre,
			   video[0],video[1],video[2]);


    liberar_tritensor_char(video,3,dimyf);

}


void prepara_foto_block( int dimx, int dimy, double block, double **datos,
			 char **foto)
{
    char buff;
    double maximo,minimo,prov;
    int xmax,ymax;
    int ix,iy,ibx,iby,beff;

    if(block>=1.)
    {
	beff=(int) block;
	maximo=datos[0][0];
	minimo=datos[0][0];
       	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    maximo=fMax(maximo,datos[iy][ix]);
	    minimo=fMin(minimo,datos[iy][ix]);
	}
	}
	if(maximo-minimo>1e-30)
	{
	    for(iy=0;iy<dimy;iy++)
	    {
	    for(ix=0;ix<dimx;ix++)
	    {
		buff=(char)((int)(255*(datos[iy][ix]-minimo)
				  /(maximo-minimo)));
		for(iby=0;iby<beff;iby++)
		{
		for(ibx=0;ibx<beff;ibx++)
		{
		    foto[iy*beff+iby][ix*beff+ibx]=buff;
		}
		}
	    }
	    }
	}
	else
	{
	    for(iy=0;iy<dimy*beff;iy++)
	    {
	    for(ix=0;ix<dimx*beff;ix++)
	    {
		foto[iy][ix]=0x7f;
	    }
	    }
	}
    }
    else
    {
	beff=(int)(1./block);
	xmax=dimx/beff;
	ymax=dimy/beff;
	maximo=-1e30;
	minimo=1e30;
	for(iy=0;iy<ymax;iy++)
	{
	for(ix=0;ix<xmax;ix++)
	{
	    prov=0.;
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		prov+=datos[iy*beff+iby][ix*beff+ibx];
	    }
	    }
	    maximo=fMax(maximo,prov);
	    minimo=fMin(minimo,prov);
	}
	}
	if(maximo>minimo)
	{
	    for(iy=0;iy<ymax;iy++)
	    {
	    for(ix=0;ix<xmax;ix++)
	    {
		prov=0.;
		for(iby=0;iby<beff;iby++)
		{
		for(ibx=0;ibx<beff;ibx++)
		{
		    prov+=datos[iy*beff+iby][ix*beff+ibx];
		}
		}
		foto[iy][ix]=(char)((int)(255*(prov-minimo)
					  /(maximo-minimo)));
	    }
	    }
	}
	else
	{
	    for(iy=0;iy<ymax;iy++)
	    {
	    for(ix=0;ix<xmax;ix++)
	    {
		foto[iy][ix]=0x7f;
	    }
	    }
	}
    }
    
}

void prepara_foto_fijo_block( int dimx, int dimy, double block, 
			      int levels, double **datos, char **foto)
{
    char buff;
    double prov;
    int xmax,ymax;
    int ix,iy,ibx,iby,beff;

    if(block>=1.)
    {
	beff=(int) block;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    prov=255.*datos[iy][ix]/((double)levels-1.);
	    if(prov>255) buff=(char) 255;
	    else if(prov<0) buff=(char) 0;
	    else buff=(char)(int)prov;
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		foto[iy*beff+iby][ix*beff+ibx]=buff;
	    }
	    }
	}
	}
    }
    else
    {
	beff=(int)(1./block);
	xmax=dimx/beff;
	ymax=dimy/beff;
	for(iy=0;iy<ymax;iy++)
	{
	for(ix=0;ix<xmax;ix++)
	{
	    prov=0.;
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		prov+=datos[iy*beff+iby][ix*beff+ibx];
	    }
	    }
	    prov=255.*prov/((double)(beff*beff*levels));
	    if(prov>255) buff=(char) 255;
	    else if(prov<0) buff=(char) 0;
	    else buff=(char)(int)prov;
	    foto[iy][ix]=buff;
	}
	}
    }
    
}

void prepara_foto_block_limites( int dimx, int dimy, double min_c, 
				 double max_c, double block, double **datos, 
				 char **foto)
{
    double **buff;
    int ix,iy;
    
    buff=reservar_matriz(dimy,dimx);
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	if(datos[iy][ix]>max_c) buff[iy][ix]=max_c;
	else if(datos[iy][ix]<min_c) buff[iy][ix]=min_c;
	else buff[iy][ix]=datos[iy][ix];
    }
    }
    prepara_foto_block(dimx,dimy,block,buff,foto);
    liberar_matriz(buff,dimy);

}


void prepara_foto_log( int dimx, int dimy, double **datos, char **foto)
{
    int ix,iy;
    double maximo,minimo,buff;

    maximo=-1e30;
    minimo=1e30;
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	buff=fabs(datos[iy][ix]);
	if(buff>1e-30)
	{
	    maximo=fMax(maximo,log(buff));
	    minimo=fMin(minimo,log(buff));
	}
    }
    }
    printf("%f   %f\n",maximo,minimo);
    if(maximo-minimo>1e-30)
    {
	maximo-=minimo;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    buff=fabs(datos[iy][ix]);
	    if(buff>1e-30) buff=log(buff)-minimo;
	    else buff=0.;
	    foto[iy][ix]=(char)((int)(255*buff/maximo));
	}
	}
    }
    else
    {
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    foto[iy][ix]=0x7f;
	}
	}
    }
}


void prepara_foto_4(int dimx, int dimy, double **datos, char **foto)
{
    double maximo,minimo,prov;
    int ix,iy;
                        
    maximo=datos[0][0];
    minimo=datos[0][0];
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	maximo=fMax(maximo,datos[iy][ix]);
	minimo=fMin(minimo,datos[iy][ix]);   
    }               
    }               
    
     	    
    if(maximo-minimo>1e-30)
    {
	for(iy=0;iy<dimy/2;iy++)
	{
	for(ix=0;ix<dimx/2;ix++)
	{
	    prov=255.*(datos[iy+dimy/2][ix+dimx/2]-minimo)/(maximo-minimo);
	    foto[iy][ix]=(char)((int)(prov));
	    
	    prov=255.*(datos[iy+dimy/2][ix]-minimo)/(maximo-minimo);
	    foto[iy][ix+dimx/2]=(char)((int)(prov));

	    prov=255.*(datos[iy][ix+dimx/2]-minimo)/(maximo-minimo);
	    foto[iy+dimy/2][ix]=(char)((int)(prov));

	    prov=255.*(datos[iy][ix]-minimo)/(maximo-minimo);
	    foto[iy+dimy/2][ix+dimx/2]=(char)((int)(prov));
	}
	}
    }
    else
    {
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    foto[iy][ix]=0x7f;
	}
	}
    }
}                       


void prepara_foto_log_4(int dimx, int dimy, double **datos, char **foto)
{
    double maximo,minimo,prov,buff;
    int ix,iy;

    maximo=-1e30;
    minimo=1e30;
    for(iy=0;iy<dimy;iy++)
    {
    for(ix=0;ix<dimx;ix++)
    {
	buff=fabs(datos[iy][ix]);
	if(buff>1e-30)
	{
	    maximo=fMax(maximo,log(buff));
	    minimo=fMin(minimo,log(buff));
	}
    }               
    }               
 
     	    
    if(maximo-minimo>1e-30)
    {
	maximo-=minimo;
	for(iy=0;iy<dimy/2;iy++)
	{
	    for(ix=0;ix<dimx/2;ix++)
	    {
		buff=fabs(datos[iy+dimy/2][ix+dimx/2]);
		if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
		else prov=0.;
		foto[iy][ix]=(char)((int)(255*prov));
	    }
	    for(ix=dimx/2;ix<dimx;ix++)
	    {
		buff=fabs(datos[iy+dimy/2][ix-dimx/2]);
		if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
		else prov=0.;
		foto[iy][ix]=(char)((int)(255*prov));
	    }
	}

	for(iy=dimy/2;iy<dimy;iy++)
	{
	    for(ix=0;ix<dimx/2;ix++)
	    {
		buff=fabs(datos[iy-dimy/2][ix+dimx/2]);
		if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
		else prov=0.;
		foto[iy][ix]=(char)((int)(255*prov));
	    }
	    for(ix=dimx/2;ix<dimx;ix++)
	    {
		buff=fabs(datos[iy-dimy/2][ix-dimx/2]);
		if(buff>1e-30) prov=(log(buff)-minimo)/maximo;
		else prov=0.;
		foto[iy][ix]=(char)((int)(255*prov));
	    }
	}
    }
    else
    {
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    foto[iy][ix]=0x7f;
	}
	}
    }
}

void prepara_char_block(int dimx, int dimy, double block, char **grayin,
			char **grayout)
{
    int sizex,sizey;
    int ix,iy,ibx,iby;
    int beff,buff,prov;

    sizex=dimx*block;
    sizey=dimy*block;
    
    if(block>=1)
    {
	beff=(int)block;
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		grayout[iby+beff*iy][ibx+beff*ix]=grayin[iy][ix];
	    }
	    }
	}
	}
    }
    else
    {
	beff=(int)(1./block);
	for(iy=0;iy<sizey;iy++)
	{
	for(ix=0;ix<sizex;ix++)
	{
	    buff=0;
	    for(iby=0;iby<beff;iby++)
	    {
	    for(ibx=0;ibx<beff;ibx++)
	    {
		prov=grayin[iy*beff+ibx][ix*beff+ibx];
		if(prov<0) prov+=256;
		buff+=prov;
	    }
	    }
	    buff=buff/(beff*beff);
	    grayout[iy][ix]=(char)buff;
	}
	}
	
    }

}

int graba_pgm( int dimx, int dimy, int bin, char *nombre_out, char **cont)
{
    FILE* canal;
    char cab_raw[46]=
	"P5\n# Creado por un programa de Antonio Turiel\n";
    char cab_asc[46]=
	"P2\n# Creado por un programa de Antonio Turiel\n";
    int ix,iy;
    int dato;


    if(bin)
    {
	canal=fopen(nombre_out,"wb");
	fwrite(cab_raw,sizeof(char),46,canal);
	fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
	for(iy=0;iy<dimy;iy++) fwrite(cont[dimy-1-iy],sizeof(char),dimx,canal);
    }
    else
    {
	canal=fopen(nombre_out,"wt");
	fwrite(cab_asc,sizeof(char),46,canal);
	fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    dato=(int) cont[dimy-1-iy][ix];
	    if(dato<0) dato=256+dato;
	    fprintf(canal,"%d\n",dato);
	}
	}
    }
    fclose(canal);
    return(0);
}


int graba_ppm( int dimx, int dimy, int bin, char *nombre_out, char **Red, 
	char **Green, char **Blue)
{
    FILE* canal;
    char cab_raw[46]=
	"P6\n# Creado por un programa de Antonio Turiel\n";
    char cab_asc[46]=
	"P3\n# Creado por un programa de Antonio Turiel\n";
    int ix,iy;
    int datoR,datoG,datoB;
    int salida;
         

    if(bin)
    {
	canal=fopen(nombre_out,"wb");
	fwrite(cab_raw,sizeof(char),46,canal);
	fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    fwrite(&Red[dimy-1-iy][ix],sizeof(char),1,canal);
	    fwrite(&Green[dimy-1-iy][ix],sizeof(char),1,canal);
	    fwrite(&Blue[dimy-1-iy][ix],sizeof(char),1,canal);
	}
	}
    }
    else
    {
	canal=fopen(nombre_out,"wt");
	fwrite(cab_asc,sizeof(char),46,canal);
	fprintf(canal,"%d %d\n%d\n",dimx,dimy,255);
	
	for(iy=0;iy<dimy;iy++)
	{
	for(ix=0;ix<dimx;ix++)
	{
	    datoR=(int) Red[dimy-1-iy][ix];
	    if(datoR<0) datoR=256+datoR;
	    datoG=(int) Green[dimy-1-iy][ix];
	    if(datoG<0) datoG=256+datoG;
	    datoB=(int) Blue[dimy-1-iy][ix];
	    if(datoB<0) datoB=256+datoB;
	    fprintf(canal,"%d  %d  %d\n",datoR,datoG,datoB);
	}
	}
    }
    fclose(canal);
    return(0);
}

int graba( int dimx, int dimy, char* nombre_out,double **data)
{
    FILE* canal;
    int ix,iy;

    canal=fopen(nombre_out,"wb");
    for(iy=0;iy<dimy;iy++) fwrite(data[iy],sizeof(double),dimx,canal);
    fclose(canal);
    return(0);
}

