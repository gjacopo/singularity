/*      Version del 2 de Septiembre, 2004   */

#define STRUCT_DEF

/*    Common definitions for char variables   */

char C0=(char)0xff; // Mask for 0
char CP=(char)0x00; // Mask for +1
char CM=(char)0x7f; // Mask for -1

// The prototype of color variable
typedef struct
{
    char Red;
    char Green;
    char Blue;
} color;


//  Function pointer type for routines in io2D
typedef int(*Read2D)(int,int,int,int,int,int,char*,double ***);
typedef int(*Read2D_mask)(int,int,int,int,int,int,char*,double ***, 
			  char **);

