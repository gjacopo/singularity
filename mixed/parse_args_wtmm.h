#ifndef   	PARSE_ARGS_WTMM_H
# define   	PARSE_ARGS_WTMM_H

int parsing_wtmm( int in0, int siflag, int *deflag, 
		  char **olarg, char **olval, char **olexp,
		  float **ptrvar_f, float **ptrval_f, 
		  int **ptrvar_i, int **ptrval_i,
		  int **ptrflag, int *type, int*olnumb);

void parse_arguments(int argc, char *argv[]);

#endif 	    /* !PARSE_ARGS_WTMM_H */
