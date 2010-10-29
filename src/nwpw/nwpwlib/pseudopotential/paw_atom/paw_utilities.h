#ifndef _PAW_GET_WORD_H_
#define _PAW_GET_WORD_H_
/*
   $Id$
*/


#define	NIL	((char *) EOF)

extern char	*paw_get_word(FILE *stream);
extern int paw_find_word(char* word, FILE *fp);
extern char*  paw_spd_Name(int );
extern void paw_lu_decompose(int N, double **c, double **a, double **b);
extern void paw_test_lu_decompose();
extern void paw_test_triang_matrix_inverse();
extern void paw_triang_matrix_inverse(char* matrix_type, int N, double **a, double **a_inv);
extern void paw_square_matrix_product(int n, double **a,double **b,double **c);
extern void paw_print_matrix(int n, double **a, char* title);
#endif


