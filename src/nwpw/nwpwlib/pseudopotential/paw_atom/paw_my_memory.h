#ifndef _PAW_MY_MEMORY_H_
#define _PAW_MY_MEMORY_H_
/*
   $Id$
*/


extern double   *paw_alloc_1d_array(int);
extern double   **paw_alloc_2d_array(int,int);
extern void     paw_dealloc_1d_array(double*);
extern void     paw_dealloc_2d_array(int,int,double**);


#endif

