#ifndef _PAW_MY_MEMORY_H_
#define _PAW_MY_MEMORY_H_
/*
   $Id: paw_my_memory.h,v 1.2 2004-10-14 22:05:03 bylaska Exp $
*/


extern double   *paw_alloc_1d_array(int);
extern double   **paw_alloc_2d_array(int,int);
extern void     paw_dealloc_1d_array(double*);
extern void     paw_dealloc_2d_array(int,int,double**);


#endif

