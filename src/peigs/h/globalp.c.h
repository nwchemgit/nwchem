/*
 $Id: globalp.c.h,v 1.4 2000-02-28 21:56:02 d3g270 Exp $

 file globalp.c.h 
 */

#ifdef  STD_INT
typedef int    Integer;
typedef unsigned int unInteger;
#else
typedef long   Integer;
typedef unsigned long unInteger;
#endif

#ifdef  STD_DBL
typedef   double         DoublePrecision;
#else
typedef   long double    DoublePrecision;
#endif

#include "blas_lapack.h"
#include "peigs_types.h"

