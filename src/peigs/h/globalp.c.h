/*
 $Id: globalp.c.h,v 1.3 1999-07-28 00:39:18 d3e129 Exp $

 file globalp.c.h 
 */

#ifdef  STD_INT
typedef int    Integer;
#else
typedef long   Integer;
#endif

#ifdef  STD_DBL
typedef   double         DoublePrecision;
#else
typedef   long double    DoublePrecision;
#endif

#include "blas_lapack.h"
#include "peigs_types.h"

