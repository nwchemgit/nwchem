/*
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

