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

#ifdef CRAY_T3D
/*
 *
 * Remap Fortran names to upper case
 * and remove trailing underscore
 *
 */

#define xstop_	XSTOP
#endif

#include "blas_lapack.h"
#include "peigs_types.h"

