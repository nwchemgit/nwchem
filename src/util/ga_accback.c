#if (defined(CRAY) || defined(ARDENT) || defined(WIN32))&& !defined(__crayx1)&&!defined(__MINGW32__)
#   define ga_access_callback_release_ GA_ACCESS_CALLBACK_RELEASE 
#endif

#include "ga.h"
#include "macdecls.h"
#if defined(CRAY)&& !defined(__crayx1) 
#define FATR
#include <fortran.h> /* Required for Fortran-C string interface on Crays */
#endif /* CRAY */
#ifdef WIN32
#include "typesf2c.h"
#endif
#ifdef USE_FAPI
#  define COPYC2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[i]=(Integer)(carr)[i];} 
#  define COPYF2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[i]=(int)(farr)[i];} 
#  define COPYF2C_64(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[i]=(int64_t)(farr)[i];} 
#  define COPYINDEX_F2C     COPYF2C
#  define COPYINDEX_F2C_64  COPYF2C_64
#  define COPYINDEX_C2F     COPYC2F
#else
#  define COPYC2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i];} 
#  define COPYF2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int)(farr)[i];} 
#  define COPYF2C_64(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int64_t)(farr)[i];} 
#  define COPYINDEX_C2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i]+1;}
#  define COPYINDEX_F2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int)(farr)[i] -1;}
#  define COPYINDEX_F2C_64(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int64_t)(farr)[i] -1;}
#define BASE_0
#endif

/**
\ingroup util_ga
@{
*/
			     
/*\ PROVIDE ACCESS TO A PATCH OF A GLOBAL ARRAY WITH CALLBACK AND RELEASE
\*/
void FATR ga_access_callback_release_(g_a, ilo, ihi, jlo, jhi, 
				     callback, 
				     arg1, arg2, arg3, arg4, arg5, arg6, arg7)
     Integer *g_a, *ilo, *ihi, *jlo, *jhi;
     Integer (*callback)(Integer *,Integer *,Integer *,Integer *,Integer *,
			 void *, Integer*, 
			 void *, void *, void *, void *, void *, void *, void *);
     void *arg1, *arg2, *arg3, *arg4, *arg5, *arg6, *arg7;
{
  Integer ndim=GA_Ndim(*g_a), lo[2], hi[2], ld[2],
          result; /* Fortran variables */
  int alo[2], ahi[2], ald[2], ag_a; /* variables for the C-interfaces */
  void *ptr;

  if(ndim != 2) 
    GA_Error("ga_access: 2D API cannot be used for array dimension",ndim);

  ag_a=*g_a;
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  COPYINDEX_F2C(lo,alo,ndim);
  COPYINDEX_F2C(hi,ahi,ndim);
  NGA_Access(ag_a,alo,ahi,&ptr,ald); /* This routine sets ald[] */
  ld[0]=ald[0];
  ld[1]=ald[1];
  result = callback(g_a,ilo,ihi,jlo,jhi,ptr,&ld[0],arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  if (result) {
    NGA_Release_update(ag_a, alo, ahi);
  } else {
    NGA_Release(ag_a, alo, ahi);
  }
} 

void FATR nga_access_callback_release_(g_a, ilo, ihi,
				     callback, 
				     arg1, arg2, arg3, arg4, arg5, arg6, arg7)
     Integer *g_a, ilo[], ihi[];
     Integer (*callback)(Integer *,Integer *,Integer *,
			 void *, Integer*, 
			 void *, void *, void *, void *, void *, void *, void *);
     void *arg1, *arg2, *arg3, *arg4, *arg5, *arg6, *arg7;
{
  Integer ndim=GA_Ndim(*g_a), ild[GA_MAX_DIM],
          result; /* Fortran variables */
  int alo[GA_MAX_DIM], ahi[GA_MAX_DIM], ald[GA_MAX_DIM],
      ag_a; /* variables for the C-interfaces */
  int ii;
  void *ptr;

  ag_a=*g_a;
  COPYINDEX_F2C(ilo,alo,ndim);
  COPYINDEX_F2C(ihi,ahi,ndim);
  for (ii = 0; ii < GA_MAX_DIM; ii++) ald[ii] = 0;
  for (ii = 0; ii < GA_MAX_DIM; ii++) ild[ii] = 0;
  NGA_Access(ag_a,alo,ahi,&ptr,ald); /* This routine sets ald[] */
  COPYC2F(ald,ild,ndim-1);
  result = callback(g_a,ilo,ihi,ptr,ild,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  if (result) {
    NGA_Release_update(ag_a, alo, ahi);
  } else {
    NGA_Release(ag_a, alo, ahi);
  }
} 

/**
@}
*/

/* $Id$ */
