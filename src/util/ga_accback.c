#include "global.h"
			     
/*\ PROVIDE ACCESS TO A PATCH OF A GLOBAL ARRAY WITH CALLBACK AND RELEASE
\*/
void FATR ga_access_callback_release_(g_a, ilo, ihi, jlo, jhi, 
				     callback, 
				     arg1, arg2, arg3, arg4, arg5, arg6)
     Integer *g_a, *ilo, *ihi, *jlo, *jhi;
     Integer (*callback)(Integer *,Integer *,Integer *,Integer *,Integer *,
			 void *, Integer*, 
			 void *, void *, void *, void *, void *, void *);
     void *arg1, *arg2, *arg3, *arg4, *arg5, *arg6;
{
  Integer lo[2], hi[2],ndim=ga_ndim_(g_a), ld[2];
  void *ptr;

  if(ndim != 2) 
    ga_error("ga_access: 2D API cannot be used for array dimension",ndim);

  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  nga_access_ptr(g_a,lo,hi,&ptr,ld);
  if (callback(g_a,ilo,ihi,jlo,jhi,ptr,&ld[0],arg1,arg2,arg3,arg3,arg5,arg6))
    nga_release_update_(g_a, lo, hi);
  else
    nga_release_(g_a, lo, hi);
} 

