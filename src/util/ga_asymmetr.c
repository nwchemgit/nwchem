/* $Id: ga_asymmetr.c,v 1.3 2004-12-08 23:49:10 bert Exp $ */
/**
 * Symmetrizes matrix A:  A := .5 * (A+A`)
 * diag(A) remains unchanged
 *    Antisymmetrizes matrix A:  A := .5 * (A-A`)
 *    copied from ga_symmetrize, A+A' replaced by A-A'

 * 
 */

#include "global.h"
#include "macdecls.h"

void FATR
gai_subtr(Integer *lo, Integer *hi, Void *a, Void *b, DoublePrecision alpha,
        Integer type, Integer nelem, Integer ndim) {

  Integer i, j, m=0;
  Integer nrow, ncol, indexA=0, indexB=0;
  DoublePrecision *A = (DoublePrecision *)a, *B = (DoublePrecision*)b;
  Integer offset1=1, offset2=1;

  nrow = hi[ndim-2] - lo[ndim-2] + 1;
  ncol = hi[ndim-1] - lo[ndim-1] + 1;

  for(i=0; i<ndim-2; i++) {
    offset1 *= hi[i] - lo[i] + 1;
    offset2 *= hi[i] - lo[i] + 1;
  }
  offset1 *= nrow;

  for(j=0; j<offset2; ++j,indexA=j,indexB=j,m=0) {
    for(i=0; i<nrow*ncol; i++, indexA += offset1, indexB += offset2) {
      if(indexA >= nelem) indexA = j + ++m*offset2;
      A[indexA] = alpha *(A[indexA] - B[indexB]);
    }
  }
}

void FATR 
ga_antisymmetrize_(Integer *g_a) {
  
  DoublePrecision alpha = 0.5;
  Integer i, me = ga_nodeid_();
  extern void * FATR ga_malloc(Integer nelem, int type, char *name);
  extern void FATR ga_free(void *ptr);
  void FATR gai_add(Integer *lo, Integer *hi, Void *a, Void *b, DoublePrecision alpha, Integer type, Integer nelem, Integer ndim);

  Integer alo[GA_MAX_DIM], ahi[GA_MAX_DIM], lda[GA_MAX_DIM], nelem=1;
  Integer blo[GA_MAX_DIM], bhi[GA_MAX_DIM], ldb[GA_MAX_DIM];
  Integer ndim, dims[GA_MAX_DIM], type;
  Logical have_data;
  Void *a_ptr, *b_ptr;

  ga_sync_();


  
  nga_inquire_internal_(g_a, &type, &ndim, dims);
  
  if (dims[ndim-1] != dims[ndim-2]) 
    ga_error("ga_sym: can only sym square matrix", 0L);
  
  /* Find the local distribution */
  nga_distribution_(g_a, &me, alo, ahi);
 
 
  have_data = ahi[0]>0;
  for(i=1; i<ndim; i++) have_data = have_data && ahi[i]>0;
  
  if(have_data) {
    nga_access_ptr(g_a, alo, ahi, &a_ptr, lda); 
    
    for(i=0; i<ndim; i++) nelem *= ahi[i]-alo[i] +1;
    b_ptr = (Void *) ga_malloc(nelem, MT_F_DBL, "v");
    
    for(i=0; i<ndim-2; i++) {bhi[i]=ahi[i]; blo[i]=alo[i]; }
    
    /* switch rows and cols */
    blo[ndim-1]=alo[ndim-2];
    bhi[ndim-1]=ahi[ndim-2];
    blo[ndim-2]=alo[ndim-1];
    bhi[ndim-2]=ahi[ndim-1];

    for (i=0; i < ndim-1; i++) 
      ldb[i] = bhi[i] - blo[i] + 1; 
    nga_get_(g_a, blo, bhi, b_ptr, ldb);
  }
  ga_sync_(); 

  if(have_data) {
    gai_subtr(alo, ahi, a_ptr, b_ptr, alpha, type, nelem, ndim);
    nga_release_update_(g_a, alo, ahi);
    ga_free(b_ptr);
  }
  ga_sync_();
}
