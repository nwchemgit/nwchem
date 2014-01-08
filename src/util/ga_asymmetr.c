/* $Id$ */
/**
 * Symmetrizes matrix A:  A := .5 * (A+A`)
 * diag(A) remains unchanged
 *    Antisymmetrizes matrix A:  A := .5 * (A-A`)
 *    copied from ga_symmetrize, A+A' replaced by A-A'

 * 
 */

#include "ga.h"
#include "macdecls.h"

/**
\ingroup util_ga
@{
*/

void FATR
gai_subtr(int *lo, int *hi, void *a, void *b, DoublePrecision alpha,
          int type, Integer nelem, int ndim) {

  Integer i, j, m=0;
  Integer nrow, ncol, indexA=0, indexB=0;
  DoublePrecision *A = (DoublePrecision *)a, *B = (DoublePrecision*)b;
  Integer offset1=1, offset2=1;

  nrow = hi[1] - lo[1] + 1;
  ncol = hi[0] - lo[0] + 1;

  for(i=2; i<ndim; i++) {
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
  int i, me = GA_Nodeid();
  extern void * FATR ga_malloc(Integer nelem, int type, char *name);
  extern void FATR ga_free(void *ptr);
  void FATR gai_subtr(int *lo, int *hi, void *a, void *b, DoublePrecision alpha,
                      int type, Integer nelem, int ndim);

  int alo[GA_MAX_DIM], ahi[GA_MAX_DIM], lda[GA_MAX_DIM];
  int blo[GA_MAX_DIM], bhi[GA_MAX_DIM], ldb[GA_MAX_DIM];
  int ndim, dims[GA_MAX_DIM], type;
  Integer nelem=1;
  Logical have_data;
  void *a_ptr, *b_ptr;

  GA_Sync();


  
  NGA_Inquire((int)(*g_a), &type, &ndim, dims);
  
  if (dims[0] != dims[1]) 
    GA_Error("ga_sym: can only sym square matrix", 0L);
  
  /* Find the local distribution */
  NGA_Distribution((int)(*g_a), me, alo, ahi);
 
 
  have_data = ahi[0]>=0;
  for(i=1; i<ndim; i++) have_data = have_data && ahi[i]>=0;
  
  if(have_data) {
    NGA_Access((int)(*g_a), alo, ahi, &a_ptr, lda); 
    
    for(i=0; i<ndim; i++) nelem *= ahi[i]-alo[i] +1;
    b_ptr = (void *) ga_malloc(nelem, MT_C_DBL, "v");
    
    for(i=2; i<ndim; i++) {bhi[i]=ahi[i]; blo[i]=alo[i]; }
    
    /* switch rows and cols */
    blo[1]=alo[0];
    bhi[1]=ahi[0];
    blo[0]=alo[1];
    bhi[0]=ahi[1];

    for (i=0; i < ndim-1; i++) 
      ldb[i] = bhi[i+1] - blo[i+1] + 1; 
    NGA_Get((int)(*g_a), blo, bhi, b_ptr, ldb);
  }
  GA_Sync(); 

  if(have_data) {
    gai_subtr(alo, ahi, a_ptr, b_ptr, alpha, type, nelem, ndim);
    NGA_Release_update((int)(*g_a), alo, ahi);
    ga_free(b_ptr);
  }
  GA_Sync();
}
/**
@}
*/
