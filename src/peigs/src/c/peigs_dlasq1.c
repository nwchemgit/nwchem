#include <stdio.h>
#include <math.h>
#include "globalp.c.h"

void peigs_dlasq1( Integer n, DoublePrecision *dplus, DoublePrecision *lplus, DoublePrecision *eval, DoublePrecision *work, Integer *info)
{   
  Integer i, j;
  DoublePrecision *dptr;
  extern void dlasq1_();
  
  for ( i = 0; i < n; i++ )
    work[i]=sqrt(dplus[i]);
  
  dptr = &work[n];
  for(i = 0;i < n-1 ;i++){
    dptr[i]=lplus[i] * work[i];
  }
  
  dlasq1_( &n, &work[0], &work[n], &work[n+n], &info );
  
  if ( info != 0 ){
    printf(" error in dlasq1 info = %d \n", info );
    fflush(stdout);
  }
  
  j = n-1;
  for(i = 0;i < n;i++){
    eval[i] = work[j] * work[j];
    --j;
  }
  return;
  
  /*
    dlasq1_( n, work, work+ *n, work+2* *n, info );
    
    j = *n-1;
    if( fabs(psgn - 1.0) < eps )
    for(i = 0;i < *n;i++){
    eval[i] = work[j] * work[j];
    --j;
    }
    else{
    for(i = 0;i < *n;i++){
    eval[i] = -(work[j] * work[j]);
    --j;
    }
    */
}

