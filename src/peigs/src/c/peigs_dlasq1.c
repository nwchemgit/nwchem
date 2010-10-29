/*
 $Id$
*/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "globalp.c.h"

void peigs_dlasq1( Integer n, DoublePrecision *dplus, DoublePrecision *lplus, DoublePrecision *eval, DoublePrecision *work, Integer *info)
{   
  Integer i, j, iii, me, nproc;
  DoublePrecision *dptr;
  extern void dlasq1_();
  FILE *file;
  char filename[40];
  extern int close();
  
  me    = mxmynd_();
  nproc = mxnprc_();
  
  for ( i = 0; i < n; i++ )
    work[i]=sqrt(dplus[i]);
  
  dptr = &work[n];
  for(i = 0;i < n-1 ;i++){
    dptr[i]=lplus[i] * work[i];
  }
  
  dlasq1_( &n, &work[0], &work[n], &work[n+n], info );
  
  if ( *info != 0 ){
    printf(" error in dlasq1 info = %d \n", *info );
    sprintf( filename, "pdspevx.%d", me);
    file = fopen(filename, "w");
    for ( iii = 0; iii < n; iii++)
      fprintf(file, " %d %20.16f %20.16f \n", iii, dplus[iii], lplus[iii]);
    fclose(file);
    fflush(file);
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

