/*
 $Id$
 *======================================================================
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with Government
 * support under Contract Number DE-AC06-76RLO-1830 awarded by the United
 * States Department of Energy.  The Government retains a paid-up
 * non-exclusive, irrevocable worldwide license to reproduce, prepare
 * derivative works, perform publicly and display publicly by or for the
 * Government, including the right to distribute to other Government
 * contractors.
 *
 *======================================================================
 *
 *  -- PEIGS  routine (version 2.1) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */
/*
  Internal PeIGS routine which computes 
  
  Q.W where Q is n1 by n2 matrix and W is an n2 by m matrix
  
  the result is stored in Z

*/


#include <stdio.h>

#include "globalp.c.h"

#define min(a,b) ((a) < (b) ? (a) : (b))

void mxm25_( n1, n2, rowQ, mapQ, m, colW, mapW, colZ, iwork, work)
     Integer *n1, *n2, *mapQ, *m, *mapW, *iwork;
     DoublePrecision *rowQ, *colW, *colZ, *work;
{
  /*
    this is the wrapper for the matrix multiplication routine: Q.W
    
    the matrix Q and W are assume to be in a column wrapped format
    
    */
  
  Integer me, nvecsQ, nvecsW, nvecsZ;
  Integer *iscrat, i;
  
  DoublePrecision **matQ, **matW, **matZ;
  DoublePrecision *dscrat;
  
  extern Integer mxmynd_();
  extern Integer fil_mapvec_();
  extern void mxm25();
  extern Integer count_list();
  
  
  me = mxmynd_();
  nvecsQ = count_list( me, mapQ, n1);
  nvecsW = count_list( me, mapW, m);
  nvecsZ = nvecsW;
  
  matQ =(DoublePrecision **) work;  
  matW =(DoublePrecision **)( work + nvecsQ);
  matZ =(DoublePrecision **)( work + nvecsQ + nvecsW);
  
  dscrat =(DoublePrecision *)( work + nvecsQ + nvecsW + nvecsZ );
  
  for( i = 0; i < nvecsQ; i++ )
    matQ[i] = &rowQ[i * *n2];
  
  for( i = 0; i < nvecsW; i++ )
    matW[i] = &colW[i * *n2];
  
  for( i = 0; i < nvecsW; i++ )
    matZ[i] = &colZ[i * *n1];
  
  iscrat  = iwork;
  
  /*
    use mxm*.c( depending on your data distribution )
    */
  
  mxm25( n1, n2, matQ, mapQ, m, matW, mapW, matZ, iscrat, dscrat );

/*
  for( i = 0; i < nvecsW; i++ )
    for( j = 0; j < *n1; j++ )
      fprintf(stderr, " i = %d j = %d Z = %g \n", j, i, matZ[i][j]);
*/
  
  return;
}



