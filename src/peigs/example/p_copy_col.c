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
#include <stdio.h>

#include "globalp.c.h"


/*
  copy columns of a lower triangular matrix A on this processor into vecD 
  */

void p_copy_col( n, mapA, vecA, vecD ) 
     Integer *n, *mapA;
     DoublePrecision **vecA, **vecD;
{
  static Integer IONE = 1;
  Integer i, me;
  Integer msize;
  Integer indx, size;
  extern Integer mxmynd_ ();
  extern void dcopy_ ();
  extern char *memcpy();
  
  msize = *n;
  me = mxmynd_ ();
  
  fprintf(stderr, " in p_copy_col me = %d \n", me);
  
  msize = *n;
  indx = 0;
  for ( i = 0; i < msize; i++ ) {
    fprintf(stderr, " me = %d just before if mapA %d = %d  \n", me, i, mapA[i]);
    if ( mapA[i] ==  me ) {
      size = msize - i;
      dcopy_ (&size, vecA[indx], &IONE, vecD[indx], &IONE);
      indx++;
    }
  }
  
  fprintf(stderr, " out p_copy_col me = %d \n", me);
  return;
}

