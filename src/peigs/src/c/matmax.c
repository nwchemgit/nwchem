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
#include <math.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

void mat_max ( n, m, colA,  mapA, norm, iwork, work)
     Integer *n, *m, *mapA, *iwork;
     DoublePrecision **colA, *norm, *work;
{
  /*
    computes the max. absolute element of matrix A
    
    n       = (input/integer) vector length
    m       = (input/integer) number of vectors
    
    colA    = (input/1-D array of pointers to doubles ) pointer to the columns
    
    mapA    = (input/1-D array of integers) mapA[i] = the processor
    id which owns column i
    
    iwork   = (scratch space/integer) length  nvecsA + m
    
    work    = DoublePrecision precision scratch space of size
              mxlbuf_() / sizeof( DoublePrecision ) + 1
    
    */
  
  static Integer IONE = 1;
  Integer ll, nprocs, i, me, nvecsA, *mapvecA;
  Integer *proclist;
  Integer *iscrat;
  DoublePrecision t, anorm;
  
  extern DoublePrecision damax_();

  extern void gmax00();
  extern Integer mxmynd_();
  extern Integer reduce_list2();
  extern Integer fil_mapvec_();
  
  me  = mxmynd_();
  
  ll = *m;
  iscrat   = iwork;
  mapvecA  = iscrat;
  nvecsA   = fil_mapvec_( &me, &ll, mapA, mapvecA);
  iscrat  += nvecsA;
  
  if ( nvecsA <= 0 )
    return;
  
  proclist = iscrat;
  nprocs   = reduce_list2( ll, mapA, proclist);
  iscrat  += nprocs;
  
  anorm = 0.0e0;
  for ( i = 0; i < nvecsA; i++ ) {
    t = damax_( n, colA[i], &IONE);
    anorm = max(anorm, t);
  }
  
  /*
    global combine max
    */

  gmax00( (char *) &anorm, 1, 5, 16, proclist[0], nprocs, proclist, work);
  
  *norm = anorm;
  return;
}


