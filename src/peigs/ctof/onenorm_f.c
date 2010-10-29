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

void one_nrm_( n, m, matA,  mapA, norm, iwork, work)
     Integer *n, *m, *mapA, *iwork;
     DoublePrecision *matA, *norm, *work;
{
  
  Integer          i, nvecs, me;

  DoublePrecision *scratch;
  DoublePrecision **colA;
  
  extern Integer mxmynd_();
  extern void    one_nrm();



  me = mxmynd_();

  colA = (DoublePrecision ** ) work;

  nvecs = -1;
  for( i = 0; i < *m; i++ ) {
    if( mapA[i] == me ) {
       nvecs++;
       colA[nvecs] = &matA[*n * nvecs];
    }
  }
  nvecs++;
  
  scratch = (DoublePrecision * ) work + nvecs + 1;

  
  one_nrm( n, m, colA, mapA, norm, iwork, scratch );
  

  return;
}

