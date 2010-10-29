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
  A simple version of quicksort; hacked from Kernighan and Ritchie, page 87
  may be faster to use a 'real data structure' but ...
  */

static void SWAPP(a,b,t) 
     Integer *a, *b, *t;
{
  t[0] = a[0];
  t[1] = a[1];
  t[2] = a[2];
  t[3] = a[3];
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
  a[3] = b[3];
  b[0] = t[0];
  b[1] = t[1];
  b[2] = t[2];
  b[3] = t[3];
  
  return;
}


Integer iqsort(v,left,right)
     Integer **v, left, right;
     /*
       quick sort in increasing order; from v[left] to v[right]
       */
{
  Integer i,last,tmp[4];
  Integer leftsz, rightsz;
  
  if ( left >= right ) return 0;
  SWAPP(v[left],v[(left+right)/2],tmp);
  last = left;
  for ( i=left+1; i<=right; i++ ) {
    leftsz = v[i][1] - v[i][0] +1;
    rightsz = v[left][1] - v[left][0] + 1;
    if (leftsz < rightsz ) {
      last++;
      SWAPP(v[last],v[i],tmp);
    }
  }
  SWAPP(v[left],v[last],tmp);
  iqsort(v,left,last-1);
  iqsort(v,last+1,right);
  return 0;
}


