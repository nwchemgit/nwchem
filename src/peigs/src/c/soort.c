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
 *  PeIGS utility routine: qqsort
 *
 *  K & R's quick sort
 *  
 */

#include <stdio.h>
#include <memory.h>

#include "globalp.c.h"

#define GG_SWAP(a,b,t) {t=a;a=b;b=t;}

Integer qqsort(v,left,right)
     Integer *v,left,right;
{
  
  Integer i,last,tmp;
  Integer qqsort();
  
  if (left >= right)
    return(0);
  
  i = (left + right)/2;
  GG_SWAP(v[left],v[i],tmp);
  
  last = left;
  for ( i=left+1; i<=right; i++ ) {
    if (v[i] < v[left] ) {
      last++;
      GG_SWAP(v[last],v[i],tmp);
    }
  }
  GG_SWAP(v[left],v[last],tmp);
  qqsort(v,left,last-1);
  qqsort(v,last+1,right);
  return(0);
}











