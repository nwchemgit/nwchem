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
  PeIGS utility routine
  shell sort given in Harbison and Steele 3rd ed.
  page. 227 of "C, A reference manual," Tartan Laboratories
  
  */

void gshellsort_ (n, v)
     Integer *n, v[];
{
  Integer gap, i, j, temp;
  Integer ndim; 

  ndim = *n;
  gap = 1;
  do ( gap = 3*gap+1); while (gap <= ndim);
  for ( gap /= 3; gap > 0; gap /= 3 )
    for ( i = gap; i < ndim; i++ ) {
      temp = v[i];
      for ( j = i-gap; (j>=0)&&(v[j]>temp); j-=gap)
	v[j+gap] = v[j];
      v[j+gap] = temp;
    }
  return;
}
void dshellsort2_(n, v, indx )
     Integer         *n, *indx;
     DoublePrecision v[]; 
{
  Integer gap, i, j, itemp;
  Integer ndim; 
  DoublePrecision temp;

  ndim = *n;
  gap = 1;
  do ( gap = 3*gap+1); while (gap <= ndim);
  for ( gap /= 3; gap > 0; gap /= 3 )
    for ( i = gap; i < ndim; i++ ) {
      temp = v[i];
      itemp = indx[i];
      for ( j = i-gap; (j>=0)&&(v[j]>temp); j-=gap) {
	v[j+gap] = v[j];
	indx[j+gap] = indx[j];
      }
      v[j+gap] = temp;
      indx[j+gap] = itemp;
    }
  return;
}
void dshellsort_(n, v )
     Integer         *n;
     DoublePrecision v[]; 
{
  Integer gap, i, j;
  Integer ndim; 
  DoublePrecision temp;

  ndim = *n;
  gap = 1;
  do ( gap = 3*gap+1); while (gap <= ndim);
  for ( gap /= 3; gap > 0; gap /= 3 )
    for ( i = gap; i < ndim; i++ ) {
      temp = v[i];
      for ( j = i-gap; (j>=0)&&(v[j]>temp); j-=gap) {
	v[j+gap] = v[j];
      }
      v[j+gap] = temp;
    }
  return;
}
