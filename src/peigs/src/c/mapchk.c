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
#include "globalp.c.h"

Integer mapchk_( n, map )

  Integer      *n, *map;
     
/*
 *  Routine to determine if map[0:n-1] is a valid list of processor
 *  ids.
 *
 *  Return 0 if 0 <= map[i] < mxnprc_(), i = 0 to n-1
 *  Otherwise return 1, i.e., an invalid map.
 *  
 */
     
{
  Integer             indx, nmax;
  Integer             *iptr;
  
  extern Integer      mxnprc_();
  

  nmax = mxnprc_() - 1;
    
  iptr = map;
  for (indx = 0; indx < *n; indx++ )

     if( *iptr < 0  || *iptr > nmax )
       return(1);

     iptr++;
      
  return(0);
}
