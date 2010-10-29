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
#include <string.h>
#include <math.h>
#include <memory.h>

#include "globalp.c.h"

void pdspev_(n, matrixA, mapA, matZ, mapZ, eval, iscratch, iscsize,
             dblptr, ibuffsz, scratch, ssize, info)
     Integer            *n, *mapA, *mapZ, *iscratch, *iscsize, *ibuffsz,
                     *ssize, *info;
     DoublePrecision         *matrixA, *matZ, *eval, *scratch, *dblptr;
{
  
/*
 *----------------------------------------------------------------------
 *
 * FORTRAN wrapper for pdspev.
 *
 * Code assumes the user is passing a DoublePrecision precision array of length
 * ibuffsz as the 10th argument (dblptr) to this subprogram.  It
 * is assumed that we can just "convert" this DoublePrecision precision FORTRAN
 * array to a C DoublePrecision pointer to DoublePrecision by just casting it as such.
 *
 * It is assumed that matrixA is a DoublePrecision precision 1-D array
 * holding the portions of A owned by this processor in packed storage.
 *
 * It is assumed that matZ is a DoublePrecision precision 1-D array.
 * On exit it will hold the eigenvectors owned by this processor in
 * packed storage.
 *
 *
 * This routine does minimal error checking.
 *----------------------------------------------------------------------
 */
  
  
  /*
   * Local variables
   * ---------------
   */
  
  Integer             i, me, indx, nvecs, nvecsA, nvecsZ, ibuffsize, linfo;
  
  char            msg[ 25 ];
  
  DoublePrecision        **buff_ptr, **vecA, **vecZ;
  
  /*
   * External procedures
   * -------------------
   */
  
  extern Integer      mxmynd_();
  
  extern Integer      ci_entry();
  extern Integer      count_list();
  extern Integer fil_mapvec_();
  extern void     l_exit_();
  
  extern void     pdspev();
  extern void     mxinit_();

  /*
   * Executable code
   * ---------------
   */
  

  
  strcpy( msg,  "Error in pdspev_f." );
  mxinit_();
  
  me = mxmynd_();
  
  /*
   *     Test the input parameters.
   */
  
  linfo = 0;
  
  if ( n == NULL )
    linfo = -1;
  else if ( matrixA == NULL )
    linfo = -2;
  else if ( mapA == NULL )
    linfo = -3;
  else if ( matZ == NULL )
    linfo = -4;
  else if ( mapZ == NULL )
    linfo = -5;
  else if ( eval == NULL )
    linfo = -6;
  else if ( iscratch == NULL )
    linfo = -7;
  else if ( iscsize == NULL )
    linfo = -8;
  else if ( dblptr == NULL )
    linfo = -9;
  else if ( ibuffsz == NULL )
    linfo = -10;
  else if ( scratch == NULL )
    linfo = -11;
  else if ( ssize == NULL )
    linfo = -12;
  else if ( info == NULL )
    linfo = -13;
  
  
  if ( linfo != 0 ) {
    if ( info != NULL )
      *info = linfo;
    
    fprintf( stderr, " %s me = %d argument %d is a pointer to NULL. \n",
             msg, me, -linfo );
    l_exit_( &linfo, msg );
    return;
  }
  
  nvecsA = count_list( me, mapA, n );
  nvecsZ = count_list( me, mapZ, n );
  
  if ( *ibuffsz < nvecsA + nvecsZ + 2 ) {
      *info = -10;
      
      fprintf( stderr, " %s me = %d ibuffsz is too small. \n", msg, me );
      l_exit_( info, msg );
      return;
    }
  
  buff_ptr  = (DoublePrecision **) dblptr;
  ibuffsize = *ibuffsz;
  
  /*
   * Initialize DoublePrecision pointers to DoublePrecision workspace.
   */
  
  vecA = buff_ptr;
  nvecs = fil_mapvec_( &me, n, mapA, iscratch );
  i = 0;
  for (indx = 0; indx < nvecs; indx++){
    vecA[indx] = &matrixA[i];
    i += ( *n - iscratch[indx]);
  }
  
  buff_ptr  += (nvecs + 1);
  ibuffsize -= (nvecs + 1);
  
  vecZ = buff_ptr;
  
  nvecs =count_list( me, mapZ, n );
  for (indx = 0; indx < nvecs; indx++){
    vecZ[indx] = &matZ[indx * *n];
  }
  
  buff_ptr  += (nvecs + 1);
  ibuffsize -= (nvecs + 1);
  
#ifdef DEBUG1
  fprintf(stderr, " just before pdspev %d \n", me);
  fflush(stderr);
#endif
  
  pdspev(n, vecA, mapA, vecZ, mapZ, eval, iscratch, iscsize,
	 buff_ptr, &ibuffsize, scratch, ssize, info);

#ifdef DEBUG1
  fprintf(stderr, " just after pdspev %d \n", me);
  fflush(stderr);
#endif
  
  return;
}





