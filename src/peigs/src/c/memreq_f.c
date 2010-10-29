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
/* ****************************************
   
   PeIGS utility memreq_ for FORTRAN routines.
   
   returns the scratch memory requirements for the PDSPGV,
   PDSPEV, and PDSPTRI routines.
   
   */

#include <stdio.h>
#include <memory.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

#define BDR sizeof(DoublePrecision)/sizeof(Integer)

void fmemreq_(type, n, mapA, mapB, mapZ, isize, rsize, ptr_size, iscratch )
     Integer *type, *n, *mapA, *mapB, *mapZ, *isize, *rsize, *ptr_size, *iscratch;
{

/*
    this subroutine computes the memory requirements for this processor to
    safely setup the eigensystem problem.  On output, isize, risize, and ptr_size
    are the sizes ( in the respective data types ) of the data needed by the eigensystem
    program.  It performs the basic error checks and assumes that the user has
    input the correct information into type, n, mapA, mapB, mapZ
    
    input:
                 memory for 

   *type = 0    the generalized eigensystem problem (pdspgvx or pdspgv )
         = 1    the standard eigensystem problem (pdspevx or pdspev )
         = 2    the tri-diagonal standard eigensystem problem (pdsptri )
 
     n: size of the matrix
 
     mapA : integer array of size n;
     mapA[i], 0<= i < n, is the processor holding this vector
 
     mapB : integer array of size n;
     mapB[i], 0<= i < n, is the processor holding this vector
    
    mapZ : integer array of size n
    mapZ[i], 0<= i < n, is the processor holding this vector
    
    output:
    
    isize = size of integer buffer
    rsize = size of buffer for DoublePrecision precision
    ptr_size= size of double precision buffer for pointers
    
    at this point one can allocate the required memory using C
    or using off-sets from Fortran arrays
*/
  
  Integer me, naproc, indx;
  Integer           linfo;

  Integer *iptr;
  char            msg[ 25 ];
  Integer nvecsA, nvecsB, nvecsZ;
  
  extern Integer count_list();
  extern Integer mxmynd_(), mxnprc_();
  extern void   l_exit_();
  extern void   memreq_(); 
  extern void     mxinit_();
  
   mxinit_();

   me     = mxmynd_();
   naproc = mxnprc_();

   strcpy( msg,  "Error in fmemreq." );  
  
   linfo = 0;

   if ( type  == NULL )
     linfo = -1;
   else if ( n == NULL )
     linfo = -2;
   else if ( mapA == NULL )
     linfo = -3;
   else if ( mapB == NULL )
     linfo = -4;
   else if ( mapZ == NULL )
     linfo = -5;
   else if ( isize == NULL )
     linfo = -6;
   else if ( rsize == NULL )
     linfo = -7;
   else if ( ptr_size == NULL )
     linfo = -8;
   else if ( iscratch == NULL )
     linfo = -9;

   if ( linfo != 0 ) {
     l_exit_( &linfo, msg );
     return;
   }
  
   if (( *type < 0 ) ||  ( *type > 2 )) {
     linfo = -1;
     l_exit_( &linfo, msg );
     return;
   }
  
   if ( *n < 1 ) {
     linfo = -2;
     l_exit_( &linfo, msg );
     return;
   }
  
  iptr = mapA;
  for ( indx = 0; indx < *n; indx++ ) {
    if ( (iptr++ ) == NULL ) {
      linfo = -3;
      l_exit_(&linfo, msg);
      return;
    }
    if ( mapA[indx] < 0  || mapA[indx] >= naproc ) {
      linfo = -3;
      l_exit_(&linfo, msg);
      return;
    }
  }
  
  iptr = mapB;
  if ( *type == 0 ) {
    for ( indx = 0; indx < *n; indx++ ) {
      if ( (iptr++ ) == NULL ) {
	linfo = -4;
	fprintf(stderr, " me = %d error in mapB in fmemreq.c \n", me);
	l_exit_(&linfo,msg);
	return;
      }
      if ( mapB[indx] < 0  || mapB[indx] >= naproc ) {
	linfo = -4;
	l_exit_(&linfo, msg);
	return;
      }
    }
  }
  
  iptr = mapZ;
  for ( indx = 0; indx < *n; indx++ ) {
    if ( (iptr++ ) == NULL ) {
      linfo = -5;
      fprintf(stderr, " me = %d error in mapZ in fmemreq.c \n",me);
      l_exit_(&linfo,msg);
      return;
    }
    if ( mapZ[indx] < 0  || mapZ[indx] >= naproc ) {
      linfo = -5;
      l_exit_(&linfo, msg);
      return;
    }
  }
  
#ifdef DEBUG0
  fprintf(stderr, " just before memreq_ in memreq_f me = %d \n", me);
#endif
  
  memreq_( type, n, mapA, mapB, mapZ, isize, rsize, ptr_size, iscratch);
  
  nvecsA = 0;
  if( *type == 1  || *type == 0 ) nvecsA = count_list( me, mapA, n);  
  
  nvecsB = 0;
  if( *type == 0  ) nvecsB = count_list( me, mapB, n);  
  
  nvecsZ = count_list( me, mapZ, n);  
  
  /*
   *  fprintf( stderr, " fmemreq ptr_size=%d, nvecsA,B,Z= %d %d %d \n",
   *                   *ptr_size, nvecsA,nvecsB,nvecsZ);
   */
  
  *isize += 0;
  *rsize += 0;
  *ptr_size += (nvecsA + nvecsB + nvecsZ + 3);
  
   return;
}
