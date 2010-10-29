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

#define max(a,b) ((a) > (b) ? (a) : (b))

void mxm25 ( n1, n2, rowQ, mapQ, m, colW, mapW, colZ, iwork, work)
     Integer *n1, *n2, mapQ[], *m, mapW[], iwork[];
     DoublePrecision **rowQ, **colW, **colZ, work[];

/**************************************************************
 *
 * Subroutine mxm25
 *
   This subroutine does the matrix multiplication Z <- Q*W
       
   where matrix Q is a n1 x n2 general matrix in packed storage format,
         distributed by rows,

         matrix W is a general n2 x m matrix in packed storage format,
         distributed by columns.

         and matrix Z is the n1 x m product Q*W, distributed the same as W.
       
   In this version we assume that Q and W share the same processors
   ( perhaps different data distributions )

   It is ok to have colZ = colW.  However, in this
   case each colZ[i]=colW[i] must point to a vector of
   length = MAX{ n1, n2 }.

   ARGUMENTS
   ---------
   In the following:

     INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"

     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsW = number of entries in mapW equal to me
                  (= count_list( me, mapW, m ))
     nvecsQ = number of entries in mapQ equal to me
                  (= count_list( me, mapQ, n1 ))
     nvecsQ_max = maximum number of entries in mapQ equal to i,
                  i = 0 to nprocs-1
                  (= max over i of count_list( i, mapQ, n1 ))
     sDP    = sizeof( DoublePrecision )

   n1 ..... (input) INTEGER
            The number of rows in matrix Q

   n2 ..... (input) INTEGER
            The number of columns in matrix Q
            and the number of rows in matrix W

   rowQ ... (input) array of pointers to DoublePrecision,
                    length (nvecsQ)
            The part of matrix Q owned by this processer stored
            in packed format, i.e., rowQ[i] points to the start
            of the i-th row of Q
            owned by this processor, i = 0 to nvecsQ-1.
                
   mapQ ... (input) INTEGER array, length (n1)
            The i-th row of Q is owned by processor
            mapQ[i], i = 0 to n1-1.

   m ...... (input) INTEGER
            The number of columns in W.
       
   colW ... (input) array of pointers to DoublePrecision,
                           length (nvecsW)

            The part of matrix W owned by this processer stored
            in packed format by columns, i.e., colW[i] points
            to the first element of the i-th column of W
            owned by this processor, i = 0 to nvecsW-1.
                
   mapW ... (input) INTEGER array, length (m)
            The i-th column of W is 
            owned by processor mapW[i], i = 0 to m-1.
                
   colZ ... (output) array of pointers to DoublePrecision,
                     length (nvecsW)

             The result matrix Z = Q * W stored identical to W.

             It is ok to have colZ = colW.  However, in this
             case each colZ[i]=colW[i] must point to a vector of
             length = MAX{ n1, n2 }.

   iwork .. (workspace) INTEGER array, length ( 3 * n1 + 2 * m )

   work ... (workspace) DOUBLE PRECISION array, 
                        length ( maximum ( n1 + m, (nvecsW + 2 * nvecsQ_max) * n2 )
*/
{
  static Integer ONE = 1;
  Integer jj, i, nvecsQ, nvecsW, l1, l2;
  Integer *mapvecQ, *mapvecW;
  Integer isize, indx, jndx, me;
  Integer *proclist, nprocs, *mapvec_in, last_proc, next_proc;
  Integer *iscrat, me_indx;
  Integer maxsz, rsize, osize;
  Integer *ijunk, idummy;
  
  DoublePrecision *ptr;
  DoublePrecision *buffer, *in_buffer, *out_buffer;
  DoublePrecision *dd_ptr, *dd_ptr2;
  
  extern Integer count_list();
  extern Integer fil_mapvec_();
  
  extern Integer mxwrit_(), mxread_();
  extern Integer menode_(), mxmynd_();
  
  extern Integer fil_mapvec_();
  extern Integer count_list();
  extern Integer indxL ();
  extern Integer mxwrit_(), mxread_();
  extern Integer mxmynd_();
  
  /*
    blas call
    */

  extern DoublePrecision ddot_();
  extern Integer reduce_list2();
  extern Integer peigs_cmod_();
  
  extern void mxsync_();
  extern void gshellsort_();
  extern void fil_dbl_lst();
  extern void daxpy_();
  extern void dcopy_();
  
  l1 = *n1;
  l2 = *n2;
  
  me = mxmynd_();
  iscrat = iwork;
  mapvecQ = iscrat;
  nvecsQ = fil_mapvec_( &me, &l1, mapQ, mapvecQ );
  
  iscrat += nvecsQ;
  mapvecW = iscrat;
  nvecsW = fil_mapvec_( &me, m, mapW, mapvecW );
  
  if ( nvecsW + nvecsQ == 0 )
    return;
  
  iscrat += nvecsW;

  proclist = iscrat;


  ijunk = (Integer *) work;
  
  for ( indx = 0; indx < l1; indx++ )
    ijunk[indx] = mapQ[indx];
  
  jndx = l1;
  for ( indx = 0; indx <  *m; indx++ ) {
    ijunk[jndx] = mapW[indx];
    jndx++;
  }
  
  /*
    the list of processors for this computation
    */
  
  indx = l1 + *m ;
  
  nprocs = reduce_list2( indx, &ijunk[0], proclist);
  gshellsort_( &nprocs, proclist);
  
  /*
    the number of distinct processors
    */
  
  iscrat += nprocs;
  mapvec_in = iscrat;
  indx = indxL ( me, nprocs, proclist);
  me_indx = indx;
  
  /*
    last_proc = ((indx - 1) % nprocs);
    */
  
  idummy = indx - 1 ;
  idummy += nprocs;
  last_proc = peigs_cmod_( &idummy, &nprocs);
  last_proc = proclist[last_proc];
  
  /*
    next_proc = (indx + 1) % nprocs;
    */
  
  idummy = indx + 1;
  idummy += nprocs;
  next_proc = peigs_cmod_( &idummy, &nprocs);
  next_proc = proclist[next_proc];
  
  maxsz = 0;
  for ( indx = 0; indx < nprocs; indx++ ) {
    isize = count_list( proclist[indx], mapQ,& l1);   /* the number Q row vectors that I own */
    maxsz = max( isize, maxsz );
  }
  
  
  /*
    
    load what I own of W (columns) into a buffer;
    for multiplication purpose
    just in case the W pointers and the Z pointers are the same
    
    */
  
  buffer = (DoublePrecision *) work;
  
  indx = 0;
  for ( jndx = 0; jndx < nvecsW; jndx++ ){
    ptr = colW[jndx];
    dcopy_( &l2, ptr, &ONE, buffer + indx, &ONE);
    indx += l2;
  }
  buffer += indx;
  
  /*
    zero out the matrix Z
    */
  
  for ( jndx = 0; jndx < nvecsW; jndx++ ){
    ptr = colZ[jndx];
    fil_dbl_lst ( l1, ptr, 0.0e0);
  }
  
  /*
    copy rows of Q into output buffer for shifting
    shift matrix Q around and multiply
    */

  indx = 0;
  for ( jndx = 0; jndx < nvecsQ; jndx++ ){
    dcopy_( &l2, rowQ[jndx], &ONE, buffer + indx, &ONE);
    indx+= l2;
  }
  
  for ( indx = 0; indx < nvecsQ; indx++ )
    mapvec_in[indx] = mapvecQ[indx];
  
  /*
    just computing matrix multiplication within our
    own vectors
    */
  
  for ( indx = 0; indx < nvecsW; indx++ ){  
    ptr = colZ[indx];
    for (jndx = 0; jndx < nvecsQ; jndx++ ){
      jj = mapvec_in[jndx];
      ptr[jj] = ddot_( &l2, &buffer[jndx*l2], &ONE, &work[indx*l2], &ONE);
    }
  }
  
  /*
    this is the block that gets send around; first multiply what I have
    */
  
  in_buffer = buffer;
  out_buffer = in_buffer + maxsz*l2;
  
  isize = nvecsQ*l2;
  osize = isize;
  
  for ( i = 0; i < nprocs-1 ; i++ ) {
    if( osize > 0 ) 
       dcopy_(&osize, in_buffer, &ONE, out_buffer, &ONE);
    
    idummy = 2;
    if ( peigs_cmod_( &me_indx, &idummy) ==  0 ) {
      
      rsize = osize*sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxwrit_( out_buffer, &rsize, &next_proc, &i );

      /*
        indx = (( me_indx - i -1) % nprocs);
        */
      
      idummy = me_indx - i -1;
      idummy += nprocs;
#ifdef DEBUG1
      fprintf(stderr, " me = %d i = %d idummy %d \n", me, i, idummy);
#endif
      indx = peigs_cmod_( &idummy, &nprocs);
      indx = proclist[indx];
      
      /*
	the current processor id number
	*/


      isize = fil_mapvec_( &indx, &l1, mapQ, mapvec_in);
      nvecsQ = isize;
      isize *= l2;
      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );

    }
    else {
      /*
        indx = (( me_indx - i -1 ) % nprocs);
        */
      
      idummy = me_indx - i - 1;
      idummy += nprocs;
      
#ifdef DEBUG1
      fprintf(stderr, "222 me = %d i = %d idummy %d \n", me, i, idummy);
#endif
      
      indx = peigs_cmod_(&idummy, &nprocs);
#ifdef DEBUG1
      fprintf(stderr, " me = %d indx = %d idummy %d \n", me, indx, idummy);
#endif

      indx = proclist[indx]; /* the current processor id number */
      isize = fil_mapvec_( &indx, &l1, mapQ, mapvec_in);
      nvecsQ = isize;
      isize *= l2;
      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
      
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxwrit_( out_buffer, &rsize, &next_proc, &i );
    }
    
    osize = isize;
    for (jndx = 0; jndx < nvecsQ; jndx++ ){
      jj = mapvec_in[jndx];
      dd_ptr = &buffer[jndx*l2];
      for ( indx = 0; indx < nvecsW; indx++ ){  
	dd_ptr2 = &work[indx*l2];
	ptr = colZ[indx];
	ptr[jj] = ddot_( &l2, dd_ptr, &ONE, dd_ptr2, &ONE);
      }
    }
  }

  /*
    c    free(proclist);
  */
  return;
}
