      subroutine ga_chol_seq(g_a, g_s,info)
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
      integer g_a ! input matrix
      integer g_s ! output matrix
      integer l_a_in, i_a_in ! ma handle and pointer
      integer ga_chol_msg
c     
c     
c     This subroutine performs a submatrix Choleski factorization of a 
c     symmetric positive matrix on one processor (debug alternative
c     to PeIGS)
c     
c     the input matrix is not destroyed
c     
c     n    - problem size
c     proc - # of processors
c     
c     
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     
      integer iproc, pan, istart, iend, k, n, npan, proc, myproc
      integer ik, info
      logical status
c     
c     
      logical oactive ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimS1, dimS2, typeS
c======================================================================
c     
c***  check environment
c     
      myproc = ga_nodeid()
      proc   = ga_nnodes()
c     
c***  check GA info
c     
      call ga_check_handle(g_a, 'ga_chol: a')
      call ga_check_handle(g_s, 'ga_chol: s')
      call ga_inquire(g_a,  typeA, dimA1, dimA2)
      call ga_inquire(g_s,  typeS, dimS1, dimS2)
c     
      if (dimA1 .ne. dimA2) then
        call ga_error('ga_chol: can only fact square matrix ', 1)
      else if (dimS1 .ne. dimS2) then
        call ga_error('ga_chol: can only fact square matrix ', 3)
      else if(dimA1.ne.dimS1 )then
        call ga_error('ga_chol: A & S dimension mismatch ', 4)
      endif
c     
      n = dimA1

      call ga_sync()
      proc=ga_nnodes()

      oactive = myproc .eq. 0

      if (oactive) then
c****
c****   allocate local storage
c****
        status = ma_push_get(mt_dbl,n*n,'a_in',l_a_in,i_a_in)        
c****
c****   collect matrix from global array
c****
        call ga_get(g_a, 1, n, 1, n, Dbl_MB(i_a_in), n )
c****
c****   cholesky factorize
c****
c****        write (6,*) "ga_chol_seq: incoming matrix"
c****        call matout(Dbl_MB(i_a_in),n,n,n)
        call dpotrf( 'LowerTriangle', n, Dbl_MB(i_a_in), n, info )
        if ( info .ne. 0) then
          write (6,*)
     &       'WARNING:  ga_chol_seq: LAPACK routine spotrf returns ',
     &       info
        endif 
c****
c****   store back to global array
c****
        call ga_put(g_s, 1, n, 1, n, Dbl_MB(i_a_in), n )
c****
c****   deallocate local storage
c****
        status = ma_pop_stack(l_a_in)
      endif 
      ga_chol_msg = 32760 - 17 
      call ga_brdcst(ga_chol_msg, info, mitob(1), 0)
      call ga_sync()
      return 
c     
      end


