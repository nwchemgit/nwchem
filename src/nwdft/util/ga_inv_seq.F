      subroutine ga_inv_seq(g_a, g_s)
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
      integer g_a               ! input matrix
      integer g_s               ! output matrix
c
c     construct matrix inverse A^-1 from the lower choleski decomposition
c     L such that A = L L^{\dag}
c
c     NOTE this differs from ga_inverse; which produces the inverse
c     OF L!!!    The functionality here is chosen to reproduce LAPACK
c
c     the input matrix is not destroyed
c
c     proc - # of processors
c     n    - problem size
c
c     
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      integer iproc, k, n,proc, myproc
      integer ik, info
      logical status
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimS1, dimS2, typeS
      integer i_a_in, l_a_in
c
c***  check environment
c
      myproc = ga_nodeid()
      proc   = ga_nnodes()
c
c***  check GA info
c
      call ga_check_handle(g_a, 'ga_inverse: a')
      call ga_check_handle(g_s, 'ga_inverse: s')
      call ga_inquire(g_a,  typeA, dimA1, dimA2)
      call ga_inquire(g_s,  typeS, dimS1, dimS2)
c
c
      if (dimA1 .ne. dimA2) then
         call ga_error('ga_inverse: can only invert square matrix ', 1)
      else if (dimS1 .ne. dimS2) then
         call ga_error('ga_inverse: can only invert square matrix ', 3)
      else if(dimA1.ne.dimS1 )then
         call ga_error('ga_inverse: A & S dimension mismatch ', 4)
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
c****   construct inverse of original matrix from its lower 
c****   Cholesky decomposition L
c****
c****        write (6,*) "ga_inv_seq: incoming matrix"
c****        call matout(Dbl_MB(i_a_in),n,n,n)

        call dpotri( 'LowerTriangle', n, Dbl_MB(i_a_in), n, info )
        if ( info .ne. 0) then
          write (6,*)
     &       'WARNING:  ga_inv_seq: LAPACK routine spotri returns ',
     &       info
        endif 
c****        write (6,*) "ga_inv_seq: output matrix"
c****        call matout(Dbl_MB(i_a_in),n,n,n)
c****
c****   store back to global array
c****
        call ga_put(g_s, 1, n, 1, n, Dbl_MB(i_a_in), n )
c****
c****   deallocate local storage
c****
        status = ma_pop_stack(l_a_in)
      endif 
      call ga_sync()
      return 
c
      end


