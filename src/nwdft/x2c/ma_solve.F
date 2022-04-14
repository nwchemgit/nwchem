      subroutine ma_solve(g_a, g_bx)
      
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"      

c     purpose: solve linear equation system AX = B with square
c     non-symmetric matrix. ga_solve relies on SCALAPACK for reliable
c     results, which we may not have available in a given build

c     routine uses 2n**2 words of extra floating point storage, with
c     input matrices of dimension n x n

      integer g_a   ! inout: coefficient matrix A
      integer g_bx ! in: rhs. B. out: solution matrix X

c     local variables

      integer g_x, g_y

      integer l_amat, k_amat
      integer l_bmat, k_bmat
      integer l_xmat, k_xmat
      integer l_ipiv, k_ipiv
      integer l_tmp, k_tmp

      integer n, n2, ntmp, inf1, inf2, dbg, type

      logical master

      character*8 pname

c     ================================================================

      dbg = 0 ! debug setting during code development
      
      pname = "ma_solve"

      master = ga_nodeid().eq.0

      if (dbg>0.and.master) write(luout,*) 'entering '//pname

      n = 0

      call ga_inquire(g_a,type,n,n)
      call ga_inquire(g_bx,type,ntmp,ntmp)

      if (ntmp.ne.n) then
        call errquit(
     &    pname//': input ga dimensions do not match',n,GA_ERR)
      end if

      if (dbg>0.and.master) write(luout,*) pname//'n = ',n

      n2 = n*n

      info1 = 0
      info2 = 0

c     allocate MA arrays and retrieve data from GAs

c     amat
      if (.not. ma_push_get(mt_dbl, n2, 'amat', l_amat, k_amat))
     &  call errquit(pname//': failed allocation amat', n2, MA_ERR)

c     tmp
      if (.not. ma_push_get(mt_dbl, n, 'tmp', l_tmp, k_tmp))
     &  call errquit(pname//': failed allocation tmp', n, MA_ERR)

c     ipiv
      if (.not. ma_push_get(mt_int, n, 'ipiv', l_ipiv, k_ipiv))
     &  call errquit(pname//': failed allocation ipiv', n, MA_ERR)      

      call dfill(n2, 0.0d0, dbl_mb(k_amat),  1)
      call dfill(n,  0.0d0, dbl_mb(k_tmp),   1)

      call ga_get(g_a, 1, n, 1, n, dbl_mb(k_amat), n)

c     invert matrix A and store in temp GA

      call dgetrf(n, n, dbl_mb(k_amat), n, int_mb(k_ipiv), inf1)
      call dgetri(n,    dbl_mb(k_amat), n, int_mb(k_ipiv),
     &  dbl_mb(k_tmp), n, inf2)

      if(dbg>0.and.master)
     &  write(luout,*) pname//': inf1, inf2 = ',inf1, inf2

c     allocate temp GA and store inverse of A
      
      if (.not.ga_duplicate(g_a, g_x, 'g_x')) call 
     &  errquit(pname//': ga_dup failed g_x ',101, GA_ERR)

      call ga_put(g_x, 1,n, 1,n, dbl_mb(k_amat), n)
      
c     MA arrays no longer needed

      if (.not.ma_chop_stack(l_amat)) call
     &  errquit(pname//': ma_chop_stack failed k_amat',l_tmp,MA_ERR)

c     allocate one more GA, perform X = A**(-1) B, and
c     store the result in g_bx
      
      if (.not.ga_duplicate(g_a, g_y, 'g_y')) call 
     &  errquit(pname//': ga_dup failed g_y ',102, GA_ERR)

      call ga_dgemm('N', 'N', n, n, n, 1.d0,
     &  g_x, g_bx, 0.0d0, g_y)

      call ga_copy(g_y, g_bx)

c     deallocate GAs

      if (.not. ga_destroy(g_y))
     &  call errquit(pname//': error destroy g_y',102, GA_ERR)

      if (.not. ga_destroy(g_x))
     &  call errquit(pname//': error destroy g_x',101, GA_ERR) 

c     all done. return

      if (dbg>0.and.master) write(luout,*) 'leaving '//pname

      return

      end
