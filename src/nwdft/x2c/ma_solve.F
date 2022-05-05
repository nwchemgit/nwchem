      subroutine ma_solve(g_a, g_bx)

      implicit none
      
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"      

c     purpose: solve linear equation system AX = B with square but
c     perhaps non-symmetric matrix. ga_solve relies on SCALAPACK for
c     reliable results, which we may not have available in a given build


      integer g_a   ! input: square coefficient matrix A
      integer g_bx  ! in-out: rhs B and solution matrix X

c     local variables

      integer g_x, g_y

      integer l_amat, k_amat
      integer l_bmat, k_bmat
      integer l_xmat, k_xmat
      integer l_ipiv, k_ipiv
      integer l_tmp, k_tmp

      integer na1, na2, nb1, nb2, ntmp, info1, info2, dbg, type
      integer n, nrhs, n2

      logical master

      character*8 pname

c     ================================================================

      dbg = 0 ! debug setting during code development
      
      pname = "ma_solve"

      master = ga_nodeid().eq.0

      if (dbg>0.and.master) write(luout,*) 'entering '//pname

      call ga_inquire(g_a,type,na1,na2)
      call ga_inquire(g_bx,type,nb1,nb2)

      if (na2.ne.nb1) then
        call errquit(
     &    pname//': input ga dimensions do not match',nb1,GA_ERR)
      end if

      if (na2.ne.na1) then
        call errquit(
     &    pname//': sorry I expect matrix A to be square',na2,GA_ERR)
      end if

c$$$      if (nb2.gt.na1) then
c$$$        call errquit(
c$$$     &    pname//': too many right hand sides',nb2,GA_ERR)
c$$$      end if

      n = na1
      nrhs = nb2

      if (dbg>0.and.master) write(luout,*) pname//': n,nrhs = ',n, nrhs

      n2 = n*n

      info1 = 0
      info2 = 0

c     allocate MA arrays and retrieve data from GAs

c     amat
      if (.not. ma_push_get(mt_dbl, n2, 'amat', l_amat, k_amat))
     &  call errquit(pname//': failed allocation amat', n2, MA_ERR)

c     bmat
      if (.not. ma_push_get(mt_dbl, n*nrhs, 'bmat', l_bmat, k_bmat))
     &  call errquit(pname//': failed allocation bmat', n*nrhs, MA_ERR)      

c     tmp
      if (.not. ma_push_get(mt_dbl, n, 'tmp', l_tmp, k_tmp))
     &  call errquit(pname//': failed allocation tmp', n, MA_ERR)

c     ipiv
      if (.not. ma_push_get(mt_int, n, 'ipiv', l_ipiv, k_ipiv))
     &  call errquit(pname//': failed allocation ipiv', n, MA_ERR)      

      call dfill(n2,     0.0d0, dbl_mb(k_amat),  1)
      call dfill(n*nrhs, 0.0d0, dbl_mb(k_bmat),  1)
      call dfill(n,      0.0d0, dbl_mb(k_tmp),   1)

      call ga_get(g_a,  1, n, 1, n,    dbl_mb(k_amat), n)
      call ga_get(g_bx, 1, n, 1, nrhs, dbl_mb(k_bmat), n)

      call dgetrf(n, n, dbl_mb(k_amat), n, int_mb(k_ipiv), info1)
      if (info1.ne.0)  call errquit(
     &  pname//': dgetrf info not zero', info1 ,UNKNOWN_ERR)      

c     OLD CODE:
c$$$c     invert matrix A and store in temp GA
c$$$      
c$$$      call dgetri(n,    dbl_mb(k_amat), n, int_mb(k_ipiv),
c$$$     &  dbl_mb(k_tmp), n, info2)
c$$$      if (info2.ne.0)  call errquit(
c$$$  &  pname//': dgetri info not zero', inf2 ,UNKNOWN_ERR)

c     new code using drgets instead of calculatinf a matrix inverse
c     explicitly

      call dgetrs (
     &  'n',
     &  n, 
     &  nrhs,
     &  dbl_mb(k_amat),
     &  n,
     &  int_mb(k_ipiv), 
     &  dbl_mb(k_bmat),
     &  n,
     &  info2) 	      

       if (info2.ne.0)  call errquit(
     &  pname//': dgetrs info not zero', info2 ,UNKNOWN_ERR)

       if(dbg>0.and.master)
     &   write(luout,*) pname//': info1, info2 = ',info1, info2

c      MA array B now has the solution. Put back into g_bx
       call ga_put(g_bx, 1,n, 1, nrhs, dbl_mb(k_bmat), n)

c$$$      OLD CODE here and commented below:       
c$$$c     allocate temp GA and store inverse of A      
c$$$      if (.not.ga_duplicate(g_a, g_x, 'g_x')) call 
c$$$     &  errquit(pname//': ga_dup failed g_x ',101, GA_ERR)
c$$$
c$$$      call ga_put(g_x, 1,n, 1,n, dbl_mb(k_amat), n)

      
c     MA arrays no longer needed

      if (.not.ma_chop_stack(l_amat)) call
     &  errquit(pname//': ma_chop_stack failed k_amat',l_tmp,MA_ERR)

c$$$c     allocate one more GA, perform X = A**(-1) B, and
c$$$c     store the result in g_bx
c$$$      
c$$$      if (.not.ga_duplicate(g_a, g_y, 'g_y')) call 
c$$$     &  errquit(pname//': ga_dup failed g_y ',102, GA_ERR)
c$$$
c$$$      call ga_dgemm('N', 'N', n, n, n, 1.d0,
c$$$     &  g_x, g_bx, 0.0d0, g_y)
c$$$
c$$$      call ga_copy(g_y, g_bx)
c$$$
c$$$c     deallocate GAs
c$$$
c$$$      if (.not. ga_destroy(g_y))
c$$$     &  call errquit(pname//': error destroy g_y',102, GA_ERR)
c$$$
c$$$      if (.not. ga_destroy(g_x))
c$$$     &  call errquit(pname//': error destroy g_x',101, GA_ERR) 

c     all done. return

      if (dbg>0.and.master) write(luout,*) 'leaving '//pname

      return

      end
