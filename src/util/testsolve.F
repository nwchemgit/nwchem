      program test
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
c
      logical status, ga_iter_lsolve, ga_iter_diag
      integer dim, n, g_a, g_b, g_soln, g_resid, i, j, iter, g_v
      integer iter_solve
      common /info/ n, g_a
      integer nroot, heap
      parameter (dim = 10, nroot=3)
      double precision col(dim), acc(nroot), tol, scale, sum
      double precision eval(nroot), used, junk
      double precision evals(dim)
      double precision drandm
      external drandm, ga_iter_lsolve, ga_iter_diag
      external test_precond, test_product
c
      n = dim
      call pbeginf
      heap = max(8, (dim-1)/nnodes()+1)
      heap = dim*(heap + 2)
      heap = heap * 1.1d0 + 128*128*3
      if (.not. ma_init(MT_DBL, heap, heap)) call errquit
     $     ('test: ma_init?', 0, MA_ERR)
      status = ma_set_auto_verify(.true.)
      call ga_initialize()
      write(6,*) ga_nodeid(), ' initialized'
      call util_flush(6)
      call ga_sync()
c

*      call testpcg

c     
c     create the arrays to test the linear equation solver
c
*ga:1:0
*     status = ga_create(MT_DBL, n, n, 'test:v', 0, 0, g_v)
*ga:1:0
      status = ga_create(MT_DBL, n, n, 'test:a', 0, 0, g_a)
      write(6,*) ga_nodeid(), ' created'
      call util_flush(6)
      call ga_sync()

*ga:1:0
*      status = ga_create(MT_DBL, n, 1, 'test:soln', 0, 0, g_soln)
*ga:1:0
*      status = ga_create(MT_DBL, n, 1, 'test:b', 0, 0, g_b)
*ga:1:0
*      status = ga_create(MT_DBL, n, 1, 'test:resid', 0, 0, g_resid)
      junk = drandm(ga_nodeid()*37*2 + 1)
      do i = ga_nodeid()+1, n, ga_nnodes()
         if (i .le. 3) then
            scale = 0.1d0
         else
            scale = 10.0d0
         endif
         col(i) = (1.0d0+0.4d0*(0.5d0 - drandm(0)))*scale
         if (i .gt. 3) scale = 1.0d0
         do j = i+1,n
            col(j) = 0.4d0*(0.5d0 - drandm(0))*scale
         enddo
         call ga_put(g_a, i, n, i, i, col(i), 1)
         call ga_put(g_a, i, i, i, n, col(i), 1)
      enddo
c
      used = -tcgtime()
      call ga_orthog(g_a, -999999, .false.)
      used = (used + tcgtime())
      if (ga_nodeid().eq.0) then
         write(6,*) ' new used ', used
         call util_flush(6)
      endif
*      call ga_dgemm('t','n',n,n,n,1.0d0,g_a, g_a, 0.0d0, g_v)
*      call ga_get_diagonal(g_v, col)
*      sum = ga_ddot(g_v, g_v)
*      if (ga_nodeid().eq.0) then
*         write(6,*) ' diagonal ', sum
*         call util_flush(6)
*      endif
      call ga_sync()
c
c     Add small amount of ranom noise in and call iterative orthog
c
      do i = ga_nodeid()+1, n, ga_nnodes()
         call ga_get(g_a, 1, n, i, i, col(1), 1)
         do j = 1, n
            col(j) = col(j) + 1.0d-8*(drandm(0)-0.5d0)
         enddo
         call ga_put(g_a, 1, n, i, i, col(1), 1)
      enddo
c
      used = -tcgtime()
      call ga_iter_orthog(-1, g_a, .false., .false.)
      used = (used + tcgtime())
      if (ga_nodeid().eq.0) then
         write(6,*) ' old used ', used
         call util_flush(6)
      endif
c      
      goto 10000
c
      if (ga_nodeid() .eq. 0) then
         do i = 1, n
            col(i) = 0.5d0 - drandm(0)
         enddo
         call ga_put(g_b, 1, n, 1, 1, col, 1)
      endif
c
      write(6,*) ' Solver WITHOUT projection'
      tol = 0.0001d0
      call ga_zero(g_soln)
      do iter = 1, 3
         if (ga_iter_lsolve(n, 40, 5, tol, test_precond,
     $        test_product, .true., .true., g_b, g_soln,
     $        acc, iter_solve)) then
            if (ga_nodeid() .eq. 0) write(6,*) ' OK'
         else
            if (ga_nodeid() .eq. 0) write(6,*) ' not OK'
         endif
         tol = tol*0.0001d0
      enddo
c
      write(6,*) ' Solver using projection'
      call ga_init_project(n, 3, 10, 10, 0.05d0, test_precond,
     $     test_product, .true., 0.1d0)
      tol = 0.0001d0
      call ga_zero(g_soln)
      do iter = 1, 3
         if (ga_iter_lsolve(n, 40, 5, tol, test_precond,
     $        test_product, .true., .true., g_b, g_soln, 
     $        acc, iter_solve)) then
            if (ga_nodeid() .eq. 0) write(6,*) ' OK'
         else
            if (ga_nodeid() .eq. 0) write(6,*) ' not OK'
         endif
         tol = tol*0.0001d0
      enddo
c
      call ga_term_project()
c
      call ga_diag_std_seq(g_a, g_v, evals)
      write(6,*) ' Correct evals'
      call output(evals, 1, n, 1, 1, n, 1, 1)
      tol = 0.01d0
c
      status = ga_destroy(g_soln)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 3, 'test:soln', 0, 0, g_soln))
     $     call errquit('gac failed', 0, GA_ERR)
c      
      call ga_zero(g_soln)
      do iter = 1, 3
         if (ga_iter_diag(n, nroot, 40, 10, tol, test_precond,
     $        test_product, .true., 0.0d0,
     $        g_soln, eval, acc, iter_solve)) then
            if (ga_nodeid() .eq. 0) write(6,*) ' OK'
         else
            if (ga_nodeid() .eq. 0) write(6,*) ' not OK'
         endif
         tol = tol*0.01d0
      enddo
c
      status = ga_destroy(g_resid)
      status = ga_destroy(g_b)
      status = ga_destroy(g_soln)
      status = ga_destroy(g_a)
      status = ga_destroy(g_v)
c      
10000 call ga_terminate
      call pend
c
      end
      subroutine test_product(acc, g_x, g_ax)
      implicit none
#include "global.fh"
      double precision acc
      integer g_x, g_ax
      integer n, g_a
      common /info/ n, g_a
c
      call ga_dgemm('n', 'n', n, 1, n, 1.0d0, g_a, g_x, 0.0d0, g_ax)
c
      end
      subroutine test_precond(g_r, shift)
      implicit none
#include "global.fh"
      integer n, g_a
      integer g_r
      common /info/ n, g_a
      double precision x, shift, denominator, numerator
      integer i
c
      call ga_sync()
      do i = ga_nodeid()+1, n, ga_nnodes()
         call ga_get(g_a, i, i, i, i, denominator, 1)
         call ga_get(g_r, i, i, 1, 1, numerator, 1)
         denominator = denominator - shift
         if (abs(denominator) .lt. 0.1d0) denominator =
     $        sign(0.1d0, denominator)
         x = numerator / denominator
         call ga_put(g_r, i, i, 1, 1, x, 1)
      enddo
c
      end
      double precision function drandm(iflag)
      implicit none
      integer iflag
      double precision util_random
      external util_random
      drandm = util_random(iflag)
      end
      subroutine testpcg
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer n, ls_max
      parameter (n = 10, ls_max=20)
      integer iwork(10)
      double precision dwork(3*ls_max+5)
      integer g_grad, g_work, g_s, g_x, g_a
      double precision step
      double precision ls_tol
      parameter (ls_tol = 1.0d-2)
      double precision eacc
      parameter (eacc = 1.0d-10)
      character*16 mode
      integer i, iter, j
      double precision x, e, local(n)
      double precision drandm
      external drandm
      logical print, conjugacy, status
      data print /.true./, conjugacy /.true./
c
*ga:1:0
      if (.not. ga_create(MT_DBL, n, n, 'testpcg: matrix', 0, 0,
     $     g_a)) call errquit('testpcg: a', n, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: grad', 0, 0,
     $     g_grad)) call errquit('testpcg: grad', n, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: work', 0, 0,
     $     g_work)) call errquit('testpcg: work', n, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: step', 0, 0,
     $     g_s)) call errquit('testpcg: step', n, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: x', 0, 0,
     $     g_x)) call errquit('testpcg: step', n, GA_ERR)
c
c     Initial guess for solution x
c
      do i = ga_nodeid()+1, n, ga_nnodes()
         x = 0.1*drandm(0)
         call ga_put(g_x, i, i, 1, 1, x, 1)
      enddo
      call ga_sync
c
c     Fill in the matrix with something positive definite
c
      do i = ga_nodeid()+1, n, ga_nnodes()
         x = 0.0d0
         do j = 1,n
            local(j) = 0.1*drandm(0)
            x = x + local(j)
         enddo
         local(i) = 1.0d0 + x
         call ga_put(g_a, 1, n, i, i, local, n)
      enddo
      call ga_symmetrize(g_a)
      do i = 1, n
         call ga_get(g_a, i, i, i, i, local(i), 1)
      enddo
c      
      mode = 'startup'
 10   continue
      call ga_pcg_minim(n, iter, e, g_grad, g_work, g_s, step, ls_tol,
     $     ls_max, eacc, conjugacy, print, iwork, dwork, mode)
      if (mode .eq. 'energy+gradient') then
         call ga_sync
         call ga_dadd(step, g_s, 1.0d0, g_x, g_x)
         call ga_dgemm('n', 'n', n, 1, n, 1.0d0, g_a, g_x, 0.0d0,g_grad)
         e = 0.5d0*ga_ddot(g_x, g_grad)
         call ga_dadd(-step, g_s, 1.0d0, g_x, g_x)
         goto 10
      else if (mode .eq. 'precondition') then
         call ga_sync
         do i = ga_nodeid()+1, n, ga_nnodes()
            call ga_get(g_grad, i, i, 1, 1, x, 1)
            x = x/(local(i) + 2.0d0)
            call ga_put(g_work, i, i, 1, 1, x, 1)
         enddo
         call ga_sync
         goto 10
      else if (mode .eq. 'accept step') then
         call ga_maxelt(g_grad, x)
         call ga_dadd(step, g_s, 1.0d0, g_x, g_x)
         if (x.gt.0.00001d0) goto 10
      else
         write(6,*) mode
         call errquit('unknown mode', 0, GA_ERR)
      endif
c
      status = ga_destroy(g_a)
      status = ga_destroy(g_x)
      status = ga_destroy(g_s)
      status = ga_destroy(g_grad)
      status = ga_destroy(g_work)
c
      end
