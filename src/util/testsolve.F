      program test
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
c
      logical status, ga_iter_lsolve, ga_iter_diag
      integer dim, n, g_a, g_b, g_soln, g_resid, i, j, iter, g_v
      integer iter_solve
      common /info/ n, g_a
      integer nroot
      parameter (dim = 100, nroot=3)
      double precision col(dim), acc(nroot), tol, scale
      double precision eval(nroot)
      double precision evals(dim)
      double precision drandm
      external drandm, ga_iter_lsolve, ga_iter_diag
      external test_precond, test_product
c
      n = dim
      call pbeginf
      if (.not. ma_init(MT_DBL, 100000, 100000)) call errquit
     $     ('test: ma_init?', 0)
      status = ma_set_auto_verify(.true.)
      call ga_initialize()
c

*      call testpcg

c     
c     create the arrays to test the linear equation solver
c
      status = ga_create(MT_DBL, n, n, 'test:v', 1, 1, g_v)
      status = ga_create(MT_DBL, n, n, 'test:a', 1, 1, g_a)
      status = ga_create(MT_DBL, n, 1, 'test:soln', 1, 1, g_soln)
      status = ga_create(MT_DBL, n, 1, 'test:b', 1, 1, g_b)
      status = ga_create(MT_DBL, n, 1, 'test:resid', 1, 1, g_resid)
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
      write(6,*) ' before '
      call ga_print(g_a)
      call ga_orthog(g_a, -999999, .false.)
      write(6,*) ' after '
      call ga_print(g_a)
      call ga_dgemm('t','n',n,n,n,1.0d0,g_a, g_a, 0.0d0, g_v)
      write(6,*) ' identity'
      call ga_print(g_v)
      call errquit('done')
c
      if (ga_nodeid() .eq. 0) then
         do i = 1, n
            col(i) = 0.5 - drandm(0)
         enddo
         call ga_put(g_b, 1, n, 1, 1, col, 1)
      endif
c
      write(6,*) ' Solver WITHOUT projection'
      tol = 0.0001d0
      call ga_zero(g_soln)
      do iter = 1, 3
         if (ga_iter_lsolve(n, 40, 5, tol, test_precond,
     $        test_product, .true., g_b, g_soln, acc, iter_solve)) then
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
     $        test_product, .true., g_b, g_soln, acc, iter_solve)) then
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
      if (.not. ga_create(MT_DBL, n, 3, 'test:soln', 1, 1, g_soln))
     $     call errquit('gac failed', 0)
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
      call ga_terminate
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
         if (abs(denominator) .lt. 0.1) denominator =
     $        sign(0.1, denominator)
         x = numerator / denominator
         call ga_put(g_r, i, i, 1, 1, x, 1)
      enddo
c
      end
      double precision function drandm(iflag)
      implicit none
      integer iflag
#ifdef KSR
      integer irandm
      external irandm
      double precision scale
      parameter(scale=1.0d0/2147483647.0d0)
      drandm = scale*irandm(iflag)
#else
      double precision drand
      external drand
      drandm = drand(iflag)
#endif
      end
      subroutine testpcg
      implicit none
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
      if (.not. ga_create(MT_DBL, n, n, 'testpcg: matrix', 1, 1,
     $     g_a)) call errquit('testpcg: a', n)
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: grad', 1, 1,
     $     g_grad)) call errquit('testpcg: grad', n)
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: work', 1, 1,
     $     g_work)) call errquit('testpcg: work', n)
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: step', 1, 1,
     $     g_s)) call errquit('testpcg: step', n)
      if (.not. ga_create(MT_DBL, n, 1, 'testpcg: x', 1, 1,
     $     g_x)) call errquit('testpcg: step', n)
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
         local(i) = 1.0 + x
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
         if (x.gt.0.00001) goto 10
      else
         write(6,*) mode
         call errquit('unknown mode', 0)
      endif
c
      status = ga_destroy(g_a)
      status = ga_destroy(g_x)
      status = ga_destroy(g_s)
      status = ga_destroy(g_grad)
      status = ga_destroy(g_work)
c
      end
