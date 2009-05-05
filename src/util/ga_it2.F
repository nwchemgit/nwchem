c $Id: ga_it2.F,v 1.15 2008-08-18 03:08:45 jhammond Exp $
      subroutine test_product(acc,g_x, g_Ax)
      implicit none
      double precision acc
      integer g_x, g_Ax
c
      integer g_a, g_b
      common /testme/ g_A, g_B
      integer n, nvec, type
c
      call ga_inquire(g_x, type, n, nvec)
      call ga_zero(g_Ax)
      call ga_dgemm('n', 'n', n, nvec, n, 1.0d0, g_A, g_x, 0.0d0, g_Ax)
c
      end
      subroutine test_nlproduct(acc,g_x, g_Ax)
      implicit none
      double precision acc
      integer g_x, g_Ax
c
      integer g_a, g_b
      common /testme/ g_A, g_B
      integer n, nvec, type
c
      call ga_inquire(g_x, type, n, nvec)
      call ga_zero(g_Ax)
      call ga_dgemm('n', 'n', n, nvec, n, 1.0d0, g_A, g_x, 0.0d0, g_Ax)
      call ga_add(1.0d0, g_AX, -1.0d0, g_B, g_AX)
c
      end
      subroutine test_precond(g_x,shift)
      implicit none
#include "global.fh"
      integer g_x
      double precision shift
      integer n, nvec, type
      integer i, ivec
      double precision x
c
      call ga_inquire(g_x, type, n, nvec)
      if (ga_nodeid() .eq. 0) then
         do ivec = 1, nvec
            do i = 1, n
               call ga_get(g_x, i, i, ivec, ivec, x, 1)
               x = x / (dble(i) - shift)
               call ga_put(g_x, i, i, ivec, ivec, x, 1)
            end do
         end do
      end if
      call ga_sync()
c
      end
      subroutine ga_lkain_test()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer n, nvec, maxsub, maxiter
      parameter (n=300, nvec=1)
      integer g_x
      external test_precond, test_product, test_nlproduct
c
      integer g_a, g_b
      common /testme/ g_a, g_b
c
      integer i
****      integer info
****      double precision a(n,n), b(n), w(n)
c
      maxsub = 4*nvec
      maxiter = 100
c
      if (.not. ga_create(MT_DBL, n, nvec, 'testx', 0, 0, g_x))
     $     call errquit('test kain', 1, GA_ERR)
      if (.not. ga_create(MT_DBL, n, nvec, 'testx', 0, 0, g_b))
     $     call errquit('test kain', 2, GA_ERR)
      if (.not. ga_create(MT_DBL, n, n, 'testx', 0, 0, g_A))
     $     call errquit('test kain', 3, GA_ERR)
c
      call ga_ran_fill(g_A, 1, n, 1, n)
      call ga_ran_fill(g_b, 1, n, 1, nvec)
      if (ga_nodeid() .eq. 0) then
         do i = 1, n
            call ga_put(g_a, i, i, i, i, 0.5*dble(i), 1)
         end do
      end if
      call ga_sync()
c
      call ga_copy(g_b, g_x)
      call test_precond(g_x,0.0d0)
c
****      call ga_get(g_a, 1, n, 1, n, a, n)
****      call ga_get(g_b, 1, n, 1, nvec, b, n)
****      call dgesv(n, nvec, a, n, w, b, n, info)
****      write(6,*) ' info ', info
****      call ga_put(g_x, 1, n, 1, nvec, b, n)
c
C     This should have something other than zero
      call ga_lkain(0,g_x, g_b,test_product,test_precond,1d-6,maxsub,
     $     maxiter,.true.,.true.)
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
c
      write(6,*)
      write(6,*)
      write(6,*)
      write(6,*) ' DOING NON LINEAR '
      write(6,*)
      write(6,*)
      write(6,*)
      call ga_copy(g_b, g_x)
      call test_precond(g_x,0.0d0)
      maxsub = 10
c
      call ga_kain(g_x, 
     $     test_nlproduct, test_precond, 
     $     1d-6, 
     $     10.0d0, 10.0d0,
     $     maxsub, maxiter, 
     $     .true.)
c
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
      call errquit('done',0, MEM_ERR)
c
      end
      subroutine ga_lkain(rtdb,g_x, g_b, product, precond, 
     $     tol, mmaxsub, maxiter, odiff, oprint)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
c
      integer g_x               ! [input/output] Initial guess/solution
      integer g_b               ! [input] Right-hand side vectors
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag
      integer rtdb
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . optional shift (not used here but used by the diagonalizer)
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.
c
      integer iter, n, nvec, nsub, isub, type, maxsub
      integer g_y, g_Ay, g_Ax, g_r, g_a, g_bb, g_c, g_xold, g_Axold
      double precision rmax, acc, ga_svd_tol
      logical converged
      logical odebug
c
      odebug = util_print('debug lsolve', print_never) .and. 
     $     ga_nodeid().eq.0
      if (.not.rtdb_get(rtdb,'cphf:acc',mt_dbl,1,acc)) acc=0.01d0*tol
      call ga_inquire(g_x, type, n, nvec)
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec
c
      if (.not.rtdb_get(rtdb,'cphf:ga_svd_tol',mt_dbl,1,
     &                  ga_svd_tol)) then
c       See comment just before the ga_svd_solve_seq call to
c       understand these choices.
        if ((100*maxsub).lt.n) then
          ga_svd_tol = 1d-7
        else
          ga_svd_tol = 1d-14
        endif 
      endif
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) n, nvec, maxsub, maxiter, tol, util_wallsec()
 1       format(//,'Iterative solution of linear equations',/,
     $        '  No. of variables', i9,/,
     $        '  No. of equations', i9,/,
     $        '  Maximum subspace', i9,/,
     $        '        Iterations', i9,/,
     $        '       Convergence', 1p,d9.1,/,
     $        '        Start time', 0p,f9.1,/)
         call util_flush(6)
      end if
c
      if (.not. ga_create(MT_DBL, n, maxsub, 'lkain: Y', 
     $     0, 0, g_y))
     $     call errquit('lkain: failed allocating subspace', maxsub,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, maxsub, 'lkain: Ay', 
     $     0, 0, g_Ay))
     $     call errquit('lkain: failed allocating subspace2', maxsub,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, nvec, 'lkain: Ax',
     $     0, 0, g_Ax))
     $     call errquit('lkain: failed allocating subspace3', nvec,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, nvec, 'lkain: r',
     $     0, 0, g_r))
     $     call errquit('lkain: failed allocating subspace4', nvec,
     &       GA_ERR)
      if (odiff) then
         if (.not. ga_create(MT_DBL, n, nvec, 'lkain: xold',
     $        0, 0, g_xold))
     $        call errquit('lkain: failed allocating subspace5', nvec,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, n, nvec, 'lkain: xold',
     $        0, 0, g_Axold))
     $        call errquit('lkain: failed allocating subspace6', nvec,
     &       GA_ERR)
         call ga_zero(g_xold)
         call ga_zero(g_Axold)
      end if
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_Ax)
      call ga_zero(g_r)
      call ga_sync()
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,2)
         call util_flush(6)
 2       format(/
     $        '   iter   nsub   residual    time',/,
     $        '   ----  ------  --------  ---------')
      end if
      nsub = 0
      converged = .false.
      do iter = 1, maxiter
         if (odiff) then
            call ga_add(1.0d0, g_x, -1.0d0, g_xold,  g_x)
         end if
         call product(acc,g_x, g_Ax)
         if (odiff) then
            call ga_add(1.0d0, g_Ax, 1.0d0, g_Axold, g_Ax)
            call ga_add(1.0d0, g_x,  1.0d0, g_xold,  g_x)
            call ga_copy(g_x, g_xold)
            call ga_copy(g_Ax, g_Axold)
         end if
         call ga_zero(g_r)
         call ga_sync()
         call ga_add(1.0d0, g_b, -1.0d0, g_Ax, g_r) ! The residual
         call ga_sync()
         call ga_maxelt(g_r, rmax)
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter, nsub+nvec, rmax, util_wallsec()
            call util_flush(6)
 3          format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1)
         end if
         if (rmax .lt. tol) then
            converged = .true.
            goto 100
         end if
         call precond(g_Ax,0.0d0)
         call precond(g_r,0.0d0)
         call ga_sync()
c
c     Copy the vectors to the subspace work area
c
         call ga_copy_patch('n', 
     $        g_Ax, 1, n, 1, nvec, 
     $        g_Ay, 1, n, nsub+1, nsub+nvec)
         call ga_copy_patch('n', 
     $        g_x, 1, n, 1, nvec, 
     $        g_y, 1, n, nsub+1, nsub+nvec)
         nsub = nsub + nvec
c
c     Form and solve the subspace equations using SVD in order
c     to manage near linear dependence in the subspace.
c     
         if (.not. ga_create(MT_DBL, nsub, nsub, 'lkain: A', 0, 0, g_a))
     $        call errquit('lkain: allocating g_a?', nsub, GA_ERR)
         if (.not. ga_create(MT_DBL, nsub, nvec, 'lkain: B', 0, 0,g_bb))
     $        call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
         if (.not. ga_create(MT_DBL, nsub, nvec, 'lkain: C', 0, 0, g_c))
     $        call errquit('lkain: allocating g_c?', nsub, GA_ERR)
         call ga_zero(g_a)
         call ga_zero(g_bb)
         call ga_zero(g_c)
         call ga_sync()
         call ga_dgemm('t','n',nsub,nsub,n,1.0d0,g_y,g_Ay,0.0d0,g_a)
         call ga_sync()
         call ga_dgemm('t','n',nsub,nvec,n,1.0d0,g_y,g_r,0.0d0,g_bb)
         call ga_sync()
         if (odebug) call ga_print(g_a)
         if (odebug) call ga_print(g_c)
c
c     The threshold used here should reflect the accuracy in the
c     products.  If very accurate products are used, then there is big
c     advantage for small cases (maxsub close to n) in using a very
c     small threshold in the SVD solve (e.g., 1e-14), but for more
c     realistic examples (maxsub << n) there is only a little
c     advantage and in the precence of real noise in the products
c     screening with a realistic threshold is important.
c
         call ga_svd_solve_seq(g_a,g_bb,g_c,ga_svd_tol)
         if (odebug) call ga_print(g_c)
c
c     Form and add the correction, in parts, onto the solution
c
         call ga_sync()
         call ga_dgemm('n','n',n,nvec,nsub,-1.0d0,g_Ay,g_c,1.0d0,g_r)
         if (odebug) then
            write(6,*) ' The update in the complement '
            call ga_print(g_r)
         end if
         call ga_sync()
         call ga_add(1.0d0, g_r, 1.0d0, g_x, g_x)
         call ga_sync()
         call ga_dgemm('n','n',n,nvec,nsub,1.0d0,g_y,g_c,0.0d0,g_r)
         if (odebug) then
            write(6,*) ' The update in the subspace '
            call ga_print(g_r)
         end if
         call ga_sync()
         call ga_add(1.0d0, g_r, 1.0d0, g_x, g_x)
         call ga_sync()
c
         if (.not. ga_destroy(g_a)) call errquit('lkain: a',0, GA_ERR)
         if (.not. ga_destroy(g_bb))call errquit('lkain: b',0, GA_ERR)
         if (.not. ga_destroy(g_c)) call errquit('lkain: c',0, GA_ERR)
c
c     Reduce the subspace as necessary
c
         if (nsub .eq. maxsub) then
            do isub = nvec+1, maxsub, nvec
               call ga_copy_patch('n', 
     $              g_Ay, 1, n, isub, isub+nvec-1, 
     $              g_Ax, 1, n, 1, nvec)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, nvec,
     $              g_Ay, 1, n, isub-nvec, isub-1)
c
               call ga_copy_patch('n', 
     $              g_y, 1, n, isub, isub+nvec-1, 
     $              g_Ax, 1, n, 1, nvec)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, nvec,
     $              g_y, 1, n, isub-nvec, isub-1)
            end do
            nsub = nsub - nvec
            call ga_sync()
         end if
c
      end do
 100  continue
c
      if (odiff) then
         if (.not. ga_destroy(g_xold)) call errquit('lkain: destroy',1,
     &       GA_ERR)
         if (.not. ga_destroy(g_Axold)) call errquit('lkain: destroy',2,
     &       GA_ERR)
      end if
      if (.not. ga_destroy(g_Ax)) call errquit('lkain: destroy',20,
     &       GA_ERR)
      if (.not. ga_destroy(g_Ay)) call errquit('lkain: destroy',3,
     &       GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit('lkain: destroy',4,
     &       GA_ERR)
      if (.not. ga_destroy(g_r)) call errquit('lkain: destroy',5,
     &       GA_ERR)
c
      if (.not. converged) call errquit('lkain: not converged',0,
     &       CALC_ERR)
c
      end
      subroutine ga_kain(
     $     g_x, 
     $     product, precond, 
     $     tol, 
     $     trustmin, trustmax, 
     $     maxsub, maxiter, 
     $     oprint)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
c
      integer g_x               ! [input/output] Initial guess/solution
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      double precision trustmin, trustmax ! [input] range to constrain trust radius
      integer maxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical oprint            ! [input] print flag
c
c     Solves the non-linear equations f(X)=0 for multiple vectors.
c
c     call product(acc, g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_kain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . optional shift (not used here but used by the diagonalizer)
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
      integer maxmaxsub
      parameter (maxmaxsub = 20)
      integer iter, n, nvec, nsub, isub, jsub, type
      integer g_y, g_Ay, g_Ax, g_delta, g_a, g_b, g_c
      double precision rmax, acc, trust
      double precision a(maxmaxsub,maxmaxsub), b(maxmaxsub), 
     $     c(maxmaxsub), csum
      logical converged
      logical odebug
c
      trust = trustmin
      acc = 0.01d0*tol
      if (maxsub .gt. maxmaxsub) maxsub = maxmaxsub
      odebug = util_print('debug lsolve', print_never) .and. 
     $     ga_nodeid().eq.0
c
      call ga_inquire(g_x, type, n, nvec)
      if (nvec .ne. 1) call errquit('kain: nvec?', nvec, GA_ERR)
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) n, maxsub, tol, trustmin, trustmax, util_wallsec()
 1       format(//,'Iterative solution of non-linear equations',/,
     $        '  No. of variables', i9,/,
     $        '  Maximum subspace', i9,/,
     $        '       Convergence', 1p,d9.1,/,
     $        '     Trust min/max', 0p,2f6.2,
     $        '        Start time', 0p,f9.1,/)
         call util_flush(6)
      end if
c
      if (.not. ga_create(MT_DBL, n, maxsub, 'kain: Y', 
     $     0, 0, g_y))
     $     call errquit('kain: failed allocating subspace', maxsub,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, maxsub, 'kain: Ay', 
     $     0, 0, g_Ay))
     $     call errquit('kain: failed allocating subspace2', maxsub,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, 1, 'kain: Ax',
     $     0, 0, g_Ax))
     $     call errquit('kain: failed allocating subspace3', 1, GA_ERR)
      call ga_zero(g_y)
      call ga_zero(g_Ay)
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,2)
 2       format(/
     $        '   iter   nsub   residual    time',/,
     $        '   ----  ------  --------  ---------')
      end if
      nsub = 0
      converged = .false.
      do iter = 1, maxiter
         call ga_zero(g_Ax)
         call product(acc, g_x, g_Ax)
         call ga_maxelt(g_Ax, rmax)
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter, nsub+1, rmax, util_wallsec()
 3          format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1)
         end if
         if (rmax .lt. tol) then
            converged = .true.
            goto 100
         end if
c
c     Copy the vectors to the subspace work area
c
         call precond(g_Ax,0.0d0)
         call ga_copy_patch('n', 
     $        g_Ax, 1, n, 1, 1, 
     $        g_Ay, 1, n, nsub+1, nsub+1)
         call ga_copy_patch('n', 
     $        g_x, 1, n, 1, 1, 
     $        g_y, 1, n, nsub+1, nsub+1)
         nsub = nsub + 1
c
c     Not converged ... make the update
c
         g_delta = g_Ax         ! A reminder that these two are aliased
         call ga_scale(g_delta, -1.0d0)
c
         if (iter .gt. 1) then
c
c     Form the reduced space matrix and RHS
c
            call ga_local_mdot(n, nsub, nsub, a, maxmaxsub, g_y, g_Ay)
            do isub = 1, nsub-1
               b(isub) = -(a(isub,nsub) - a(nsub,nsub))
            end do
            do isub = 1, nsub-1
               do jsub = 1, nsub-1
                  a(isub,jsub) = a(isub,jsub)
     $                 - a(nsub,jsub) - a(isub,nsub) + a(nsub,nsub)
               end do
            end do
c
c     Solve the subspace equations (lazily using existing GA routine)
c     
            if (.not. ga_create(MT_DBL,nsub-1,nsub-1,'kain: A',
     $           nsub-1,nsub-1,g_a))
     $           call errquit('kain: allocating g_a?', nsub, GA_ERR)
            if (.not. ga_create(MT_DBL,nsub-1,1,'kain: B',nsub-1,1,g_b))
     $           call errquit('kain: allocating g_bb?', nsub, GA_ERR)
            if (.not. ga_create(MT_DBL,nsub-1,1,'kain: C',nsub-1,1,g_c))
     $           call errquit('kain: allocating g_c?', nsub, GA_ERR)
            if (ga_nodeid() .eq. 0) then
               call ga_put(g_a, 1, nsub-1, 1, nsub-1, a, maxmaxsub)
               call ga_put(g_b, 1, nsub-1, 1, 1, b, 1)
            end if
            call ga_sync
c
            call ga_svd_solve_seq(g_a,g_b,g_c,1d-14)
c
            if (odebug) call ga_print(g_c)
            if (ga_nodeid() .eq. 0) 
     $           call ga_get(g_c, 1, nsub-1, 1, 1, c, 1)
            call ga_brdcst(1, c, mdtob(nsub-1), 0)
            write(6,*) ' KAIN SUBSPACE COEFFS'
            call output(c, 1, nsub-1, 1, 1, nsub-1, 1, 1)
            call ga_sync
            if (.not. ga_destroy(g_a)) call errquit('kain: a',0, GA_ERR)
            if (.not. ga_destroy(g_b)) call errquit('kain: b',0, GA_ERR)
            if (.not. ga_destroy(g_c)) call errquit('kain: c',0, GA_ERR)
c     
c     Form the correction
c     
            csum = 0.0d0
            do isub = 1, nsub-1
               csum = csum + c(isub)
               call ga_add_patch( c(isub),  g_y, 1, n, isub, isub, 
     $              1.0d0, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
               call ga_add_patch(-c(isub), g_Ay, 1, n, isub, isub, 
     $              1.0d0, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
            end do
            call ga_add_patch(-csum,  g_y, 1, n, nsub, nsub, 
     $           1.0d0, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
            call ga_add_patch( csum,  g_Ay, 1, n, nsub, nsub, 
     $           1.0d0, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
         endif
c
c     Step restriction
c
         call ga_maxelt(g_delta, rmax)
         if (rmax .gt. trust) then
            if (oprint) write(6,*) ' RESTRICTION ', rmax, trust
            call ga_scale(g_delta, trust/rmax)
         end if
c
         call ga_add(1.0d0, g_delta, 1.0d0, g_x, g_x)
c     
c     Reduce the subspace as necessary (note g_delta=g_Ax destroyed)
c     
         if (nsub .eq. maxsub) then
            do isub = 2, maxsub
               call ga_copy_patch('n', 
     $              g_Ay, 1, n, isub, isub, 
     $              g_Ax, 1, n, 1, 1)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, 1,
     $              g_Ay, 1, n, isub-1, isub-1)
c     
               call ga_copy_patch('n', 
     $              g_y, 1, n, isub, isub, 
     $              g_Ax, 1, n, 1, 1)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, 1,
     $              g_y, 1, n, isub-1, isub-1)
            end do
            nsub = nsub - 1
         end if
c
      end do
 100  continue
c
      if (.not. ga_destroy(g_Ax)) call errquit('kain: destroy',20,
     &       GA_ERR)
      if (.not. ga_destroy(g_Ay)) call errquit('kain: destroy',3,
     &       GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit('kain: destroy',4, GA_ERR)
c
      if (.not. converged) call errquit('kain: not converged',0,
     &       CALC_ERR)
c
      end
      subroutine ga_svd_solve_seq(g_a, g_b, g_x, tol)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer g_a, g_b, g_x
      double precision tol
c
c     Solve for X from the linear equations
c
c     A*X = B
c
c     A(m,n)*X(n,nvec) = B(m,nvec)
c
c     Where A is a general real matrix (not necessarily square, or
c     symmetric, or full rank) and X and B are matrices with one or more
c     columns representing the solutions and right hand sides.  Singular
c     values of A less than tol are neglected.  X is returned.
c
c     If the SVD of A is U*values*VT, then the solution
c     is of the form
c
c     V*(1/values)*UT*B
c
c     where the reciprocal of values less than tol are neglected.
c     
      integer m,n,nn,type,nvec,nsing,l_val, k_val,g_u,g_vt,i,g_tmp
      logical oprint
c
      oprint = util_print('debug svdsolve', print_high) .and.
     $     ga_nodeid().eq.0
c
      call ga_inquire(g_a, type, m, n)
      call ga_inquire(g_b, type, nn, nvec)
      if (nn .ne. n) call errquit('gasvdsol: b does not conform',nn,
     &       GA_ERR)
      nsing = min(m,n)
      if (.not. ma_push_get(MT_DBL, nsing, 'gasvdsol', l_val, k_val))
     $     call errquit('gasvdsol: val',nsing, MA_ERR)
      if (.not. ga_create(MT_DBL,m,nsing,'gasvd',0,0,g_u))
     $     call errquit('gasvdsol: u',m*nsing, GA_ERR)
      if (.not. ga_create(MT_DBL,nsing,n,'gasvd',0,0,g_vt))
     $     call errquit('gasvdsol: u',nsing*n, GA_ERR)
      if (.not. ga_create(MT_DBL,nsing,nvec,'gasvd',0,0,g_tmp))
     $     call errquit('gasvdsol: tmp',nsing*nvec, GA_ERR)
      call ga_zero(g_tmp)
c
      call ga_svd_seq(g_a, g_u,g_vt,dbl_mb(k_val))
c
      do i = 0, nsing-1
         if (dbl_mb(k_val+i) .lt. tol) then
            if (ga_nodeid() .eq. 0 .and. oprint) then
              write(6,*) ' neglecting ', i+1, dbl_mb(k_val+i)
            endif
            dbl_mb(k_val+i) = 0.0d0
         else
            dbl_mb(k_val+i) = 1.0d0/dbl_mb(k_val+i)
         end if
      end do
c
      call ga_dgemm('t','n',nsing,nvec,m,1.0d0,g_u,g_b,0.0d0,g_tmp)
      call ga_scale_lh(g_tmp,dbl_mb(k_val))
      call ga_zero(g_x)
      call ga_dgemm('t','n',n,nvec,nsing,1.0d0,g_vt,g_tmp,0.0d0,g_x)
c
      if (.not. ga_destroy(g_tmp)) call errquit('gasvdsol: des',1,
     &       GA_ERR)
      if (.not. ga_destroy(g_u)) call errquit('gasvdsol: des',2,
     &       GA_ERR)
      if (.not. ga_destroy(g_vt)) call errquit('gasvdsol: des',3,
     &       GA_ERR)
      if (.not. ma_pop_stack(l_val)) call errquit('gasvdsol: pop',4,
     &       GA_ERR)
c
      end
      subroutine ga_svd_seq(g_a, g_u, g_vt, values)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_a, g_u, g_vt
      double precision values(*)
c     
c     Perform SVD on rectangular matrix
c
c     nsing = min(n,m)
c     g_a(m,n)      --- input matrix
c     g_u(m,nsing)  --- left singular vectors (output)
c     g_vt(nsing,n) --- right singular vectors transposed (output)
c     values(nsing) --- singular values (output)
c
c     A = U*values*VT
c
c     A possible parallel algorithm is to diagonalize ATA to get
c     V and AAT to get U --- both have values**2 as eigenvalues.
c
      integer n, m, type, l_a, k_a, l_u, k_u, l_vt, k_vt, 
     $     l_work, k_work, lwork, info, nsing
c     
      call ga_inquire(g_a, type, m, n)
      nsing = min(m,n)
      if (ga_nodeid() .eq. 0) then
         lwork = 10*max(m,n)
         if (.not. ma_push_get(MT_DBL, m*n, 'gasvd', l_a, k_a))
     $        call errquit('gasvd: a',m*n, MA_ERR)
         if (.not. ma_push_get(MT_DBL, m*nsing, 'gasvd', l_u, k_u))
     $        call errquit('gasvd: u',m*nsing, MA_ERR)
         if (.not. ma_push_get(MT_DBL, nsing*n, 'gasvd', l_vt, k_vt))
     $        call errquit('gasvd: vt',nsing*n, MA_ERR)
         if (.not. ma_push_get(MT_DBL, lwork, 'gasvd', l_work, k_work))
     $        call errquit('gasvd: work',lwork, MA_ERR)
c
         call ga_get(g_a, 1, m, 1, n, dbl_mb(k_a), m)
c
         call dgesvd('s','s',m,n,dbl_mb(k_a),m,values,
     $        dbl_mb(k_u),m,dbl_mb(k_vt),nsing,
     $        dbl_mb(k_work),lwork,info)
         if (info .ne. 0) call errquit('gasvd: failed', info, MEM_ERR)
c
         call ga_put(g_u,  1, n,     1, nsing, dbl_mb(k_u),  n)
         call ga_put(g_vt, 1, nsing, 1, m,     dbl_mb(k_vt), n)
c
         if (.not. ma_chop_stack(l_a)) call errquit('gasvd ma',0,
     &       MA_ERR)
      end if
      call ga_sync()
      call ga_brdcst(1,values,n*8,0)
      call ga_sync()
c     
      end
