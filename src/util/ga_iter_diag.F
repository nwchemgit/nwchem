      logical function ga_iter_diag(n, nroot, maxiter, maxsub, tol,
     $     precond, product, oprint, eval0, g_evec, eval, rnorm, iter)
C$Id: ga_iter_diag.F,v 1.8 1999-06-16 00:10:28 d3e129 Exp $
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
      integer n                 ! Matrix dimension
      integer nroot             ! No. of eigen vectors sought
      integer maxiter           ! Maximum no. of iterations
      integer maxsub            ! Maximum dimension of iterative subspace
      double precision tol      ! Required norm of residual
      external precond          ! Preconditioner
      external product          ! Matrix-vector product
      logical oprint            ! True for printing of convergence to unit 6
      double precision eval0    ! Estimate of lowest eval
      integer g_evec            ! n by nroot GA for guess and final
      double precision eval(nroot) ! Returns eigen values
      double precision rnorm(nroot) ! Returns norm of residual for each root
      integer iter              ! Returns actual no. of iterations used
c     
c     Solve the eigenvalue equation Ax = ex with x in GA and a routine
c     (product) to form a matrix vector product to a required precision.
c     Return true if converged, false otherwise. Rnorm returns the actual
c     attained precision for each root.
c     
c     The block-Davidson-like algorithm solves for the best solution for
c     each eigenroottor in the iterative subspace (x_i, i = 1, k) with
c     _     _          _                     _
c     A y = Sye, where A_ij = x_i*A*x_j, and S_ij = x_i * x_j
c     
c     The matrix vector products A*x_i are performed by the user provided
c     routine product to a precision specified by this routine
c     (currently products are performed one at a time).
c     
c     The best solution within the iterative subspace is then
c     
c     x = sum(i) x_i * y_i
c     
c     New expansion vectors are added by multiplying the residual
c     
c     r = (A - shift*I) x
c     
c     with some approximation (P) to the inverse of (A-shift*I).  This 
c     preconditioning is performed by the user provided routine precond().
c     If eval0 is a hard zero then shift is chosen as the current estimate
c     for the eigenvalue that the next update strives to improve.
c     Otherwise shift is fixed as eval0 which is appropriate for
c     convergence to a known energy spectrum from some poor initial guess.
c     
c     The program cyles through the lowest nroot roots updating each that
c     does not yet satisfy the convergence criterion which is
c     
c     rnorm(root) = ||r|| < tol 
c     
c     On input the global array x(n,nroot) should contain either an initial
c     guess at the eigen vectors or zeroes.  If any vector is zero then
c     random numbers are used.
c     
c     subroutine product(precision, g_x, g_ax)
c     
c     computes the product A*x to the specified precision (absolute
c     magnitude error in any element of the product) returning the result
c     in the g_ax.
c     
c     subroutine precond(g_r, shift)
c     
c     Apply an approximation (P) to the inverse of (A - shift I) to the 
c     vector in g_r overwriting g_r with the result.
c     
c     If the initial guess is zero no redundant matrix product is formed
c     
      integer g_x               ! GA to hold expansion vectors
      integer g_p               ! GA to hold product vectors
      integer g_v               ! GA to hold temporary vector
      integer g_w               ! GA to hold temporary product vector
c
      double precision test_tol
      double precision shift
      double precision prod_acc
      logical converged
      integer maxdim
      parameter (maxdim = 10)   ! Maximum size of iterative subspace
      double precision aa(maxdim,maxdim), ss(maxdim,maxdim),
     $     yy(maxdim,maxdim), ee(maxdim), fv1(maxdim), fv2(maxdim)
      double precision xnorm
      integer i, j, nsub, info, cur_root, iroot
      logical odebug
      data odebug /.false./
c     
      if (maxsub .gt. maxdim) call errquit
     $     ('ga_iter_diag: maxsub > maxdim', maxdim)
      if (maxsub .lt. nroot+2) call errquit
     $     ('ga_iter_diag: subspace too small', nsub)
c     
c     Allocate and zero the temporary global arrays
c     
*ga:1:0
      if (.not. ga_create(MT_DBL, n, maxsub, 'ga_iter_diag: x', 0,
     $     maxsub, g_x))
     $     call errquit('ga_iter_diag: ga_create of x', maxsub)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, maxsub, 'ga_iter_diag: p', 0,
     $     maxsub, g_p))
     $     call errquit('ga_iter_diag: ga_create of p', maxsub)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 0, 'ga_iter_diag: v', 1, 1,
     $     g_v)) call errquit('ga_iter_diag: ga_create of v', n)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 0, 'ga_iter_diag: w', 1, 1,
     $     g_w)) call errquit('ga_iter_diag: ga_create of w', n)
      call ga_zero(g_x)
      call ga_zero(g_p)
c     
c     Compute matrix-vector products to slightly beyond the
c     required tolerance
c     
      prod_acc = tol * 0.01d0
c     
c     Form initial guess vectors (if necessary), orthonormalize
c     and perform initial matrix-vector products
c     
      do i = 1, nroot
         call ga_copy_patch('n', g_evec, 1, n, i, i, g_v, 1, n, 1, 1)
         xnorm = sqrt(ga_ddot(g_v, g_v))
         if (xnorm .eq. 0.0d0) then
            call ga_ran_fill(g_v, 1, n, 1, 1)
            do iter = 1,4
               call precond(g_v, eval0) ! Eliminate high-energy modes
            enddo
         endif
         call ga_orthog_vec(n, i-1, g_x, g_v, 1)
         xnorm = sqrt(ga_ddot(g_v, g_v))
         call ga_dscal(g_v, 1.0d0/xnorm)
c     
         call ga_zero(g_w)
         call product(prod_acc, g_v, g_w)
         call ga_copy_patch('n', g_v, 1, n, 1, 1, g_x, 1, n, i, i)
         call ga_copy_patch('n', g_w, 1, n, 1, 1, g_p, 1, n, i, i)
      enddo
c     
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) tol, prod_acc
 1       format(/,12x,'----- iterative diag (tol:',1p,e8.1,' prod-acc:',
     $        e8.1,') -----'/)
         write(6,2) 
 2       format(
     $        7x, ' iter nsub  cur root   eigenvalue      residual',
     $        '    time'/
     $        7x, ' ---- ---- ---- ---- --------------- -----------',
     $        ' -------')
         call util_flush(6)
      endif
c     
      test_tol = 0.5d0
      cur_root  = 0             ! Current vector being updated
      nsub = nroot              ! No. of vectors within subspace
      do iter = 1, maxiter
c     
c     Now have nsub vectors and products in g_x and g_p.
c     Form the reduced space matrices.
c     
         do i = 1, nsub
            do j = 1, i
               aa(j,i) = ga_ddot_patch(
     $              g_x, 'n', 1, n, i, i, 
     $              g_p, 'n', 1, n, j, j)
               aa(i,j) = aa(j,i)
               ss(j,i) = ga_ddot_patch(
     $              g_x, 'n', 1, n, i, i, 
     $              g_x, 'n', 1, n, j, j)
               ss(i,j) = ss(j,i)
            enddo
         enddo
         if (odebug .and. ga_nodeid().eq.0) then
            write(6,*) ' ga_iter_diag: Reduced space matrix'
            call output(aa, 1, nsub, 1, nsub, maxdim, nsub, 1)
            write(6,*) ' ga_iter_diag: Reduced space overlap'
            call output(ss, 1, nsub, 1, nsub, maxdim, nsub, 1)
         endif
c     
c     Solve the reduced space equations (Ay = Sye)
c     
         call rsg(maxdim, nsub, aa, ss, ee, 1, yy, fv1, fv2, info)
         if (info .ne. 0) call errquit
     $        ('ga_iter_diag: rsg failed', info)
         if (odebug .and. ga_nodeid().eq.0) then
            write(6,*) ' ga_iter_diag: Reduced space eigenvectors'
            call output(yy, 1, nsub, 1, nsub, maxdim, nsub, 1)
            write(6,*) ' ga_iter_diag: Reduced space eigenvalues'
            call output(ee, 1, nsub, 1, 1, nsub, 1, 1)
         endif
         call dcopy(nroot, ee, 1, eval, 1)
c     
c     Resum the products and vectors to the eigen-space
c     x(i,t) = sum(s) x(i,s) * y(s,t)
c     
         call ga_mix(g_x, n, nsub, yy, maxdim)
         call ga_mix(g_p, n, nsub, yy, maxdim)
c     
c     Form all the residuals in g_evec
c     
         call ga_zero(g_evec)
         do iroot = 1, nroot
            call ga_dadd_patch(
     $           1.0d0, g_p, 1, n, iroot, iroot, 
     $           1.0d0, g_evec, 1, n, iroot, iroot,
     $           g_evec, 1, n, iroot, iroot)
            call ga_dadd_patch(
     $           -eval(iroot), g_x, 1, n, iroot, iroot,
     $           1.0d0, g_evec, 1, n, iroot, iroot,
     $           g_evec, 1, n, iroot, iroot)
            rnorm(iroot) = sqrt(ga_ddot_patch(
     $           g_evec, 'n', 1, n, iroot, iroot,
     $           g_evec, 'n', 1, n, iroot, iroot))
         enddo
c     
c     Print
c     
         if (oprint .and. ga_nodeid().eq.0) then
            do iroot = 1, nroot
               write(6,3) iter, nsub, cur_root, iroot, eval(iroot),
     $              rnorm(iroot), tcgtime()
 3             format(7x,4i5,1p,e16.9,e12.4,0p,f8.1)
            enddo
            call util_flush(6)
         endif
c     
c     Determine the next vector to update
c     
 20      converged = .true.
         cur_root = 0           ! Comment this out to cycle
         do iroot = 1, nroot
            cur_root = mod(cur_root,nroot) + 1
            if (rnorm(iroot) .gt. test_tol) then
               converged = .false.
               goto 10
            endif
         enddo
         
 10      if (converged) then
            if (test_tol .le. tol) then
               goto 1000
            else
               test_tol = max(test_tol * 0.25d0, tol)
               goto 20
            endif
         endif
c     
         if (iter.eq.maxiter) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' ga_iter_diag: exceeded max iters ... ',
     $              'aborting solve'
               call util_flush(6)
            endif
            goto 1000
         endif
c     
c     Form a new update vector by preconditioning and orthonormalizing
c     (have to normalize to ensure precision since don't have assurance
c     of monotonic convergence unless initial vectors are close to
c     the correct eigenvectors)
c     
         call ga_copy_patch('n', g_evec, 1, n, cur_root, cur_root, 
     $        g_v, 1, n, 1, 1)
         if (eval0 .eq. 0.0d0) then
            shift = eval(iroot)
         else
            shift = eval0
         endif
         call precond(g_v, shift)
         call ga_orthog_vec(n, nsub, g_x, g_v, 1)
         xnorm = sqrt(ga_ddot(g_v, g_v))
         call ga_dscal(g_v, 1.0d0/xnorm)
c     
c     Form new product and shove vectors into the array of vectors
c     
         call ga_zero(g_w)
         call product(prod_acc, g_v, g_w)
         if (nsub .lt. maxsub) nsub = nsub + 1
         call ga_copy_patch('n', g_v, 1, n, 1, 1, g_x, 1, n, nsub, nsub)
         call ga_copy_patch('n', g_w, 1, n, 1, 1, g_p, 1, n, nsub, nsub)
c     
      enddo
c     
c     Should never actually fall out of the bottom
c     
      call errquit('ga_iter_diag: fell out of loop bottom?', 0)
c     
 1000 continue
c     
c     Copy solutions into output arrays
c     
      call ga_copy_patch('n', 
     $     g_x, 1, n, 1, nroot, 
     $     g_evec, 1, n, 1, nroot)
      call dcopy(nroot, ee, 1, eval, 1)
c     
      ga_iter_diag = converged
c     
      if (.not. ga_destroy(g_v)) call errquit
     $     ('ga_iter_diag: failed to destory v', 0)
      if (.not. ga_destroy(g_w)) call errquit
     $     ('ga_iter_diag: failed to destory w', 0)
      if (.not. ga_destroy(g_x)) call errquit
     $     ('ga_iter_diag: failed to destory x', 0)
      if (.not. ga_destroy(g_p)) call errquit
     $     ('ga_iter_diag: failed to destory p', 0)
c     
      end
