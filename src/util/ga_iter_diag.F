      logical function ga_iter_diag(n, maxiter, maxsub, tol,
     $     precond, product, oprint, g_evec, eval, acc, iter)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
      integer n                 ! Problem dimension
      integer maxiter           ! Maximum no. of iterations
      integer maxsub            ! Maximum dimension of iterative subspace
      double precision tol      ! Required norm of residual
      external precond          ! Preconditioner
      external product          ! Matrix-vector product
      logical oprint            ! True for printing of convergence to unit 6
      integer g_evec            ! GA with initial guess and returns vector
      double precision eval     ! Initial guess for eigen value and returns
      double precision acc      ! Returns precision actually achieved
      integer iter              ! Returns actual no. of iterations used
c
c     To be extended soon to multiple eigen vectors.
c
c     Solve the eigenvalue equation Ax = ex with x in GA and a routine
c     (product) to form a matrix vector product to a required precision.
c     Return true if converged, false otherwise. Acc returns the actual
c     attained precision.
c     
c     The Davidson-like algorithm solves for the best solution in the 
c     iterative subspace (x_i, i = 1, k) with
c     _     _          _                     _
c     A y = Sye, where A_ij = x_i*A*x_j, and S_ij = x_i * x_j
c
c     The matrix vector products A*x_i are performed by the user provided
c     routine product to a precision specified by this routine.
c
c     The best solution within the iterative subspace is then
c
c     x = sum(i) x_i * y_i
c
c     New expansion vectors are added by multiplying the residual
c
c     r = (A*x - ex) 
c
c     with some approximation (P) to the inverse of A.  This preconditioning
c     is performed by the user provided routine precond().
c
c     The convergence condition is that
c
c     acc = ||r|| / abs(initial guess eval) < tol 
c
c     (if abs(initial guess eval) is zero then just use ||r||).
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
      integer g_w               ! GA to hold temporary vector
      double precision prod_acc
      logical converged
      integer maxdim
      parameter (maxdim = 50)   ! Maximum size of iterative subspace
      double precision aa(maxdim,maxdim), ss(maxdim,maxdim),
     $     yy(maxdim,maxdim), ee(maxdim), fv1(maxdim), fv2(maxdim) 
      double precision rx, xx
      double precision rnorm, rnorm_prev, bnorm, convg_ratio
      integer i, j, nsub, info
      logical odebug
      data odebug /.false./
c
      if (maxsub .gt. maxdim) call errquit
     $     ('ga_iter_diag: maxsub > maxdim', maxdim)
c
c     Right now there must be an initial guess ... normalize it.
c
      rnorm = ga_ddot(g_evec, g_evec)
      if (rnorm .eq. 0.0d0) call errquit
     $     (' ga_iter_diag: must have non-zero initial guess', 0)
      call ga_dscal(g_evec, 1.0d0/sqrt(rnorm))
c
c     Required precision in matrix-vector products is probably the
c     eigen-value * the tolerance + 3 extra guard digits (to try
c     to prevent variational collapse).
c
      if (eval .ne. 0.0d0) then
         prod_acc = abs(eval) * tol * 0.001d0
         bnorm = abs(eval)
      else
         prod_acc = tol * 0.1d0
         bnorm = 1.0d0
      endif
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) tol, prod_acc
 1       format(/,12x,'----- iterative diag (tol:',1p,e8.1,' prod-acc:',
     $        e8.1,') -----'/)
         write(6,2) 
 2       format(7x, ' iter nbas    eigenvalue      rnorm   ',
     $        '  rnorm-scal   ratio     time'/
     $          7x, ' ---- ---- --------------- -----------',
     $        ' ----------- --------- -------')
         call util_flush(6)
      endif
c
c     Allocate and zero the temporary global arrays
c
      if (.not. ga_create(MT_DBL, n, maxsub, 'ga_iter_diag: x', 1,
     $     maxsub, g_x))
     $     call errquit('ga_iter_diag: ga_create of x', maxsub)
      if (.not. ga_create(MT_DBL, n, maxsub, 'ga_iter_diag: p', 1,
     $     maxsub, g_p))
     $     call errquit('ga_iter_diag: ga_create of p', maxsub)
      if (.not. ga_create(MT_DBL, n, 1, 'ga_iter_diag: w', 1, 1,
     $     g_w)) call errquit('ga_iter_diag: ga_create of w', maxsub)
      call ga_zero(g_x)
      call ga_zero(g_p)
c
      nsub = 1                  ! Points to next free subspace vector
      do iter = 1, maxiter
c
c     Arrive here with the arrays g_x and g_p containing nsub-1 vectors
c     and corresponding products.  The current trial vector is in g_evec.
c
c     Form A * g_evec in g_w
c
         call ga_zero(g_w)
         call product(prod_acc, g_evec, g_w)
c     
c     Copy the vector and product into the arrays
c     
         call ga_copy_patch('n', 
     $        g_w, 1, n, 1, 1, 
     $        g_p, 1, n, nsub, nsub) 
         call ga_copy_patch('n', 
     $        g_evec, 1, n, 1, 1, 
     $        g_x, 1, n, nsub, nsub) 
c
c     Now have nsub vectors and products in g_x and g_p
c     and the rhs back in g_w. Form the reduced space matrices
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
c
         eval = ee(1)
c
c     Form the residual in g_evec
c
         call ga_zero(g_evec)
         do i = 1, nsub
            call ga_dadd_patch(
     $           yy(i,1), g_p, 1, n, i, i, 
     $           1.0d0, g_evec, 1, n, 1, 1,
     $           g_evec, 1, n, 1, 1)
            call ga_dadd_patch(
     $           -eval*yy(i,1), g_x, 1, n, i, i, 
     $           1.0d0, g_evec, 1, n, 1, 1,
     $           g_evec, 1, n, 1, 1)
         enddo
c
c     Print and determine convergence or too many iterations
c     
         rnorm = sqrt(ga_ddot(g_evec, g_evec))
         acc = rnorm/bnorm
         converged = acc .lt. tol
         if (iter.gt.1) then
            convg_ratio = rnorm/rnorm_prev
         else
            convg_ratio = 0.0d0
         endif
         rnorm_prev = rnorm
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter,nsub,eval,rnorm,acc,convg_ratio,tcgtime()
 3          format(7x,2i5,1p,e16.9,2e12.4,0p,f10.4,f8.1)
            call util_flush(6)
         endif
         if (converged) goto 1000
         if (iter.eq.maxiter) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' ga_iter_diag: exceeded max iters ... ',
     $              'aborting solve'
               call util_flush(6)
            endif
            goto 1000
         endif
c
c     Form a new update vector by preconditioning and orthogonalizing
c     to the previous expansion space
c
         call precond(g_evec, eval)
 20      rnorm = sqrt(ga_ddot(g_evec,g_evec))
         call ga_dscal(g_evec, 1.0d0/rnorm)
         do i = 1, nsub
            rx = ga_ddot_patch(
     $           g_evec, 'n', 1, n, 1, 1,
     $           g_x, 'n', 1, n, i, i)
            xx = ga_ddot_patch(
     $           g_x, 'n', 1, n, i, i, 
     $           g_x, 'n', 1, n, i, i)
            call ga_dadd_patch(
     $           1.0d0, g_evec, 1, n, 1, 1,
     $           -rx/xx, g_x, 1, n, i, i,
     $           g_evec, 1, n, 1, 1)
         enddo
         rnorm = sqrt(ga_ddot(g_evec,g_evec))
         call ga_dscal(g_evec, 1.0d0/rnorm)
         if (rnorm .lt. 0.9d0) goto 20
c
         if (nsub .lt. maxsub) then
            nsub = nsub + 1
         else
c
c     If have run out of space have to compress the list by summing
c     the oldest two vectors and shuffling down.
c     
            rnorm = sqrt(yy(1,1)**2 + yy(2,1)**2)
            call ga_dscal_patch(g_x, 1, n, 1, 1, yy(1,1)/rnorm)
            call ga_dscal_patch(g_p, 1, n, 1, 1, yy(1,1)/rnorm)
            call ga_dadd_patch(
     $           yy(2,1)/rnorm, g_x, 1, n, 2, 2,
     $           1.0d0, g_x, 1, n, 1, 1,
     $           g_x, 1, n, 1, 1)
            call ga_dadd_patch(
     $           yy(2,1)/rnorm, g_p, 1, n, 2, 2,
     $           1.0d0, g_p, 1, n, 1, 1,
     $           g_p, 1, n, 1, 1)
            do i = 3,nsub
               call ga_copy_patch('n', 
     $              g_x, 1, n, i, i, 
     $              g_w, 1, n, 1, 1)
               call ga_copy_patch('n', 
     $              g_w, 1, n, 1, 1, 
     $              g_x, 1, n, i-1, i-1)
               call ga_copy_patch('n', 
     $              g_p, 1, n, i, i, 
     $              g_w, 1, n, 1, 1)
               call ga_copy_patch('n', 
     $              g_w, 1, n, 1, 1, 
     $              g_p, 1, n, i-1, i-1)
            enddo
         endif
      enddo
c
c     Should never actually fall out of the bottom
c
      call errquit('ga_iter_diag: fell out of loop bottom?', 0)
c
 1000 continue
      call ga_zero(g_evec)
      do i = 1, nsub
         call ga_dadd_patch(
     $        yy(i,1), g_x, 1, n, i, i, 
     $        1.0d0, g_evec, 1, n, 1, 1,
     $        g_evec, 1, n, 1, 1)
      enddo
      eval = ee(1)
c
      ga_iter_diag = converged
c
      if (.not. ga_destroy(g_w)) call errquit
     $     ('ga_iter_diag: failed to destory w', 0)
      if (.not. ga_destroy(g_x)) call errquit
     $     ('ga_iter_diag: failed to destory x', 0)
      if (.not. ga_destroy(g_p)) call errquit
     $     ('ga_iter_diag: failed to destory p', 0)
c
      end
