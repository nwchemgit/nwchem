      logical function ga_iter_lsolve(n, maxiter, maxsub, tol,
     $     precond, product, ostall, oprint, g_b, g_soln, acc, iter)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer n                 ! Problem dimension
      integer maxiter           ! Maximum no. of iterations
      integer maxsub            ! Maximum dimension of iterative subspace
      double precision tol      ! Required norm of residual / norm of rhs
      external precond          ! Preconditioner
      external product          ! Matrix-vector product
      logical ostall            ! If true abort if convergence is too slow
      logical oprint            ! True for printing of convergence to unit 6
      integer g_b               ! GA with RHS vector
      integer g_soln            ! GA with initial guess and returns solution
      double precision acc      ! Returns relative precision actually achieved
      integer iter              ! Returns actual no. of iterations used
c
c     To be extended soon to multiple RHS vectors ... Miles has this?
c
c     Solve the linear equation Ax = b with x and b in GAs and a routine
c     (product) to form a matrix vector product to a required precision.
c     Return true if converged, false otherwise. Acc returns the actual
c     attained precision.
c     
c     The least-squares algorithm solves for the best solution in the 
c     iterative subspace (x_i, i = 1, k) with
c     _     _         _                            _
c     A c = b , where A_ij = (A*x_i).(A*x_j) , and b_i = (A*x_i).b
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
c     r = (A*x - b) 
c
c     with some approximation (P) to the inverse of A.  This preconditioning
c     is performed by the user provided routine precond().
c
c     The convergence condition is that
c
c     acc = ||r|| / ||b|| < tol
c
c     subroutine product(precision, g_x, g_ax)
c
c     computes the product A*x to the specified precision (absolute
c     magnitude error in any element of the product) returning the result
c     in the g_ax.
c
c     subroutine precond(g_r, shift)
c
c     Apply an approximation (P) to the inverse of (A - shift*I) to the
c     vector in g_r overwriting g_r with the result.
c
c     If the initial guess is zero no redundant matrix product is formed
c
c     If (ostall) then if the reduction of error between sucessive 
c     iterations is too slow then terminate.
c
      integer g_x               ! GA to hold expansion vectors
      integer g_p               ! GA to hold product vectors
      integer g_w               ! GA vector for resummation only.
      double precision prod_acc
      logical converged
      integer maxdim
      parameter (maxdim = 100)   ! Maximum size of iterative subspace
      double precision aa(maxdim, maxdim), bb(maxdim), yy(maxdim)
      double precision rx, xx
      double precision rnorm, rnorm_prev, bnorm, convg_ratio
      double precision mone, one, zero
      parameter (mone = -1.0d0, one = 1.0d0, zero = 0.0d0)
      integer temp(maxdim), info ! For dgesv
      integer i, nsub
c     integer j
      logical odebug
      logical owarn
c
      odebug = util_print('debug lsolve',print_never)
      owarn  = util_print('warnings',print_low) .and.
     $     (ga_nodeid() .eq. 0)
c
      if (maxsub .gt. maxdim) call errquit
     $     ('ga_iter_lsolve: maxsub > maxdim', maxdim, INPUT_ERR)
c
c     Required precision in matrix-vector products is about the maximum
c     element in the rhs vector * the tolerance + 3 extra guard digits
c     (1 digit only causes allied test to have problems)
c     (2 digits causes 6-31g* zolite clusters to have problems)
c
c     Experiment with rms value rather than max
c
      bnorm = sqrt(ga_ddot(g_b, g_b))
      prod_acc = bnorm / sqrt(dble(n))
*     call ga_maxelt(g_b, prod_acc)
cold      prod_acc = prod_acc * tol * 0.001d0
      prod_acc = prod_acc * tol * 0.00001d0
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) tol, prod_acc
    1    format(/,12x,'----- linear solve (tol:',1p,e8.1,' prod-acc:',
     $        e8.1,') -----'/)
         write(6,2) 
    2    format(14x, ' iter nbas    rnorm     rnorm-scal   ratio',
     $        '     time'/
     $        14x, ' ---- ---- ----------- ----------- ---------',
     $        ' -------')
         call util_flush(6)
      end if
c
c     Allocate and zero the temporary global arrays
c
*ga:1:0
      if (.not. ga_create(MT_DBL, n, 1, 'ga_iter_lsolve: w', 0,
     $     0, g_w))
     $     call errquit('ga_iter_lsolve: ga_create of w', 1, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, maxsub, 'ga_iter_lsolve: x', 0,
     $     maxsub, g_x))
     $     call errquit('ga_iter_lsolve: ga_create of x', maxsub,
     &       GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, n, maxsub, 'ga_iter_lsolve: p', 0,
     $     maxsub, g_p))
     $     call errquit('ga_iter_lsolve: ga_create of p', maxsub,
     &       GA_ERR)
      call ga_zero(g_x)
      call ga_zero(g_p)
c
c     If the initial guess is zero then use Pb
c
      rnorm = ga_ddot(g_soln, g_soln)
      if (rnorm .eq. zero) then
         call ga_copy(g_b, g_soln)
*         call precond(g_soln, zero)
         call ga_precond_project(g_soln, zero, precond, g_w)
      end if
c
      converged = .false.       ! get rid of compiler warning
      nsub = 1                  ! Points to next free subspace vector
      do iter = 1, maxiter
c
c     Arrive here with the arrays g_x and g_p containing nsub-1 vectors
c     and corresponding products.  The current trial vector is in g_soln
c     and the RHS vector is in g_b
c
c     Use g_b to hold the product vector to save space ... keep b where
c     the current solution vector would live in g_x
c
         call ga_mat_to_vec(g_b, 1, n, 1, 1, g_x, 1, nsub, 1.0d0, '=')
*         call ga_copy_patch('n', 
*     $        g_b, 1, n, 1, 1, 
*     $        g_x, 1, n, nsub, nsub) 
c     
         call ga_zero(g_b)
         call product(prod_acc, g_soln, g_b)
c     
c     Copy the product into the array of product vectors, b back to g_b
c     and the current trial vector back into the array of trial vectors
c     
         call ga_mat_to_vec(g_b, 1, n, 1, 1, g_p, 1, nsub, 1.0d0, '=')
*         call ga_copy_patch('n', 
*     $        g_b, 1, n, 1, 1, 
*     $        g_p, 1, n, nsub, nsub) 
         call ga_mat_to_vec(g_x, 1, n, nsub, nsub, g_b, 1, 1,1.0d0,'=')
*         call ga_copy_patch('n', 
*     $        g_x, 1, n, nsub, nsub,
*     $        g_b, 1, n, 1, 1)
         call ga_mat_to_vec(g_soln, 1, n, 1, 1, g_x, 1, nsub,1.0d0,'=')
*         call ga_copy_patch('n', 
*     $        g_soln, 1, n, 1, 1, 
*     $        g_x, 1, n, nsub, nsub) 
c
c     Now have nsub vectors and products in g_x and g_p
c     and the rhs back in g_b. Form the reduced space matrices
c
         call ga_local_mdot(n, nsub, nsub, aa, maxdim, g_p, g_p)
         call ga_local_mdot(n, nsub, 1,    bb, nsub,   g_p, g_b)
*         do i = 1, nsub
*            do j = 1, i
*               aa(j,i) = ga_ddot_patch(
*     $              g_p, 'n', 1, n, i, i, 
*     $              g_p, 'n', 1, n, j, j)
*               aa(i,j) = aa(j,i)
*            end do
*            bb(i) = ga_ddot_patch(
*     $           g_p, 'n', 1, n, i, i,
*     $           g_b, 'n', 1, n, 1, 1)
*         end do
         if (odebug .and. ga_nodeid().eq.0) then
            write(6,*) ' ga_iter_lsolve: Reduced space matrix '
            call output(aa, 1, nsub, 1, nsub, maxdim, nsub, 1)
            write(6,*) ' ga_iter_lsolve: Reduced space rhs '
            call output(bb, 1, nsub, 1, 1, nsub, 1, 1)
         end if
c
c     Solve the reduced space equations (AAyy = bb) ... solution in yy
c
         call dgesv(nsub, 1, aa, maxdim, temp, bb, maxdim, info)
         if (info .ne. 0) call errquit
     $        ('ga_iter_lsolve: dgesv failed', info, GA_ERR)
         call dcopy(nsub, bb, 1, yy, 1)
         if (odebug .and. ga_nodeid().eq.0) then
            write(6,*) ' ga_iter_lsolve: Reduced space solution '
            call output(yy, 1, nsub, 1, 1, nsub, 1, 1)
         end if
c
c     Form the residual in g_soln
c
         call ga_zero(g_soln)
         do i = 1, nsub
            call ga_dadd_patch(
     $           yy(i), g_p, 1, n, i, i, 
     $           one, g_soln, 1, n, 1, 1,
     $           g_soln, 1, n, 1, 1)
         end do
         call ga_dadd_patch(
     $         mone, g_b, 1, n, 1, 1, 
     $         one, g_soln, 1, n, 1, 1,
     $        g_soln, 1, n, 1, 1)
c
c     Print and determine convergence or too many iterations
c     
         rnorm = sqrt(ga_ddot(g_soln, g_soln))
         acc = rnorm/bnorm
         converged = acc .lt. tol
         if (iter.gt.1) then
            convg_ratio = rnorm/rnorm_prev
         else
            convg_ratio = zero
         end if
         rnorm_prev = rnorm
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter,nsub,rnorm,acc,convg_ratio,util_cpusec()
    3       format(14x,2i5,1p,2e12.4,0p,f10.4,f8.1)
            call util_flush(6)
         end if
         if (converged) goto 1000
         if (iter.eq.maxiter) then
            if (owarn) then
               write(6,*) ' ga_iter_lsolve: exceeded max iters ... ',
     $              'aborting solve'
               call util_flush(6)
            end if
            goto 1000
         end if
c
c     If things are converging very slowly punt back to 
c     the caller for he/she to decide what should happen
c
         if (ostall .and. convg_ratio .gt. 0.91d0) then
            if (owarn) then
               write(6,*) ' ga_iter_lsolve: convergence stagnant ... ',
     $              'aborting solve'
               call util_flush(6)
            end if
            goto 1000
         end if
c
c     Form a new update vector by preconditioning and orthogonalizing
c     approximately to the previous expansion space
c
*         call precond(g_soln, zero)
         call ga_precond_project(g_soln, zero, precond, g_w)
         do i = 1, nsub
            rx = ga_ddot_patch(
     $           g_soln, 'n', 1, n, 1, 1,
     $           g_x, 'n', 1, n, i, i)
            xx = ga_ddot_patch(
     $           g_x, 'n', 1, n, i, i, 
     $           g_x, 'n', 1, n, i, i)
            call ga_dadd_patch(
     $           one, g_soln, 1, n, 1, 1,
     $           -rx/xx, g_x, 1, n, i, i,
     $           g_soln, 1, n, 1, 1)
         end do
c
         if (nsub .lt. maxsub) then
            nsub = nsub + 1
         else
c
c     If have run out of space have to compress the list by summing
c     the oldest two vectors and shuffling down.
c
            call ga_dscal_patch(g_x, 1, n, 1, 1, yy(1))
            call ga_dscal_patch(g_p, 1, n, 1, 1, yy(1))
            call ga_dadd_patch(
     $           yy(2), g_x, 1, n, 2, 2,
     $           one, g_x, 1, n, 1, 1,
     $           g_x, 1, n, 1, 1)
            call ga_dadd_patch(
     $           yy(2), g_p, 1, n, 2, 2,
     $           one, g_p, 1, n, 1, 1,
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
            end do
         end if
      end do
c
c     Should never actually fall out of the bottom
c
      call errquit('ga_iter_lsolve: fell out of loop bottom?', 0,
     &       GA_ERR)
c
 1000 continue
      call ga_zero(g_soln)
      do i = 1, nsub
         call ga_dadd_patch(
     $        yy(i), g_x, 1, n, i, i, 
     $        one, g_soln, 1, n, 1, 1,
     $        g_soln, 1, n, 1, 1)
      end do
c
      ga_iter_lsolve = converged
c
      if (.not. ga_destroy(g_x)) call errquit
     $     ('ga_iter_lsolve: failed to destory x', 0, GA_ERR)
      if (.not. ga_destroy(g_p)) call errquit
     $     ('ga_iter_lsolve: failed to destory p', 0, GA_ERR)
      if (.not. ga_destroy(g_w)) call errquit
     $     ('ga_iter_lsolve: failed to destroy w', 0, GA_ERR)
c
      end
