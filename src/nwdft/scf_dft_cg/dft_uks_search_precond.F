      subroutine dft_uks_search_precond(rtdb, g_grad, g_work)
      implicit none
#include "cuhf.fh"      
#include "cscfps.fh"
#include "global.fh"
#include "cscf.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      integer g_grad
      integer g_work
c
c     Precondition the gradient with an approximation to the Hessian
c     return the result in g_work.
c
c     Old version only divided by the level-shifted diagonals which 
c     is in rohf the preconditioner to the preconditioner!
c
      logical oprint_conv, oprint_eval, oprint_debug, oprint_solve
      double precision gnorm
      double precision max_grad
      double precision min_shift
      integer max_precond_iter  ! Max. no. of iters for preconditioner
      integer max_precond_sub   ! Max. dim. of precond. iterative space
      parameter (max_precond_iter=20)
c
      double precision precond_acc ! Accuracy attained by solver
      integer precond_iter       ! No. of iterations used by solver
c
      integer nvir, ivoff
      double precision maxelt
c
      external uhf_precond, dft_uks_hessv
      logical ga_iter_lsolve
      external ga_iter_lsolve
c
      oprint_conv = util_print('convergence', print_default)
      oprint_eval = util_print('intermediate evals', print_debug)
      oprint_debug= util_print('uhf_debug', print_debug)
      oprint_solve= util_print('solve', print_high)
c
      if (.not.rtdb_get(rtdb, 'scf:maxsub', MT_INT, 1, 
     $     max_precond_sub)) max_precond_sub = 10
c
      call ga_maxelt(g_grad, max_grad)
      gnorm = sqrt(ga_ddot(g_grad, g_grad))
      eprec = max(1d-7,gnorm*0.01d0, tol2e*100.0d0) ! Approx precision
                                                    ! available
c
c     Determine shift necessary to make the preconditioner
c     positive definite
c
 30   continue
      call uhf_make_shift(min_shift)
c
c     Figure out which preconditioner to use, convergence
c     thresholds, shifts, ...
c
      call scf_precond_select(gnorm, max_grad, min_shift,
     $     oprint_conv)
c
c     Do a diagonalization if it will make the Hessian more
c     diagonally dominant
c
      if (ododiag) then
         call uhf_canon(.false., oprint_eval)
         if (oadapt) then
           call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     $           0, 'after canon', .true., int_mb(k_irs))
           call scf_movecs_sym_adapt(basis, g_movecs(2), .false.,
     $           0, 'after canon', .true., int_mb(k_irs+nmo))
         endif
         call uhf_get_grad(g_grad)
         call uhf_make_shift(min_shift) ! Recompute
      end if
c
c     After doing the diag make sure that the shift is still OK
c
      if (lshift .lt. min_shift) then
         lshift = min_shift + 2.0d0
         if (ga_nodeid().eq.0 .and. oprint_conv) then
            write(6,3131) lshift
 3131       format('  Setting level-shift to ', f6.2,
     $           ' to force positive preconditioner')
            call util_flush(6)
         end if
      end if
c
c     Next two lines precondition with just the eigenvalue differences
c
*      call ga_copy(g_grad,g_work)
*      call uhf_precond(g_work, 0.0d0)
c
      call ga_zero(g_work)      ! ESSENTIAL
c     
c     Attempt to solve the damned equations
c
 20   if (.not. ga_iter_lsolve(cuhf_vlen, max_precond_iter,
     $     max_precond_sub, precond_tol, 
     $     uhf_precond, dft_uks_hessv, .true., oprint_solve, 
     $     g_grad, g_work, precond_acc, precond_iter)) then
c     
c        Iterative solution did not converge.  If we got at least
c        1 sig. fig. then just continue.  Otherwise for the 1-e
c        approximation increase the level-shift and restart with
c        the current vector as guess.  If the 2-e approximation did
c        not converge disable it forever and resort to the 1-e method.
c     
         if (precond_acc .gt. 0.1d0) then
            if (pflg .eq. 2) then
               odisable_nr = .true.
               maxiter = maxiter + 10
               if (ga_nodeid() .eq. 0 .and. oprint_conv) then
                  write(6,22) maxiter
 22               format(/' Disabled NR: increased maxiter to ',i3/)
                  call util_flush(6)
               endif
               goto 30
            else
               if (lshift .eq. 0.0d0) then
                  lshift = 2.0d0
                  call ga_zero(g_work)
               else
                  lshift = lshift*2.0d0
               endif
            endif
            if (ga_nodeid() .eq. 0 .and. oprint_conv) then
               write(6,2) lshift
 2             format(' Increased level shift to ', f8.2)
               call util_flush(6)
            endif
            goto 20
         endif
      endif
c
      if(oadapt) then
        nvir = nmo - nalpha
        call scf_sym_screen(nalpha, 0, nvir, int_mb(k_irs),
     $       g_work, 0, .true., maxelt)
        nvir = nmo - nbeta
        ivoff = nalpha*(nmo-nalpha)
        call scf_sym_screen(nbeta, 0, nvir, int_mb(k_irs+nmo),
     $       g_work, ivoff, .true., maxelt)
      endif
c
      end
c $Id$
