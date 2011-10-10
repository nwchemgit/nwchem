      subroutine dft_roks_search_precond(rtdb, g_grad, g_work)
*
* $Id$
*
      implicit none
#include "crohf.fh"      
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
      integer max_precond_iter  ! Max. no. of iters for preconditioner
      integer max_precond_sub   ! Max. dim. of precond. iterative space
      parameter (max_precond_iter=20)
      double precision max_grad, gnorm ! Measures of gradient
      double precision min_shift ! Shift forcing 1-e precond to be
                                 ! positive
      double precision precond_acc ! Accuracy attained by solver
      integer precond_iter       ! No. of iterations used by solver
      double precision prodmax
      integer nvir
      double precision maxelt
      logical oprint_conv, oprint_eval, oprint_debug,
     $  oprint_solve
      logical ga_iter_lsolve
      external rohf_hessv_precond, dft_roks_hessv, ga_iter_lsolve
c
      if (oscfps) call pstat_on(ps_solve)
      oprint_conv = util_print('convergence', print_default)
      oprint_eval = util_print('intermediate evals', print_debug)
      oprint_debug= util_print('rohf_debug', print_debug)
      oprint_solve= util_print('solve', print_high)
c
      if (.not.rtdb_get(rtdb, 'scf:maxsub', MT_INT, 1, 
     $     max_precond_sub)) max_precond_sub = 10
c
      call ga_maxelt(g_grad, max_grad)
      gnorm = sqrt(ga_ddot(g_grad, g_grad))
c
 30   continue
c
c     Determine shift necessary to make the preconditioner
c     positive definite
c
      call rohf_make_shift(min_shift)
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
         call rohf_canon(.false., oprint_eval)
         if (oadapt) call scf_movecs_sym_adapt(basis,g_movecs,.false.,
     $       0, 'after canon',.true.,int_mb(k_irs))
         call rohf_get_grad(g_grad)
         call rohf_make_shift(min_shift) ! Recompute
      endif
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
         endif
      endif
c
      call ga_zero(g_work)      ! ESSENTIAL
c     
c     Attempt to solve the damned equations
c
 20   if (.not. ga_iter_lsolve(crohf_vlen, max_precond_iter,
     $     max_precond_sub, precond_tol, 
     $     rohf_hessv_precond, dft_roks_hessv, .true., oprint_solve, 
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
c     Screen out symmetry junk
c
      nvir = nmo - nclosed - nopen
      if(oadapt) call scf_sym_screen(nclosed, nopen, nvir,int_mb(k_irs),
     $  g_work, 0, .true., maxelt)
c
c     Screen small rotations to keep out noise.
c     Cannot screen with a large number (greater than accuracy in
c     integrals) since in a non-abelian group we could break symmetry...
c     either screen with a small threshold or screen with a large one
c     and then enforce symmetry (which is expensive in high order
c     groups)
c     Also, tol must tend to zero for tight convergence.
c
      call ga_maxelt(g_work, prodmax)
      prodmax = min(prodmax**2,prodmax*1d-2,1d-12)
      call ga_screen(g_work, prodmax)
      if (oscfps) call pstat_off(ps_solve)
c
      end
