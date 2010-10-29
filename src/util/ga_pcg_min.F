      subroutine ga_pcg_minim(n, iter, e, g_grad, g_work, g_s, step,
     $     ls_tol, ls_max, eacc, conjugacy, oprint, oprint_ls,
     $     iwork, dwork, mode)
C$Id$
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer n                 ! No. of parameters [input]
      integer iter              ! Current macro iteration no. [input/output]
      double precision e        ! Energy of current point [input]
      integer g_grad            ! GA with current gradient [input]
      integer g_work            ! GA for work space [input]
      integer g_s               ! GA with current search direction [output]
      double precision step     ! Step to take along the line search direction [output]
      double precision ls_tol   ! Tolerance for current line search [input]
      integer ls_max            ! Max. no. of steps for line search [input]
      double precision eacc     ! Absolution precision in energy [input]
      logical conjugacy         ! True if conjugacy is employed [input]
      logical oprint            ! True if should print convergence [input]
      logical oprint_ls         ! True if should print line search [input]
      integer iwork(1000)       ! Integer work space [input]
      double precision dwork(3*ls_max) ! Double precision work space [input]
      character*16 mode         ! Used to communicate with calling routine
c
c     Use preconditioned conjugate gradient to minimize the value of
c     a function.  It is assumed that the preconditioning is pretty
c     good so that steps are generally in the vicinity of unity
c     (this assumption is from the line-search).
c
c     Have a look in testsolve.F (routine testpcg) for sample usage.
c
c     g_grad holds the current gradient
c     g_s    holds the current search direction
c     g_work holds the previous preconditioned gradient
c     
      integer action, act_startup, act_ls, act_precond, act_conj,
     $     act_accept
      parameter (act_startup=0, act_precond=1, act_conj=2, act_ls=3,
     $     act_accept=4)
c     
c     Pointers into iwork and dwork for info
c     
      integer iwork_action      ! Action
      integer iwork_ls_nstep    ! Current no. of line search steps
      parameter (iwork_action = 1, iwork_ls_nstep = 2)
c     
      integer ls_e, ls_g, ls_s, sk_gk, sk_gk1, gk1_c_gk, gk1_c_gk1
      parameter (sk_gk=1, sk_gk1=2, gk1_c_gk=3, gk1_c_gk1 = 4, ls_e = 5)
c     
      logical ls_print          ! oprint_ls restricted to node 0
      integer nstep
      logical line_search
      double precision beta, numerator, denominator, gnorm, gmax
      double precision mone
      parameter (mone = -1.0d0)
      external line_search
c     
      if (mode .eq. 'startup') iwork(iwork_action) = act_startup
c     
c     Set up pointers and look inside iwork to see what we are
c     currently doing
c     
      action = iwork(iwork_action)
      nstep  = iwork(iwork_ls_nstep)
      ls_g = ls_e + ls_max
      ls_s = ls_g + ls_max
      ls_print = .false.
c     
      if (action .eq. act_startup) then
c     
c     Startup calculation ... get energy and gradient and the
c     next thing needed is a preconditioning
c     
         iter = 1
         step = 0.0d0
         call dfill((3*ls_max), 0.0d0, dwork, 1)
         call ifill(2, 0, iwork, 1)
         call ga_zero(g_work)
         call ga_zero(g_s)
         mode = 'energy+gradient'
         action = act_accept
      else if (action .eq. act_precond) then
c     
c     Preconditioning required.  Next thing to do is the conjugacy
c     
         dwork(sk_gk1) = ga_ddot(g_s, g_grad)
         dwork(gk1_c_gk) = ga_ddot(g_grad, g_work)
         mode = 'precondition'
         action = act_conj
      else if (action .eq. act_conj) then
c     
c     Apply the conjugacy condition and then commence the line search
c     
         dwork(gk1_c_gk1) = ga_ddot(g_grad, g_work)
         if (.not. conjugacy) then
            beta = 0.0d0
         else
            numerator = (dwork(gk1_c_gk1) - dwork(gk1_c_gk))
            denominator = (dwork(sk_gk1)  - dwork(sk_gk))
            if (denominator .ne. 0.0d0) then
               beta = numerator/denominator
            else
               beta = 0.0d0
            endif
         endif
         if (beta .lt. 0.0d0) then
*            if (ga_nodeid() .eq. 0) then
*               write(6,*) ' pcg: negative beta = ', beta
*               call util_flush(6)
*            endif
            beta = 0.0d0
         endif
         call ga_dadd(mone, g_work, beta, g_s, g_s)
         dwork(sk_gk) = ga_ddot(g_s, g_grad) ! For the next iteration
c     
         step = 0.0d0
         nstep = 1
         mode = 'energy+gradient'
         action = act_ls
      else if (action .eq. act_accept) then
c     
c     Accept the current point ... next will be a precond
c     
         mode = 'accept step'
         action = act_precond
      endif
c     
      if (action .eq. act_ls) then
c     
c     We have a new energy+gradient for a line search underway.  Determine
c     if the LS has converged or what the next step should be
c     
         dwork(ls_e+nstep-1) = e
         dwork(ls_g+nstep-1) = ga_ddot(g_grad, g_s)
         dwork(ls_s+nstep-1) = step
         ls_print = oprint_ls .and. (ga_nodeid().eq.0)
         if (line_search(nstep, ls_max, dwork(ls_e), dwork(ls_g),
     $        dwork(ls_s), eacc, ls_tol, ls_print)) then
c     
c     The LS has converged.  Return the current step as accepted.
c     Next action will be to precondition the gradient.
c     
            mode = 'accept step'
            action = act_precond
         else
c     
c     Just continue with the line search
c     
            step = dwork(ls_s+nstep-1)
         endif
      endif
c
      if (mode .eq. 'accept step') then
         if (oprint) then
            gnorm = sqrt(ga_ddot(g_grad, g_grad))
            call ga_maxelt(g_grad, gmax)
            if (ga_nodeid() .eq. 0) then
c     
c     Don't write out the header unless are also printing LS info
c     except on the first iteration ... keeps output compact
c
               if (ls_print .or. iter.eq.1) write(6,1)
 1             format(/,
     $              13x,' iter       energy          gnorm     gmax   ',
     $              '    time'/
     $              13x,'----- ------------------- --------- ---------',
     $              ' --------')
               write(6,2) iter, e, gnorm, gmax, util_cpusec()
 2             format(13x,i5,f20.10,1p,2d10.2,0p,f9.1)
               call util_flush(6)
            endif
         endif
         iter = iter + 1
      endif
c
      iwork(iwork_action) = action
      iwork(iwork_ls_nstep) = nstep
c     
      end

