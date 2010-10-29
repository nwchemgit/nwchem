      logical function line_search(nstep, max_nstep, energies,
     $     gradients, steps, eacc, tol, oprint)
C$Id$
      implicit none
#include "errquit.fh"
#include "util.fh"
c
      integer nstep             ! Current no. of steps ... initially 1
      integer max_nstep         ! Max. no. of steps
      double precision energies(max_nstep) ! Energies of points in LS
      double precision gradients(max_nstep) ! Gradients of points in LS
      double precision steps(max_nstep) ! Steps of points in LS
      double precision eacc     ! Absolute precision in energies
      double precision tol      ! Conv if tol > |grad_final|/|grad_init|
      logical oprint            ! Printing enabled if true
c
      double precision gs0, gs1, e0, e1, step0, step1, ediff, ls_ratio
      double precision ls_hessian, ls_grad, ls_step, ls_step_max
      double precision e_low, step_low, g_low
      double precision step_neg_h
      integer i, nstep_neg_h
c
c     Line search with gradients.  Current version only uses the
c     last two points in the LS but subsequent versions can use
c     all of the available data.
c
c     Also assumes good preconditioning so that first step is always
c     unity.
c
      if (nstep.lt.1) then
         call errquit('line_search: nstep ivalid', nstep, INPUT_ERR)
      else if (nstep.eq.1) then
         if (oprint) then
            write(6,1)
 1          format(/,30x,'----- Line search ------'/)
            write(6,2)
 2          format(1x,
     $           '       step       energy         lsg     num-lsg ',
     $           '  num-lsh   lsrat    time'/1x,
     $           '     ------- ---------------- --------- ---------',
     $           ' --------- ------- --------')
            write(6,3) steps(1), energies(1), gradients(1), 
     $           0.0, 0.0, 0.0, util_cpusec()
            call util_flush(6)
         endif
         if (steps(1).ne.0.0d0) call errquit
     $        ('line_search: step(1) .ne. 0.0?', 0, INPUT_ERR)
         nstep = nstep + 1
         steps(nstep) = 1.0d0
         line_search = .false.
         return
      endif
c     
c     determine current lowest energy and step
c     
      e_low    = energies(1)
      step_low = steps(1)
      g_low    = gradients(1)
      do i = 2, nstep
         if (energies(i).lt.e_low) then
            e_low    = energies(i)
            step_low = steps(i)
            g_low    = gradients(i)
         endif
      enddo
c     
c     numerical estimates of hessian and gradient
c     
      gs1 = gradients(nstep)
      gs0 = gradients(nstep-1)
      step1 = steps(nstep)
      step0 = steps(nstep-1)
      e1 = energies(nstep)
      e0 = energies(nstep-1)
      ediff = e1 - e0
      ls_hessian = (gs1 - gs0) / (step1 - step0)
      ls_grad = (e1 - e0) / (step1 - step0) +
     $     0.5d0*ls_hessian*(step1 - step0)
      ls_ratio = abs(gs1 / gradients(1))
c     
      if (oprint) then
         write(6,3) step1, e1, gs1,
     $        ls_grad, ls_hessian, ls_ratio, util_cpusec()
 3       format(5x,f8.4,f17.8,1p,3d10.2,d8.1,0p,f9.1)
         call util_flush(6)
      endif
c     
c     Test for convergence of the line search
c     
      if ( (ls_ratio .lt. tol) .and.
     $     (ediff .lt. eacc) .and.
     $     (e1-e_low) .lt. eacc) then
         line_search = .true.
         return
      endif
c     
c     Now determine where we are going next
c     
      if (((e1-e_low).gt.eacc)) then
*     .and. abs(step1-step_low).gt.0.1d0) then
c     
c     Have gone up-hill ... backtrack optimistically making sure that
c     we are on the downhill side of the lowest energy point
c     
         if (oprint) then
            write(6, 4) e1 - e_low
            call util_flush(6)
 4          format(20x,'Energy increase :',f23.8,' backtrack')
         endif
c
c     The next factor was 0.5 but this can lead to a cycle being
c     formed by constantly backtracking and then restricting 
c     the step (it also has a factor of 0.5).  Hence, 1/3
c     which is relatively prime.
c
         ls_step = (step_low - step1) / 3d0
         if ((step1 - step_low + ls_step)*g_low .gt. 0) then
            if (oprint) then
               write(6, 41) g_low
               call util_flush(6)
 41            format(20x,'On wrong side of minimum energy point ',f6.2)
            endif
            ls_step = step_low - step1 - 0.1d0*g_low
         endif
      else
c     
c     Newton Raphson predicted step
c     
         ls_step = (gs1*step0 - gs0*step1)/(gs1 - gs0) - step1
c     
c     if the hessian is positive and the numerical gradient agrees
c     well with the computed gradient (indicating that anharmonic
c     terms may be small) then can take a bigger step.
c     
         ls_step_max = min(4.0d0,max(1.0d0, 0.5d0*step1))
         if (abs((ls_grad-gs1)/gs1) .lt. 0.01) then
            ls_step_max = 4.0d0*ls_step_max
         endif
c     
         if (ls_hessian .lt. 0.0d0) then
c     
c     If have negative curvature simply step a small amount forward
c     following the gradient downhill.  If the previous step also
c     had negative curvature, then double the step size.
c     
            if (oprint) then
               write(6,912) ls_hessian
               call util_flush(6)
 912           format(20x,' Negative curvature: hessian=',1p,d9.2)
            endif
            step_neg_h = 0.501d0  ! 0.501 identifies neg H step
            ls_step = sign(step_neg_h,-gs1)
            nstep_neg_h = nint(abs(step0-step1)/step_neg_h)
            if (abs(step_neg_h*nstep_neg_h-abs(step0-step1)).lt.1d-10)
     $          ls_step = ls_step*nstep_neg_h*2
c
         else if ((ls_step+step1) .lt. 0.0) then
c     
c     overall step cannot be negative since we have followed the gradient
c     
            if (oprint) then
               write(6,*) ' restricting step sign ', ls_step
               call util_flush(6)
            endif
            ls_step = -0.1d0*step1
         else if (abs(ls_step).gt.ls_step_max) then
c     
c     Enforce maximum step
c     
            if (oprint) then
               write(6,*) ' restricting step size ', ls_step
               call util_flush(6)
            endif
            ls_step = sign(ls_step_max,ls_step)
         endif
      endif
c     
      if (nstep .lt. max_nstep) then
         nstep = nstep + 1
         steps(nstep) = steps(nstep-1) + ls_step
      else
         if (oprint) then
            write(6,*) ' !! Warning ... line search did not converge !!'
            call util_flush(6)
         endif
         line_search = .true.
         return
*         call errquit('line_search: exhausted max steps',max_nstep)
      endif
      line_search = .false.
c
      end
      
