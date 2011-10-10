      subroutine uccsd(rtdb)
      implicit none
      integer rtdb
c
c     initial guess for amplitudes
c     zero for now.

c     open product file

c     initialize solver

      do iter = 1, maxiter

c     lambda
         call uccsdt_lambda(d_amp)

c     bterm

c     cterm

c     aterm and rest

c     ENERGY WHERE?

c     solver

c     update and save for restart

      end do

c     analysis

      end
c $Id$
