C
C $Id: hess_pxfy.F,v 1.1 2000-05-11 11:16:08 bjohnson Exp $
C
      subroutine hess_pxfy(hess, g_sol, g_rhs, oactive, ncent, ipert,
     &     restr)
c
c     Accumulates dP/dx * dF/dy contribution to hessian, where dF/dy is
c     the explicit contribution to the derivative of the Fock matrix, i.e.
c     the RHS of the CPSCF.
c
c     Set up for current approach of doing one CPSCF solution at a time.
c     The current solution is contracted with all the RHS's.
c
c     BGJ (5/00)
c
      implicit none
c
#include "global.fh"
c
      integer ncent             ! [input] number of nuclear centers
      double precision hess(3*ncent,3*ncent) ! [updated] hessian matrix
      integer g_sol             ! [input] GA handle for current solution
      integer g_rhs(3*ncent)    ! [input] GA handles for RHS's of same
c                               !         spin as current solution
      logical oactive(ncent)    ! [input] key to active atoms
      integer ipert             ! [input] linearized index of current pert
      logical restr             ! [input] calculation is restricted
c
      integer jcent, jdir, jpert
      double precision deltah
c
c     Loop over derivative Fock matrices and contract with current
c     CPSCF solution
c
      do jcent = 1, ncent
         if (oactive(jcent)) then
            jpert = 3*(jcent-1)
            do jdir = 1, 3
               jpert = jpert + 1
               deltah = ga_ddot(g_sol, g_rhs(jpert))
               if (restr) then
                  deltah = deltah * 2.d0
               endif
               hess(ipert,jpert) = hess(ipert,jpert) + deltah
               if (ipert .ne. jpert) then
                  hess(jpert,ipert) = hess(jpert,ipert) + deltah
               endif
            enddo
         endif
      enddo
c
      return
      end
