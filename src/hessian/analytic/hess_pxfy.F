C
C $Id: hess_pxfy.F,v 1.3 2001-12-28 18:30:08 windus Exp $
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
c              write(6,*) 'deltah for i,j ',deltah,ipert,jpert
c              if (restr) then
c                 deltah = deltah * 2.d0
c              endif
               if (ipert .ne. jpert) then
                  deltah = deltah * 0.5d+00
                  hess(ipert,jpert) = hess(ipert,jpert) - deltah
                  hess(jpert,ipert) = hess(jpert,ipert) - deltah
               else
                  hess(ipert,jpert) = hess(ipert,jpert) - deltah
               endif
            enddo
         endif
      enddo
c
      return
      end
C
      subroutine hess_pxfy2(hess, g_sol, g_rhs, oactive, ncent, ipert,
     &     restr, nocc, nvirt, ndens)
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
      integer g_rhs             ! [input] GA handles for RHS's
      logical oactive(ncent)    ! [input] key to active atoms
      integer ipert             ! [input] linearized index of current pert
      logical restr             ! [input] calculation is restricted
      integer nocc              ! [input] number of occupieds
      integer nvirt             ! [input] number of virtuals
      integer ndens             ! [input] spin of current solution
c
      integer jcent, jdir, jpert
      integer alo(2), ahi(2), blo(3),bhi(3)
      double precision deltah
c
c     Loop over derivative Fock matrices and contract with current
c     CPSCF solution
c
      alo(1) = 1
      ahi(1) = nocc
      alo(2) = 1
      ahi(2) = nvirt
      blo(2) = 1
      bhi(2) = nocc
      blo(3) = 1
      bhi(3) = nvirt
      do jcent = 1, ncent
         if (oactive(jcent)) then
            jpert = 3*(jcent-1)
            do jdir = 1, 3
               jpert = jpert + 1
               blo(1) = jpert + (ndens-1)*ncent*3
               bhi(1) = jpert + (ndens-1)*ncent*3
               deltah = nga_ddot_patch(g_sol, 'n', alo, ahi,
     $                                 g_rhs, 'n', blo, bhi)
c              write(6,*) 'deltah for i,j ',deltah,ipert,jpert
c              if (restr) then
c                 deltah = deltah * 2.d0
c              endif
               if (ipert .ne. jpert) then
                  deltah = deltah * 0.5d+00
                  hess(ipert,jpert) = hess(ipert,jpert) - deltah
                  hess(jpert,ipert) = hess(jpert,ipert) - deltah
               else
                  hess(ipert,jpert) = hess(ipert,jpert) - deltah
               endif
            enddo
         endif
      enddo
c
      return
      end
