c $Id: hess_pxfy.F,v 1.5 2001-12-28 21:54:02 windus Exp $
      subroutine hess_pxfy(hess, g_sol, g_rhs, oactive, ncent,
     &     restr, nocc, nvirt, ndens)
c
c     Accumulates dP/dx * dF/dy contribution to hessian, where dF/dy is
c     the explicit contribution to the derivative of the Fock matrix, i.e.
c     the RHS of the CPSCF.
c
c     Set up for current approach of doing all CPSCF solutions at a time.
c     The solutions are contracted with all the RHS's.
c
      implicit none
c
#include "global.fh"
c
      integer ncent             ! [input] number of nuclear centers
      double precision hess(3*ncent,3*ncent) ! [updated] hessian matrix
      integer g_sol             ! [input] GA handle for solutions
      integer g_rhs             ! [input] GA handles for RHS's
      logical oactive(ncent)    ! [input] key to active atoms
      logical restr             ! [input] calculation is restricted
      integer ndens             ! [input] number of densities
      integer nocc(ndens)       ! [input] number of occupieds
      integer nvirt(ndens)      ! [input] number of virtuals
c
      integer icent, jcent, idir, jdir, ipert, jpert
      integer index, jindex, ndns
      integer alo(2), ahi(2), blo(3),bhi(3)
      double precision deltah
c
c     Loop over CPSCF RHS and contract with current
c     CPSCF solution
c
      do ndns = 1, ndens
       alo(1) = 1
       ahi(1) = nocc(ndns)*nvirt(ndns)
       blo(1) = 1
       bhi(1) = nocc(ndns)*nvirt(ndns)
       do icent = 1, ncent
        if (oactive(icent)) then
         ipert = 3*(icent-1)
         do idir = 1, 3
          ipert = ipert + 1
          index = ipert + (ndns-1)*ncent*3
          alo(2) = index
          ahi(2) = index
          do jcent = 1, icent
           if (oactive(jcent)) then
            jpert = 3*(jcent-1)
            do jdir = 1, 3
              jpert = jpert + 1
              if (jpert.le.ipert) then ! take advantage of symmetry of hessian
               jindex = jpert + (ndns-1)*ncent*3
               blo(2) = jindex
               bhi(2) = jindex
               deltah = nga_ddot_patch(g_sol, 'n', alo, ahi,
     $                                 g_rhs, 'n', blo, bhi)
               write(6,*) 'deltah for i,j ',deltah,ipert,jpert
c              if (restr) then
c                 deltah = deltah * 2.d0
c              endif
               if (ipert .ne. jpert) then
                  hess(ipert,jpert) = hess(ipert,jpert) - deltah
                  hess(jpert,ipert) = hess(jpert,ipert) - deltah
               else
                  hess(ipert,jpert) = hess(ipert,jpert) - deltah
               endif
              endif
            enddo
           endif
          enddo
         enddo
        endif
       enddo
      enddo
c
      return
      end
