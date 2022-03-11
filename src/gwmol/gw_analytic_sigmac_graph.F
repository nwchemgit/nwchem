      subroutine gw_analytic_sigmac_graph(sigma,omegain,Eold,Omegam,
     $                              wmn,eta,npoles,nmo,plus,me,npoints)
c
c     Obtains the diagonal matrix elements of the self-energy operator
c     (Equation (37)). This version computes many omegas at once.
c
      implicit none
#include "gw.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      integer,intent(in) :: wmn,npoles,nmo,plus,me,npoints
      double precision,intent(in) :: omegain(npoints),Eold(nmo)
      double precision,intent(in) :: eta,Omegam(npoles)
      double precision,intent(out) :: sigma(npoints)

      integer :: ilo, ihi, jlo, jhi, ipole, imo, ld, adrWMN
      double precision :: omega(npoints),w,delta,denom, eta3

      call ga_distribution(wmn, me, ilo, ihi, jlo , jhi)

      sigma(:) = 0d0

      call ga_access(wmn, ilo, ihi ,jlo, jhi, adrWMN, ld)
      eta3 = 3d0*eta

      !
      ! Perform the summation over the occupied orbitals
      !
!$omp parallel private(omega,w)
!$omp do collapse(2) reduction(+: sigma )
      do imo=1,plus
        do ipole=jlo,jhi
          omega(:) = omegain(:) - Eold(imo) + omegam(ipole)
          w = dbl_mb(adrWMN + (ipole-jlo)*ld + imo - 1)
          sigma(:) = sigma(:) + w*omega(:)/(omega(:)**2+eta)
        enddo
      enddo
!$omp end do      
c
      !
      ! Perform the summation over the unoccupied orbitals
      !
!$omp do collapse(2) reduction(+: sigma )
      do imo=plus+1,nmo
        do ipole=jlo,jhi
          omega(:) = omegain(:) - Eold(imo) - omegam(ipole)
          w = dbl_mb(adrWMN + (ipole-jlo)*ld + imo - 1)
          sigma(:) = sigma(:) + w*omega(:)/(omega(:)**2+eta)
        enddo
      enddo
!$omp end  do      
!$omp end parallel
c
      call ga_release(wmn, ilo, ihi ,jlo, jhi)

      !
      ! Sum contributions from all MPI ranks.
      !
      call ga_dgop(1038,sigma,npoints,'+')

      end subroutine
