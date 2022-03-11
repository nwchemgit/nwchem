      subroutine gw_analytic_sigmac(sigma,dsigma,d2sigma,omegain,Eold,
     $                              Omegam,wmn,eta,npoles,nmo,
     $                              plus,me)
c
c     Obtains the diagonal matrix elements of the self-energy operator
c     (Equation (37)). 
c
      implicit none
#include "gw.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      integer,intent(in) :: wmn,npoles,nmo,plus,me
      double precision,intent(in) :: omegain,Eold(nmo)
      double precision,intent(in) :: eta,Omegam(npoles)
      double precision,intent(out) :: sigma,dsigma,d2sigma

      integer :: ilo, ihi, jlo, jhi, ipole, imo, ld, adrWMN
      double precision :: omega,w,delta,denom, eta3
      double precision  :: resp(3),resm(3)

      call ga_distribution(wmn, me, ilo, ihi, jlo , jhi)

      resp(:) = 0d0
      resm(:) = 0d0

      call ga_access(wmn, ilo, ihi ,jlo, jhi, adrWMN, ld)
      eta3 = 3d0*eta

      !
      ! Perform the summation over the occupied orbitals
      !
!$omp parallel private(delta,omega,w,denom)
!$omp do collapse(2) reduction(+: resp )
      do imo=1,plus
        do ipole=jlo,jhi
          omega = omegain - Eold(imo) + omegam(ipole)
          w = dbl_mb(adrWMN + (ipole-jlo)*ld + imo - 1)
          denom = omega**2 + eta
          resp(1) = resp(1) + w*omega/denom
          resp(2) = resp(2) + w*(eta-omega**2)/denom**2
          resp(3) = resp(3) + w*2d0*omega*(omega**2-eta3)/denom**3
        enddo
      enddo
!$omp end do      
c
      !
      ! Perform the summation over the unoccupied orbitals
      !
!$omp do collapse(2) reduction(+: resm )
      do imo=plus+1,nmo
        do ipole=jlo,jhi
          omega = omegain - Eold(imo) - omegam(ipole)
          w = dbl_mb(adrWMN + (ipole-jlo)*ld + imo - 1)
          denom = omega**2 + eta
          resm(1) = resm(1) + w*omega/denom
          resm(2) = resm(2) + w*(eta-omega**2)/denom**2
          resm(3) = resm(3) + w*2d0*omega*(omega**2-eta3)/denom**3
        enddo
      enddo
!$omp end  do      
!$omp end parallel
c
      resp = resp + resm

      !
      ! Sum contributions from all MPI ranks
      !
      call ga_dgop(1038,resp,3,'+')
      sigma = resp(1)
      dsigma = resp(2)
      d2sigma = resp(3)

      call ga_release(wmn, ilo, ihi ,jlo, jhi)

      end subroutine
