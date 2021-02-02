      subroutine gw_analytic_sigmac(sigma,dsigma,d2sigma,omegain,Eold,
     $                              Omegam,wmn,scratch,eta,npoles,nmo,
     $                              plus,me)
      implicit none
#include "gw.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      integer,intent(in) :: wmn,npoles,nmo,plus,me
      double precision,intent(in) :: omegain,Eold(nmo)
      double precision,intent(in) :: eta,Omegam(npoles)
      double precision,intent(inout) :: scratch(nmo,*)
      double precision,intent(out) :: sigma,dsigma,d2sigma

      integer :: ilo, ihi, jlo, jhi, ipole
      double precision :: res(3)

      call ga_distribution(wmn, me, ilo, ihi, jlo , jhi)

      res(:) = 0d0
      do ipole=jlo,jhi
        scratch(1:plus,1) = omegain - Eold(1:plus) + omegam(ipole)
        scratch(plus+1:nmo,1) = omegain - Eold(plus+1:nmo) - 
     $                          omegam(ipole)

        scratch(:,2) = scratch(:,1)/(scratch(:,1)**2 + eta)
        scratch(:,3) = (eta-scratch(:,1)**2)/(scratch(:,1)**2+eta)**2
        scratch(:,4) = 2d0*scratch(:,1)*(scratch(:,1)**2 - 3d0*eta)/
     $                 (scratch(:,1)**2 + eta)**3

        call ga_get(wmn,1,nmo,ipole,ipole,scratch,nmo)
        call dgemv('t',nmo,3,1d0,scratch(1,2),nmo,scratch,1,1d0,res,1)
      enddo

      call dgop(mt_dbl,res,3,'+')
      sigma = res(1)
      dsigma = res(2)
      d2sigma = res(3)

      end subroutine
