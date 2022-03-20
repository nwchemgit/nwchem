      subroutine gw_cdgw_buildin_graph(iW,Evals,x,w,iWidx,omega,output,
     $                            eta,ngl,myld,nqp)
C**********************************************************************
C*                                                                    *
C*    Computes the contribution from the numerical integral over the  *
C*    imaginary axes to the self-energy matrix elements               *
C*                                                                    *
C*    This subroutine uses the screened Coulomb matrix elements       *
C*    stored in iW (corresponding to the nominator of Equation (40)), *
C*    in order to compute I_{nn} in Equation (40) of                  *
C*                                                                    *
C*    JCTC 17, 7504 (2021)                                            *
C*                                                                    *
C*    and the corresponding derivative with respect to omega.         *
C*                                                                    *
C*                                                                    *
C**********************************************************************
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: ngl, myld, nqp
      integer,intent(in),dimension(0:myld,nqp) :: iWidx
      double precision, intent(in) :: omega,eta
      double precision, intent(in) :: Evals(*)
      double precision, intent(in) :: x(*), w(*)
      double precision, intent(in), dimension(myld,ngl,nqp) :: iW
      double precision, intent(out),dimension(2,nqp) :: output

      integer igl,jmo,mynmo, iqp
      double precision w2
      double precision,dimension(myld) :: vector
      double complex,dimension(2,nqp) :: res
      double complex,dimension(myld) :: etavec,factor,dfactor
      double complex :: imag = dcmplx(0d0,1d0)
      double complex :: ieta,iomega,afactor
      
      ieta = dcmplx(0d0,eta)
      res(:,:) = 0d0

      !
      ! Loop over all requested quasiparticle indices
      !
      do iqp=1,nqp

        mynmo = iwidx(0,iqp)
        if (mynmo.eq.0) cycle

        !
        ! Precompute arrays appearing in the denominator of 
        ! Equation (40)
        do jmo=1,mynmo
          vector(jmo) = omega - Evals(iwidx(jmo,iqp))
          etavec(jmo) = sign(eta,evals(iwidx(jmo,iqp)))
        enddo

        !
        ! Loop over all imaginary grid points. The values of all I_{nn}
        ! will be stored in res(:,1), while the derivatives with respect
        ! to omega will be storedin res(:,2)
        !
        do igl=1,ngl
          iomega = imag*dsqrt(x(igl))
          do jmo=1,mynmo
            factor(jmo) = 1d0/(vector(jmo) + iomega + etavec(jmo)*0d0)
            dfactor(jmo) = factor(jmo)**2
          enddo
          res(1,iqp) = res(1,iqp) - w(igl)*
     $                 dot_product(iw(1:mynmo,igl,iqp),factor(1:mynmo))
          res(2,iqp) = res(2,iqp) + w(igl)*
     $                dot_product(iw(1:mynmo,igl,iqp),dfactor(1:mynmo))
        enddo
      enddo

      !
      ! Sum contributions from all MPI ranks.
      !
      call ga_zgop(1038,res,2*nqp,'+')

      !
      ! Copy the results to the appropriate arrays
      !
#ifdef GWCMPLX
      output(:,:) = res(:,:)
#else
      output(:,:) = dble(res(:,:))
#endif

      end subroutine
