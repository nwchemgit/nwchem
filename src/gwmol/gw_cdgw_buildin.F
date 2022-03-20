#ifdef GWCMPLX
#define GWTYP double complex
#else
#define GWTYP double precision
#endif

      subroutine gw_cdgw_buildin(iW,Evals,x,w,iWidx,omega,eta,I,dI,ngl,
     $                           myld)
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
C*    and the corresponding derivative with respect to omega          *
C*                                                                    *
C**********************************************************************

      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: ngl, myld
      integer,intent(in),dimension(0:myld) :: iWidx
      double precision, intent(in) :: omega
      double precision, intent(in) :: Evals(*)
      double precision, intent(in) :: x(*), w(*)
      double precision, intent(in) :: eta
      GWTYP, intent(in), dimension(myld,ngl) :: iW
      GWTYP, intent(out) :: I,dI

      integer igl,jmo,mynmo
      double complex,dimension(iwidx(0)) :: etavec
      double precision          :: vector(iwidx(0))

      double complex             :: iomega, factor, res(2), tmp(2)
      double complex, parameter  :: imag = dcmplx(0d0,1d0)
      double complex  :: ieta

      ieta = dcmplx(0d0,eta)

      mynmo = iwidx(0)
      res(:) = dcmplx(0d0,0d0)
      tmp(:) = dcmplx(0d0,0d0)
      if (mynmo.eq.0) goto 100

      !
      ! Precompute arrays appearing in the denominator of
      ! Equation (40)
      !
      do jmo=1,mynmo
        vector(jmo) = omega - Evals(iwidx(jmo))
        etavec(jmo) = sign(1d0,evals(iwidx(jmo)))*ieta
      enddo

      !
      ! Loop over all points in the imaginary grid to
      ! obtain the value of I_{nn} (res(1)), and its derivative
      ! with respect to omega (res(2))
      !
!$omp parallel do private(iomega,jmo,factor) reduction(+:res)
      do igl=1,ngl
        iomega = imag*dsqrt(x(igl))
        do jmo=1,mynmo
          factor = 1d0/(vector(jmo) + iomega + etavec(jmo)*0d0)
          res(1) = res(1) - w(igl)*iw(jmo,igl)*factor
          res(2) = res(2) + w(igl)*iw(jmo,igl)*factor**2
        enddo
      enddo
!$omp end parallel do

      ! 
      ! Sum results from all MPI ranks.
      !
 100  call ga_zgop(1038,res,2,'+')


      !
      ! Put results back in the appropriate variables
      !
#ifdef GWCMPLX
      I = res(1)
      dI = res(2)
#else
      I = dble(res(1))
      dI = dble(res(2))
#endif
     
      end subroutine
