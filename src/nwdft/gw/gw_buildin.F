#ifdef GWCMPLX
#define GWTYP double complex
#else
#define GWTYP double precision
#endif

      subroutine gw_buildin(iW,Evals,x,w,iWidx,omega,I,dI,ngl,myld)
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: ngl, myld
      integer,intent(in),dimension(0:myld) :: iWidx
      double precision, intent(in) :: omega
      double precision, intent(in) :: Evals(*)
      double precision, intent(in) :: x(*), w(*)
      GWTYP, intent(in), dimension(myld,ngl) :: iW
      GWTYP, intent(out) :: I,dI

      integer igl,jmo,mynmo
      GWTYP,dimension(iwidx(0)) :: etavec
      double precision          :: vector(iwidx(0))

      double complex             :: iomega, factor, res(2), tmp(2)
      double complex, parameter  :: imag = dcmplx(0d0,1d0)
#ifdef GWCMPLX
      double precision,parameter :: eta = 0.001d0
      double complex, parameter  :: ieta = dcmplx(0d0,eta)
#else
      double precision,parameter  :: eta = 0.0000d0
      double precision, parameter :: ieta = eta
#endif


      mynmo = iwidx(0)
      res(:) = dcmplx(0d0,0d0)
      tmp(:) = dcmplx(0d0,0d0)
      if (mynmo.eq.0) goto 100

      do jmo=1,mynmo
        vector(jmo) = omega - Evals(iwidx(jmo))
        etavec(jmo) = sign(ieta,evals(iwidx(jmo)))
      enddo

!$omp parallel do private(iomega,jmo,factor) reduction(+:res)
      do igl=1,ngl

        iomega = imag*dsqrt(x(igl))

        do jmo=1,mynmo
          factor = 1d0/(vector(jmo) + iomega + etavec(jmo))
          res(1) = res(1) - w(igl)*iw(jmo,igl)*factor
          res(2) = res(2) + w(igl)*iw(jmo,igl)*factor**2
        enddo

      enddo

 100  call ga_zgop(mt_dcpl,res,2,'+')

#ifdef GWCMPLX
      I = res(1)
      dI = res(2)
#else
      I = dble(res(1))
      dI = dble(res(2))
#endif
     
      end subroutine
