      subroutine gw_cdgw_buildin_graph(iW,Evals,x,w,iWidx,omega,output,
     $                            eta,ngl,myld,nqp)
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

      do iqp=1,nqp

        mynmo = iwidx(0,iqp)
        if (mynmo.eq.0) cycle

        do jmo=1,mynmo
          vector(jmo) = omega - Evals(iwidx(jmo,iqp))
          etavec(jmo) = sign(eta,evals(iwidx(jmo,iqp)))
        enddo

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

      call ga_zgop(1038,res,2*nqp,'+')

#ifdef GWCMPLX
      output(:,:) = res(:,:)
#else
      output(:,:) = dble(res(:,:))
#endif

      end subroutine
