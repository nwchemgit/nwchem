      subroutine gw_cdgw_buildin_graph(iW,Evals,x,w,iWidx,omega,output,
     $                            ngl,myld,nqp)
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: ngl, myld, nqp
      integer,intent(in),dimension(0:myld,nqp) :: iWidx
      double precision, intent(in) :: omega
      double precision, intent(in) :: Evals(*)
      double precision, intent(in) :: x(*), w(*)
      double precision, intent(in), dimension(myld,ngl,nqp) :: iW
      double precision, intent(out),dimension(2,nqp) :: output

      integer igl,jmo,mynmo, iqp
      double precision w2,a
      double precision,dimension(myld) :: vector,etavec
      double precision,dimension(myld) :: factor,dfactor
      double precision,parameter :: eta = 0.001d0

      output(:,:) = 0d0

      do iqp=1,nqp

        mynmo = iwidx(0,iqp)
        if (mynmo.eq.0) cycle

        do jmo=1,mynmo
          vector(jmo) = omega - Evals(iwidx(jmo,iqp))
          etavec(jmo) = sign(eta,-evals(iwidx(jmo,iqp)))
        enddo

        do igl=1,ngl
          w2 = dsqrt(x(igl))
          do jmo=1,mynmo
            a = 1d0/(vector(jmo)**2 + (w2-etavec(jmo))**2)
            factor(jmo) =  vector(jmo)*a
            dfactor(jmo) = ((w2-etavec(jmo))**2 - vector(jmo)**2)*a**2
          enddo
          output(1,iqp) = output(1,iqp) - w(igl)*
     $                 dot_product(iw(1:mynmo,igl,iqp),factor(1:mynmo))
          output(2,iqp) = output(2,iqp) - w(igl)*
     $                dot_product(iw(1:mynmo,igl,iqp),dfactor(1:mynmo))
        enddo
      enddo

      call ga_dgop(mt_dbl,output,2*nqp,'+')

      end subroutine
