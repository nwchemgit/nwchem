      subroutine gw_buildin(iW,Evals,x,w,iWidx,omega,I,dI,ngl,myld)
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: ngl, myld
      integer,intent(in),dimension(0:myld) :: iWidx
      double precision, intent(in) :: omega
      double precision, intent(in) :: Evals(*)
      double precision, intent(in) :: x(*), w(*)
      double precision, intent(in), dimension(myld,ngl) :: iW
      double precision, intent(out) :: I,dI

      integer igl,jmo,mynmo,istep
      double precision w2,a
      double precision,dimension(iwidx(0)) :: vector,wmo,etavec
      double precision,dimension(iwidx(0)) :: factor, dfactor
      double precision,dimension(2)     :: res
      double precision,parameter :: eta = 0.001d0

      integer,external :: ga_nodeid

      mynmo = iwidx(0)
      res(:) = 0d0
      if (mynmo.eq.0) goto 100

      do jmo=1,iwidx(0)
        vector(jmo) = omega - Evals(iwidx(jmo))
        etavec(jmo) = sign(eta,-evals(iwidx(jmo)))
      enddo


      do igl=1,ngl

        w2 = dsqrt(x(igl))
        wmo(:) = iW(1:iwidx(0),igl)

        do jmo=1,iwidx(0)
          a = 1d0/(vector(jmo)**2 + (w2-etavec(jmo))**2)
          factor(jmo) =  vector(jmo)*a
          dfactor(jmo) = ((w2-etavec(jmo))**2 - vector(jmo)**2)*a**2
        enddo
        res(1) = res(1) - w(igl)*dot_product(wmo,factor)
        res(2) = res(2) - w(igl)*dot_product(wmo,dfactor)

      enddo

 100  call ga_dgop(mt_dbl,res,3,'+')

      I = res(1)
      dI = res(2)
     
      end subroutine
