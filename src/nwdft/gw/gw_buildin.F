      subroutine gw_buildin(iW,Evals,x,w,iWidx,omega,I,dI,d2I,fdstep,
     $                      ngl,myld)
      implicit none
#include "mafdecls.fh"      
      integer,intent(in) :: ngl, myld
      integer,intent(in),dimension(0:myld) :: iWidx
      double precision, intent(in) :: fdstep, omega
      double precision, intent(in) :: Evals(*)
      double precision, intent(in) :: x(*), w(*)
      double precision, intent(in), dimension(myld,ngl) :: iW
      double precision, intent(out) :: I,dI,d2I

      integer igl,jmo,mynmo,istep
      double precision w2,a
      double precision,dimension(iwidx(0)) :: vector,wmo,etavec
      double precision,dimension(iwidx(0)) :: factor, dfactor, d2factor
      double precision,dimension(0:2)     :: res
      double precision,parameter :: eta = 0.001d0

      integer,external :: ga_nodeid

      mynmo = iwidx(0)
      do jmo=1,iwidx(0)
        vector(jmo) = omega - Evals(iwidx(jmo))
        etavec(jmo) = sign(eta,-evals(iwidx(jmo)))
      enddo

      res(:) = 0d0

      do igl=1,ngl
        if (iwidx(0).eq.0) exit

        w2 = dsqrt(x(igl))
        wmo(:) = iW(1:iwidx(0),igl)

        do jmo=1,iwidx(0)
          a = 1d0/(vector(jmo)**2 + (w2-etavec(jmo))**2)
          factor(jmo) =  vector(jmo)*a
          dfactor(jmo) = ((w2-etavec(jmo))**2 - vector(jmo)**2)*a**2
          d2factor(jmo) = 2d0*vector(jmo)*(vector(jmo)**2 - 
     $                    3d0*(w2-etavec(jmo))**2)*a**3
        enddo
        res(0) = res(0) - w(igl)*dot_product(wmo,factor)
        res(1) = res(1) - w(igl)*dot_product(wmo,dfactor)
        res(2) = res(2) - w(igl)*dot_product(wmo,d2factor)

      enddo

      call ga_dgop(mt_dbl,res,3,'+')

      I = res(0)
      dI = res(1)
      d2I = res(2)
     
      end subroutine
