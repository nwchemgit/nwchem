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
      double precision w2
      double precision,dimension(iwidx(0)) :: vector,wmo
      double precision,dimension(iwidx(0)) :: factor
      double precision,dimension(-2:2)     :: res

      double precision,external :: ddot

      mynmo = iwidx(0)
      do jmo=1,iwidx(0)
        vector(jmo) = omega - Evals(iwidx(jmo))
      enddo

      res(:) = 0d0

      do igl=1,ngl
        w2 = x(igl)
        wmo(:) = iW(1:iwidx(0),igl)

        do istep=-2,2
          do jmo=1,iwidx(0)
            factor(jmo) =  (vector(jmo)+istep*fdstep)/
     $                     ((vector(jmo)+istep*fdstep)**2 + w2)
          enddo
          res(istep) = res(istep) - w(igl)*dot_product(wmo,factor)
        enddo
      enddo

      call ga_dgop(mt_dbl,res,5,'+')


      I = res(0)
      dI = (res(-2) - 8d0*res(-1) + 8d0*res(1) - res(2))/(12d0*fdstep)
      d2I = -res(-2)+16d0*res(-1)-30d0*res(0)+16d0*res(1)-res(2)
      d2I = d2I/(12d0*fdstep**2)
     
      end subroutine
