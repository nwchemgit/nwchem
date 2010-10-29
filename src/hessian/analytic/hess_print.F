c
C $Id$
c
      subroutine hess_print
      implicit none
#include "hess_info.fh"
#include "mafdecls.fh"
      integer nn
      nn=n3xyz/3
      call hess_print0(dbl_mb(k_exy),nn)
      return
      end
      subroutine hess_print0(hessian,na)
      implicit real*8 (a-h,o-z)
#include "stdio.fh"
      dimension hessian(3,na,3,na)
c
c
      write(luout,100)
  100 format(/61('-'))
c
      write(luout,*)'            The Analytical Hessian Matrix '

      write(luout,*)'  '
c
c     do iat=1,na
c        do jat=iat,na
      do iat=1,na
         do jat=1  ,na
             write(luout,*) 'Atoms =', iat,jat
                xx=hessian(1,iat,1,jat)
                xy=hessian(1,iat,2,jat)
                xz=hessian(1,iat,3,jat)
                yx=hessian(2,iat,1,jat)
                yy=hessian(2,iat,2,jat)
                yz=hessian(2,iat,3,jat)
                zx=hessian(3,iat,1,jat)
                zy=hessian(3,iat,2,jat)
                zz=hessian(3,iat,3,jat)
                write(luout,103) 'xx xy xz : ',xx,xy,xz
                write(luout,103) 'yx yy yz : ',yx,yy,yz
                write(luout,103) 'zx zy zz : ',zx,zy,zz
c            endif
         enddo
      enddo
c
      write(luout,100)
c
  103 format(a10,3(f15.6,1x) )
  101 format(a10,26x,1(f15.6,1x) )
c
      end
c==============================================================
