      subroutine argos_diana_rotate(xs,xw,nw)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer nw
c
      real*8 rv(3),rw(3),rotang,rx(3),ry(3),pifac
      integer i,j,k
c
      pifac=tan(1.0d0)/90.0d0
c
      do 1 k=1,nrot
      rv(1)=0.0d0
      rv(2)=0.0d0
      rv(3)=0.0d0
      rw(1)=0.0d0
      rw(2)=0.0d0
      rw(3)=0.0d0
      rw(irot(k))=1.0d0
      rotang=arot(k)*pifac
      do 2 i=1,nsa
      rx(1)=xs(i,1)
      rx(2)=xs(i,2)
      rx(3)=xs(i,3)
      ry(1)=xs(i,1)
      ry(2)=xs(i,2)
      ry(3)=xs(i,3)
      call rotate(rv,rw,rotang,ry,rx)
      xs(i,1)=rx(1)
      xs(i,2)=rx(2)
      xs(i,3)=rx(3)
    2 continue
      do 3 i=1,nw
      do 4 j=1,nwa
      rx(1)=xw(i,j,1)
      rx(2)=xw(i,j,2)
      rx(3)=xw(i,j,3)
      ry(1)=xw(i,j,1)
      ry(2)=xw(i,j,2)
      ry(3)=xw(i,j,3)
      call rotate(rv,rw,rotang,ry,rx)
      xw(i,j,1)=rx(1)
      xw(i,j,2)=rx(2)
      xw(i,j,3)=rx(3)
    4 continue
    3 continue
    1 continue
c
      return
      end
c $Id$
