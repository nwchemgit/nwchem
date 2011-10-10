      subroutine argos_diana_order(x)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3)
c
      integer i,ia,n,igrp,j
      real*8 xsijx,xsijy,xsijz,xskjx,xskjy,xskjz
      real*8 rsij2,rskj2,rsij2i,rskj2i,rsikji,cphi
c
      real*8 xa(3),xb(3),xc(3)
c
      do 1 i=1,nord
      igrp=idord(i,1)
      n=ldef(igrp)
      do 2 j=1,n
      ia=idef(igrp,j)
      xa(1)=x(ia,1)
      xa(2)=x(ia,2)
      xa(3)=x(ia,3)
    2 continue
      xa(1)=xa(1)/dble(n)-x(idord(i,3),1)
      xa(2)=xa(2)/dble(n)-x(idord(i,3),2)
      xa(3)=xa(3)/dble(n)-x(idord(i,3),3)
      igrp=idord(i,2)
      n=ldef(igrp)
      do 3 j=1,n
      ia=idef(igrp,j)
      xb(1)=x(ia,1)
      xb(2)=x(ia,2)
      xb(3)=x(ia,3)
    3 continue
      xb(1)=xb(1)/dble(n)-x(idord(i,3),1)
      xb(2)=xb(2)/dble(n)-x(idord(i,3),2)
      xb(3)=xb(3)/dble(n)-x(idord(i,3),3)
      xc(1)=x(idord(i,4),1)-x(idord(i,3),1)
      xc(2)=x(idord(i,4),2)-x(idord(i,3),2)
      xc(3)=x(idord(i,4),3)-x(idord(i,3),3)
c
      xsijx=xa(1)-xb(1)
      xskjx=xc(1)-xb(1)
      xsijy=xa(2)-xb(2)
      xskjy=xc(2)-xb(2)
      xsijz=xa(3)-xb(3)
      xskjz=xc(3)-xb(3)
c
      rsij2=xsijx*xsijx+xsijy*xsijy+xsijz*xsijz
      rskj2=xskjx*xskjx+xskjy*xskjy+xskjz*xskjz
      cphi=xsijx*xskjx+xsijy*xskjy+xsijz*xskjz
      rsij2i=one/rsij2
      rskj2i=one/rskj2
      rsikji=one/sqrt(rsij2*rskj2)
      cphi=cphi*rsikji
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
c
cc      rord(i,1)=rord(i,1)+acos(cphi)
cc     rord(i,2)=rord(i,2)+cphi*cphi
c
    1 continue
c
      return
      end
c $Id$
