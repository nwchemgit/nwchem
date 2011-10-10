      subroutine argos_diana_grpcogang(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),igr,imol(msa)
      real*8 x(nsa,3)
c
      integer igrp,jgrp,kgrp
      integer i,ia,ja,ka
      real*8 dx
      real*8 cogi(3),cogj(3),cogk(3),factor
      real*8 boxh(3)
      real*8 xsijx,xskjx,xsijy,xskjy,xsijz,xskjz,rsij2,rskj2
      real*8 cphi,rsij2i,rskj2i,rsikji,phi
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      igrp=igroups(igr,1)
      jgrp=igroups(igr,2)
      kgrp=igroups(igr,3)
c
      if(ldef(igrp).lt.0) return
      if(ldef(jgrp).lt.0) return
      if(ldef(kgrp).lt.0) return
c
      do 1 i=1,3
      cogi(i)=zero
      cogj(i)=zero
      cogk(i)=zero
    1 continue
c
      do 2 i=1,ldef(igrp)
      ia=idef(igrp,i)
      cogi(1)=cogi(1)+x(ia,1)
      cogi(2)=cogi(2)+x(ia,2)
      cogi(3)=cogi(3)+x(ia,3)
    2 continue
      factor=one/dble(ldef(igrp))
      cogi(1)=cogi(1)*factor
      cogi(2)=cogi(2)*factor
      cogi(3)=cogi(3)*factor
c
      do 3 i=1,ldef(jgrp)
      ja=idef(jgrp,i)
      cogj(1)=cogj(1)+x(ja,1)
      cogj(2)=cogj(2)+x(ja,2)
      cogj(3)=cogj(3)+x(ja,3)
    3 continue
      factor=one/dble(ldef(jgrp))
      cogj(1)=cogj(1)*factor
      cogj(2)=cogj(2)*factor
      cogj(3)=cogj(3)*factor
c
      do 4 i=1,ldef(kgrp)
      ka=idef(kgrp,i)
      cogk(1)=cogk(1)+x(ka,1)
      cogk(2)=cogk(2)+x(ka,2)
      cogk(3)=cogk(3)+x(ka,3)
    4 continue
      factor=one/dble(ldef(kgrp))
      cogk(1)=cogk(1)*factor
      cogk(2)=cogk(2)*factor
      cogk(3)=cogk(3)*factor
c
      if(igroups(igr,6).gt.0) then
      do 6 i=1,igroups(igr,6)
      dx=cogi(i)-cogj(i)
      if(dx.lt.-boxh(i)) cogi(i)=cogi(i)+box(i)
      if(dx.gt.boxh(i)) cogi(i)=cogi(i)-box(i)
      dx=cogk(i)-cogk(i)
      if(dx.lt.-boxh(i)) cogk(i)=cogk(i)+box(i)
      if(dx.gt.boxh(i)) cogk(i)=cogk(i)-box(i)
    6 continue
      endif
c
      xsijx=cogi(1)-cogj(1)
      xskjx=cogk(1)-cogj(1)
      xsijy=cogi(2)-cogj(2)
      xskjy=cogk(2)-cogj(2)
      xsijz=cogi(3)-cogj(3)
      xskjz=cogk(3)-cogj(3)
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
      phi=acos(cphi)
c
      write(lfngrp,1001) igr,igroups(igr,5),phi
 1001 format(2i5,4f12.6)
c
      return
      end
 
c $Id$
