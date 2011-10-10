      subroutine argos_diana_grpcogdis(sgmnam,imol,isel,wt,x,igr)
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
      integer igrp,jgrp
      integer i,ia,ja
      real*8 dist,dx,dy,dz
      real*8 cogi(3),cogj(3),factor
      real*8 boxh(3)
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      igrp=igroups(igr,1)
      jgrp=igroups(igr,2)
c
      if(ldef(igrp).lt.0) return
      if(ldef(jgrp).lt.0) return
c
      do 1 i=1,3
      cogi(i)=zero
      cogj(i)=zero
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
      dx=abs(cogi(1)-cogj(1))
      dy=abs(cogi(2)-cogj(2))
      dz=abs(cogi(3)-cogj(3))
c
      if(igroups(igr,6).eq.1) then
      if(dz.gt.boxh(3)) dz=dz-box(3)
      elseif(igroups(igr,6).eq.2) then
      if(dx.gt.boxh(1)) dx=dx-box(1)
      if(dy.gt.boxh(2)) dy=dy-box(2)
      elseif(igroups(igr,6).eq.3) then
      if(dx.gt.boxh(1)) dx=dx-box(1)
      if(dy.gt.boxh(2)) dy=dy-box(2)
      if(dz.gt.boxh(3)) dz=dz-box(3)
      endif
      dist=sqrt(dx*dx+dy*dy+dz*dz)
c
      write(lfngrp,1001) igr,igroups(igr,5),dist,dx,dy,dz
 1001 format(2i5,4f12.6)
c
      return
      end
c $Id$
