      subroutine argos_diana_grpdis(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "argos_diana_common.fh"
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
      integer i,j,jfrom,ito,ia,ja,number
      real*8 dist,dx,dy,dz
c
      igrp=igroups(igr,1)
      jgrp=igroups(igr,2)
c
      number=0
c
      if(ldef(igrp).lt.0) return
c
      write(lfngrp,1000) igr,igroups(igr,5)
 1000 format(2i5)
c
      ito=ldef(igrp)
      if(igrp.eq.jgrp) ito=ito-1
      do 1 i=1,ito
      ia=idef(igrp,i)
      jfrom=1
      if(igrp.eq.jgrp) jfrom=i+1
      do 2 j=jfrom,ldef(jgrp)
      ja=idef(jgrp,j)
      dx=abs(x(ia,1)-x(ja,1))
      dy=abs(x(ia,2)-x(ja,2))
      dz=abs(x(ia,3)-x(ja,3))
      if(igroups(igr,6).eq.1) then
      if(dz.gt.box(3)) dz=dz-box(3)
      elseif(igroups(igr,6).eq.2) then
      if(dx.gt.box(1)) dx=dx-box(1)
      if(dy.gt.box(2)) dy=dy-box(2)
      elseif(igroups(igr,6).eq.3) then
      if(dx.gt.box(1)) dx=dx-box(1)
      if(dy.gt.box(2)) dy=dy-box(2)
      if(dz.gt.box(3)) dz=dz-box(3)
      endif
      dist=sqrt(dx*dx+dy*dy+dz*dz)
      if(dist.ge.rgroups(igr,1).and.dist.le.rgroups(igr,2)) then
      write(lfngrp,1001)
     + imol(ia),sgmnam(ia)(11:16),sgmnam(ia)(1:5),sgmnam(ia)(6:10),
     + imol(ja),sgmnam(ja)(11:16),sgmnam(ja)(1:5),sgmnam(ja)(6:10),dist
 1001 format(2(i5,a6,' ',a5,':',a5,' '),f12.6)
      number=number+1
      endif
    2 continue
    1 continue
c
      write(lfngrp,1002) 0
 1002 format(i5)
c
      return
      end
c $Id$
