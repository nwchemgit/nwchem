      subroutine argos_diana_grpvectors(sgmnam,imol,isel,wt,x,igr)
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
      real*8 x(nsa,3),vi(3),vj(3)
c
      integer igrp,jgrp
      integer i,j,jfrom,ito,ia,ja,number
      real*8 dist,dx,dy,dz,dot
c
      igrp=igroups(igr,1)
      jgrp=igroups(igr,2)
      print*,'igrp,jgrp=',igrp,jgrp
c
      number=0
c
      if(ldef(igrp).lt.0) return
c
      write(lfngrp,1000) igr,igroups(igr,5)
 1000 format(2i5)
c
      if(igrp.eq.jgrp)
     + call md_abort('vectors: single atom list',0)
      if(ldef(igrp).ne.ldef(jgrp))
     + call md_abort('vectors: unequal atom lists',0)
      ito=ldef(igrp)
      do 1 i=1,ito
      ia=idef(igrp,i)
      ja=idef(jgrp,i)
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
c      write(lfngrp,1001)
c     + imol(ia),sgmnam(ia)(11:16),sgmnam(ia)(1:5),sgmnam(ia)(6:10),
c     + imol(ja),sgmnam(ja)(11:16),sgmnam(ja)(1:5),sgmnam(ja)(6:10),dist
c 1001 format(2(i5,a6,' ',a5,':',a5,' '),f12.6)
      write(lfngrp,1003) i,(x(ia,j),j=1,3),(x(ja,j)-x(ia,j),j=1,3),dist
 1003 format(i5,7f12.6)
      number=number+1
    1 continue
      do 2 i=1,ito
      ia=idef(igrp,i)
      ja=idef(jgrp,i)
      vi(1)=x(ia,1)-x(ja,1)
      vi(2)=x(ia,2)-x(ja,2)
      vi(3)=x(ia,3)-x(ja,3)
      do 3 j=i,ito
      ia=idef(igrp,j)
      ja=idef(jgrp,j)
      vj(1)=x(ia,1)-x(ja,1)
      vj(2)=x(ia,2)-x(ja,2)
      vj(3)=x(ia,3)-x(ja,3)
      dot=(vi(1)*vj(1)+vi(2)*vj(2)+vi(3)*vj(3))/
     + (sqrt(vi(1)**2+vi(2)**2+vi(3)**2)*
     +  sqrt(vj(1)**2+vj(2)**2+vj(3)**2))
      write(lfngrp,1004) i,j,dot
 1004 format(2i5,f12.6)
    3 continue
    2 continue
c    1 continue
c
      write(lfngrp,1002) 0
 1002 format(i5)
c
      return
      end
c $Id$
