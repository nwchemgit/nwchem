      subroutine argos_diana_grploc(sgmnam,imol,isel,wt,x,igr,iwrk)
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
      integer isel(nsa),igr,imol(msa),iwrk(mxdef,mxnum,maxgrp)
      real*8 x(nsa,3)
c
      integer i,j,igrp,ia,ja,k,ntmp,ito
      real*8 dx,dy,dz,dist
c
      integer num,ndex(100)
      real*8 rd(100)
c
      igrp=igroup(igr,1)
      ito=ldef(igrp)
c
c      write(lfngrp,1000) time
c 1000 format('Atom distances for selected atoms at time ',f12.6)
c
      do 1 i=1,ito
      ia=idef(igrp,i)
      num=0
      do 2 ja=1,nsa
      if(ia.ne.ja.and.isel(ja).ne.0) then
      dx=abs(x(ia,1)-x(ja,1))
      dy=abs(x(ia,2)-x(ja,2))
      dz=abs(x(ia,3)-x(ja,3))
      if(igroup(igr,3).eq.1) then
      if(dz.gt.box(3)) dz=dz-box(3)
      elseif(igroup(igr,3).eq.2) then
      if(dx.gt.box(1)) dx=dx-box(1)
      if(dy.gt.box(2)) dy=dy-box(2)
      elseif(igroup(igr,3).eq.3) then
      if(dx.gt.box(1)) dx=dx-box(1)
      if(dy.gt.box(2)) dy=dy-box(2)
      if(dz.gt.box(3)) dz=dz-box(3)
      endif
      dist=sqrt(dx*dx+dy*dy+dz*dz)
      if(dist.ge.rgroup(igr,1).and.dist.le.rgroup(igr,2)) then
      if(num.lt.100) then
      num=num+1
      ndex(num)=ja
      rd(num)=dist
      endif
      endif
      endif
    2 continue
      if(num.gt.0) then
      do 3 j=1,num-1
      do 4 k=j+1,num
      if(ndex(j).gt.ndex(k)) then
      ntmp=ndex(j)
      ndex(j)=ndex(k)
      ndex(k)=ntmp
      endif
    4 continue
    3 continue
      do 5 j=1,num
      if(ndex(j).eq.0) goto 6
      if(ndex(j).ne.iwrk(igrp,j,i)) goto 6
    5 continue
      if(iwrk(igrp,num+1,i).eq.0) goto 1
    6 continue
      iwrk(igrp,num+1,i)=0
      do 7 j=1,num
      iwrk(igrp,j,i)=ndex(j)
    7 continue
      if(num.lt.11) then
      write(lfnloc,1001) time,ia,(ndex(j),j=1,num)
 1001 format(f12.6,11i6)
      else
      write(lfnloc,1002) time,ia,(ndex(j),j=1,num)
 1002 format(f12.6,11i6,/,(18x,10i6))
      endif
c      write(lfnloc,1001)
c     + imol(ia),sgmnam(ia)(11:16),sgmnam(ia)(1:5),sgmnam(ia)(6:10),
c     + (imol(ndex(j)),sgmnam(ndex(j))(11:16),sgmnam(ndex(j))(1:5),
c     + sgmnam(ndex(j))(6:10),j=1,num)
c 1001 format(2(i5,a6,' ',a5,':',a5,' '),/,(t25,i5,a6,' ',a5,':',a5,' '))
      endif
    1 continue
c
      return
      end
c $Id$
