      subroutine ana_groups(card,sgmnam,imol,isel,wt,x)
c
c $Id: ana_groups.F,v 1.6 2001-05-09 18:14:54 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*80 card
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),imol(msa)
      real*8 x(nsa,3)
c
      integer i
c
      ngroups=ngroups+1
      if(ngroups.gt.maxgrp) call errquit('Increase maxgrp',maxgrp)
c
      read(card(8:53),1000) (igroups(ngroups,i),i=1,4),
     + (rgroups(ngroups,i),i=1,2)
 1000 format(2i7,i5,i3,2f12.6)
c
      filgrp=card(54:80)
      if(filgrp(1:1).ne.' ') then
      open(unit=lfngrp,file=filgrp(1:index(filgrp,' ')-1),
     + form='formatted',status='unknown')
      if(igroups(ngroups,3).eq.1)
     + call ana_grpdis(sgmnam,imol,isel,wt,x,ngroups)
      ngroups=ngroups-1
      close(unit=lfngrp)
      endif
c
      return
      end
      subroutine ana_group(card,sgmnam,imol,isel,wt,x,iwrk)
c
c $Id: ana_groups.F,v 1.6 2001-05-09 18:14:54 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*80 card
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),imol(msa),iwrk(mxdef,mxnum,maxgrp)
      real*8 x(nsa,3)
c
      integer i,j
c
      ngroup=ngroup+1
      if(ngroup.gt.maxgrp) call errquit('Increase maxgrp',maxgrp)
c
      read(card(8:46),1000) (igroup(ngroup,i),i=1,3),
     + (rgroup(ngroup,i),i=1,2)
 1000 format(i7,i5,i3,2f12.6)
c
      do 1 i=1,mxdef
      do 2 j=1,mxnum
      iwrk(i,j,ngroup)=0
    2 continue
    1 continue
c
      return
      end
      subroutine ana_grpdis(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "ana_common.fh"
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
      integer i,j,ito,jfrom,ia,ja,number
      real*8 dist,dx,dy,dz
c
      igrp=igroups(igr,1)
      jgrp=igroups(igr,2)
c
      number=0
c
      if(ldef(igrp).lt.0) return
c
      write(lfngrp,1000) time
 1000 format('Atom distances for selected atoms at time ',f12.6)
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
      if(igroups(igr,4).eq.1) then
      if(dz.gt.box(3)) dz=dz-box(3)
      elseif(igroups(igr,4).eq.2) then
      if(dx.gt.box(1)) dx=dx-box(1)
      if(dy.gt.box(2)) dy=dy-box(2)
      elseif(igroups(igr,4).eq.3) then
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
      return
      end
      subroutine ana_lochdr(sgmnam,imol,isel)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer imol(msa),isel(nsa),igr
c
      integer i,num,igrp,ito,ia,j
c
      num=0
      do 1 i=1,nsa
      if(isel(i).ne.0) num=num+1
    1 continue
c
      write(lfnloc,1000) ngroup
 1000 format(2i7)
      do 4 j=1,ngroup
      igrp=igroup(j,1)
      ito=ldef(igrp)
      write(lfnloc,1000) ito,nsa
      do 2 i=1,ito
      ia=idef(igrp,i)
      write(lfnloc,1001) ia,imol(ia),sgmnam(ia)(11:16),sgmnam(ia)(1:5),
     + sgmnam(ia)(6:10)
    2 continue
    4 continue
      write(lfnloc,1000) num,nsa
      do 3 i=1,nsa
      if(isel(i).ne.0) then
      write(lfnloc,1001) i,imol(i),sgmnam(i)(11:16),sgmnam(i)(1:5),
     + sgmnam(i)(6:10)
 1001 format(i7,i5,a6,' ',a5,':',a5)
      endif
    3 continue
c
      return
      end
      subroutine ana_grploc(sgmnam,imol,isel,wt,x,igr,iwrk)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),igr,imol(msa),iwrk(mxdef,mxnum,maxgrp)
      real*8 x(nsa,3)
c
      integer i,j,igrp,ia,ja,ito,k,ntmp
      real*8 dx,dy,dz,dist
c
      integer num,ndx(100)
      real*8 rd(100)
c
      igrp=igroup(igr,1)
      ito=ldef(igrp)
c
c      write(lfnloc,1000) time
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
      ndx(num)=ja
      rd(num)=dist
      endif
      endif
      endif
    2 continue
      if(num.gt.0) then
      do 3 j=1,num-1
      do 4 k=j+1,num
      if(ndx(j).gt.ndx(k)) then
      ntmp=ndx(j)
      ndx(j)=ndx(k)
      ndx(k)=ntmp
      endif
    4 continue
    3 continue
      do 5 j=1,num
      if(ndx(j).eq.0) goto 6
      if(ndx(j).ne.iwrk(igrp,i,j)) goto 6
    5 continue
      if(iwrk(igrp,i,num+1).eq.0) goto 1
    6 continue
      iwrk(igrp,i,num+1)=0
      do 7 j=1,num
      iwrk(igrp,i,j)=ndx(j)
    7 continue
      if(num.lt.11) then
      write(lfnloc,1001) time,ia,(ndx(j),j=1,num)
 1001 format(f12.6,11i6)
      else
      write(lfnloc,1002) time,ia,(ndx(j),j=1,num)
 1002 format(f12.6,11i6,/,(18x,10i6))
      endif
c      write(lfnloc,1001)
c     + imol(ia),sgmnam(ia)(11:16),sgmnam(ia)(1:5),sgmnam(ia)(6:10),
c     + (imol(ndx(j)),sgmnam(ndx(j))(11:16),sgmnam(ndx(j))(1:5),
c     + sgmnam(ndx(j))(6:10),j=1,num)
c 1001 format(2(i5,a6,' ',a5,':',a5,' '),/,(t25,i5,a6,' ',a5,':',a5,' '))
      endif
    1 continue
c
      return
      end
      subroutine ana_histo(x,w,ih)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ih
      real*8 x(nsa,3),w(mwm,mwa,3)
c
      integer i,ia,ito,j,ndx
c
      dhis=1.2d0*box(3)/real(idhis(ih,2))
c
      ito=ldef(idhis(ih,1))
      if(ito.gt.0) then
      do 1 i=1,ito
      ia=idef(idhis(ih,1),i)
      ndx=(x(ia,3)-rhis)/dhis
      ihis(ndx,ih)=ihis(ndx,ih)+1
    1 continue
      else
      do 2 j=1,nwm
      do 3 i=1,-ito
      ia=idef(idhis(ih,1),i)
      ndx=(w(j,ia,3)-rhis)/dhis
      ihis(ndx,ih)=ihis(ndx,ih)+1
    3 continue
    2 continue
      endif
c
      return
      end      
      subroutine ana_order(x)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
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
      rord(i,1)=rord(i,1)+acos(cphi)
      rord(i,2)=rord(i,2)+cphi*cphi
c
    1 continue
c
      return
      end

