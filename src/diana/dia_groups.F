      subroutine dia_groups(sgmnam,imol,isel,wt,x)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),imol(msa)
      real*8 x(nsa,3)
c
      integer i
c
      ngroups=ngroups+1
      if(ngroups.gt.maxgrp) call md_abort('Increase maxgrp',maxgrp)
c
      if(card(1:6).eq.'groups') then
      read(card(8:67),1001) (igroups(ngroups,i),i=1,6),
     + (rgroups(ngroups,i),i=1,2)
 1001 format(4i7,i5,i3,2f12.6)
      endif
      print*,'++++++++++',igroups(ngroups,5) 
c
      filgrp=card(68:80)
      if(filgrp(1:1).ne.' ') then
      open(unit=lfngrp,file=filgrp(1:index(filgrp,' ')-1),
     + form='formatted',status='unknown')
      if(igroups(ngroups,5).eq.1)
     + call dia_grpcogdis(sgmnam,imol,isel,wt,x,ngroups)
      if(igroups(ngroups,5).eq.2)
     + call dia_grpdis(sgmnam,imol,isel,wt,x,ngroups)
      if(igroups(ngroups,5).eq.4)
     + call dia_grpcogang(sgmnam,imol,isel,wt,x,ngroups)
      if(igroups(ngroups,5).eq.5)
     + call dia_grpvectors(sgmnam,imol,isel,wt,x,ngroups)
      ngroups=ngroups-1
      close(unit=lfngrp)
      endif
c
      return
      end
      subroutine dia_group(sgmnam,imol,isel,wt,x,iwrk)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),imol(msa),iwrk(mxdef,mxnum,maxgrp)
      real*8 x(nsa,3)
c
      integer i,j
c
      ngroup=ngroup+1
      if(ngroup.gt.maxgrp) call md_abort('Increase maxgrp',maxgrp)
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
      subroutine dia_grpdis(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "dia_common.fh"
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
      subroutine dia_grpvectors(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "dia_common.fh"
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
      subroutine dia_lochdr(sgmnam,imol,isel)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer imol(msa),isel(nsa)
c
      integer i,num,igrp,ia,j,ito
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
      subroutine dia_grploc(sgmnam,imol,isel,wt,x,igr,iwrk)
c
      implicit none
c
#include "dia_common.fh"
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
      subroutine dia_histo(x,w,ihi)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ihi,ito
      real*8 x(nsa,3),w(mwm,mwa,3)
c
      integer i,ia,j,ihndx
c
      dhis=1.2d0*box(3)/real(idhis(ihi,2))
c
      ito=ldef(idhis(ihi,1))
      if(ito.gt.0) then
      do 1 i=1,ito
      ia=idef(idhis(ihi,1),i)
      ihndx=(x(ia,3)-rhis)/dhis
      ihis(ihndx,ihi)=ihis(ihndx,ihi)+1
    1 continue
      else
      do 2 j=1,nwm
      do 3 i=1,-ito
      ia=idef(idhis(ihi,1),i)
      ihndx=(w(j,ia,3)-rhis)/dhis
      ihis(ihndx,ihi)=ihis(ihndx,ihi)+1
    3 continue
    2 continue
      endif
c
      return
      end      
      subroutine dia_order(x)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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
      subroutine dia_grpcogdis(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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
      subroutine dia_grpcogang(sgmnam,imol,isel,wt,x,igr)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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

