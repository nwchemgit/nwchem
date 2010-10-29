      logical function argos_prepare_wrtpdb(lfnout,lfnpdb,lrgpdb,filpdb,
     + iopt,box,
     + num,amass,mat,
     + csa,isat,isgm,imol,ifra,xs,vs,msa,nsa,
     + cwa,iwat,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc,slvnam,iropt,irrand,nxrep,nyrep,nzrep,drep,
     + msb,nsb,idsb,rdist,nskip,iskip,lang,lfnmrg,nmerge,xmerge,filmrg,
     + irenum,invert,ihop,ips)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
      character*2 argos_prepare_atnam
      real*8 argos_prepare_atsiz
      integer argos_prepare_merge
      external argos_prepare_atnam,argos_prepare_atsiz,
     + argos_prepare_merge
c
      integer lfnout,lfnpdb,iopt,msa,nsa,mwm,nwm,mwa,nwa,mat,mwmc,nwmc
      integer lfnmrg,lrgpdb
      character*255 filpdb,filmrg(100)
      integer num(mat),isat(msa),isgm(msa),imol(msa),ifra(msa),iwat(mwa)
      integer ihop(msa),ips(msa)
      character*16 cwa(mwa),csa(msa)
      character*3 slvnam
      real*8 amass(mat),xs(3,msa),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),box(3)
      integer iropt,irrand,nxrep,nyrep,nzrep,nskip
      real*8 rdist,drep(3),xmerge(3,100)
      integer msb,nsb,nmols,irenum,invert,nmerge
      integer idsb(2,msb),iskip(3,nskip),lang(*)
c
      integer length,ioff
c
      integer i,j,k,l,nats,nums,numslv,isoff,iaoff,natom,isegm,nsegm
      integer ixrep,iyrep,izrep
      real*8 temp,dbx,dby,dbz,dist,dbz0,angle
      real*8 o(3),px(3),py(3),pz(3),x(3),y(3),v(3)
      character*3 sname
      character*1 charr
      integer icharr,ilang,inum
c
      real*8 dumdst
c
      dumdst=0.025d0
c
      isoff=0
      iaoff=0
      inum=-1
      if(irenum.ne.0) inum=0
      charr=' '
      icharr=ichar(charr)
      dbz0=0.0d0
c
      o(1)=0.0d0
      o(2)=0.0d0
      o(3)=0.0d0
      px(1)=1.0d0
      px(2)=0.0d0
      px(3)=0.0d0
      py(1)=0.0d0
      py(2)=1.0d0
      py(3)=0.0d0
      pz(1)=0.0d0
      pz(2)=0.0d0
      pz(3)=1.0d0
c
c     open PDB file
c
      length=index(filpdb,' ')-1
      open(unit=lfnpdb,file=filpdb(1:length),form='formatted',
     + status='unknown',err=9999)
c
      if(nxrep.lt.1) nxrep=1
      if(nyrep.lt.1) nyrep=1
      if(nzrep.lt.-2) nzrep=1
      nums=0
c
      if(nxrep.le.1.and.nyrep.le.1.and.iabs(nzrep).le.1) then
c
      if(lrgpdb.eq.0) then
      write(lfnpdb,1000) (1.0d1*box(i),i=1,3),90.0,90.0,90.0
 1000 format(
     + 'HEADER',/,
     + 'TITLE     ',/,
     + 'TITLE    2',/,
     + 'TITLE    3',/,
     + 'REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.1, 25-OCT-1996',/,
     + 'CRYST1',3f9.3,3f7.2)
      else
      write(lfnpdb,2000) (1.0d1*box(i),i=1,3),90.0,90.0,90.0
 2000 format(
     + 'HEADER',/,
     + 'TITLE     ',/,
     + 'TITLE    2',/,
     + 'TITLE    3',/,
     + 'REMARK   4 YYYY DOES NOT COMPLY WITH FORMAT V. 2.1,',
     + ' 25-OCT-1996',/,
     + 'CRYST1',3f9.3,3f7.2,/,
     + 'LRGPDB')
      endif
c
      ioff=0
c      if(filmrg(1:1).ne.' ') then
      if(nmerge.gt.0) then
      do 6060 k=1,nmerge
      ioff=ioff+argos_prepare_merge(lfnmrg,xmerge(1,k),filmrg(k),
     + lfnout,lfnpdb,slvnam,inum,lrgpdb)
 6060 continue
      ioff=inum
c      ioff=argos_prepare_merge(lfnmrg,filmrg,lfnout,lfnpdb,slvnam,inum)
      if(ioff.lt.0) call md_abort('Error in argos_prepare_merge',ioff)
      if(inum.ge.0) ioff=ioff+1-isgm(1)
      else
      if(irenum.ne.0) ioff=-isgm(1)+1
      endif
c
      nums=0
      do 11 k=1,nsb
      i=idsb(1,k)
      j=idsb(2,k)
c      write(*,'(5i5)') i,j,k,isgm(i),isgm(j)
      if(isgm(i).ne.isgm(j)) then
      dist=sqrt((xs(1,j)-xs(1,i))*(xs(1,j)-xs(1,i))+
     + (xs(2,j)-xs(2,i))*(xs(2,j)-xs(2,i))+
     + (xs(3,j)-xs(3,i))*(xs(3,j)-xs(3,i)))
      if(dist.gt.argos_prepare_atsiz(num(isat(j)))+
     + argos_prepare_atsiz(num(isat(i)))) then
      write(lfnpdb,1101) csa(i)(11:14),csa(i)(1:3),' ',isgm(i)+ioff,
     + csa(j)(11:14),csa(j)(1:3),isgm(j)+ioff
 1101 format('LINK',t13,a4,1x,a3,1x,a1,i4,16x,a4,1x,a3,2x,i4)
      endif
      endif
   11 continue
      do 1 i=1,nsa
      if(ips(i).ge.0.or.ihop(i).eq.0) then
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      else
      x(1)=xs(1,i)-xs(1,i+ihop(i))
      x(2)=xs(2,i)-xs(2,i+ihop(i))
      x(3)=xs(3,i)-xs(3,i+ihop(i))
      dist=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      dist=dumdst/dist
      x(1)=xs(1,i+ihop(i))+dist*x(1)
      x(2)=xs(2,i+ihop(i))+dist*x(2)
      x(3)=xs(3,i+ihop(i))+dist*x(3)
      endif
      temp=4.009103873d1*amass(isat(i))*
     + (vs(1,i)*vs(1,i)+vs(2,i)*vs(2,i)+vs(3,i)*vs(3,i))
      if(temp.gt.999.0) temp=999.989
      sname=csa(i)(1:3)
      if(sname(1:1).eq.'_') sname(1:1)=' '
      if(sname(2:2).eq.'_') sname(2:2)=' '
      if(sname(3:3).eq.'_') sname(3:3)=' '
      if(lrgpdb.eq.0) then
      write(lfnpdb,1001) i,csa(i)(11:14),sname,' ',
     + isgm(i)+ioff,(1.0d1*x(k),k=1,3),temp,
     + argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,2001) i,csa(i)(11:14),sname,
     + isgm(i)+ioff,(1.0d1*x(k),k=1,3),temp,
     + argos_prepare_atnam(num(isat(i)))
      endif
 1001 format('ATOM',i7,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,6x,f6.2,4x,a2)
 2001 format('ATOM',i7,1x,a4,1x,a3,i6,4x,3f8.3,6x,f6.2,4x,a2)
      if(isgm(i).gt.nums) nums=isgm(i)
    1 continue
      else
      dbz0=0.0d0
      charr='A'
      icharr=ichar(charr)
      if(nzrep.lt.0) then
      if(nsa.gt.0.and.nwm.eq.0) dbz0=xs(3,1)
      if(nsa.eq.0.and.nwm.gt.0) dbz0=xw(3,1,1)
      if(nsa.gt.0.and.nwm.gt.0) dbz=min(xs(3,1),xw(3,1,1))
      if(invert.eq.0) then
      do 31 i=1,nsa
      dbz0=min(dbz0,xs(3,i))
   31 continue
      do 32 i=1,nwm
      do 33 j=1,nwa
      dbz0=min(dbz0,xw(3,j,i))
   33 continue
   32 continue
      dbz0=dbz0-0.5d0*rdist
      else
      do 41 i=1,nsa
      dbz0=max(dbz0,xs(3,i))
   41 continue
      do 42 i=1,nwm
      do 43 j=1,nwa
      dbz0=max(dbz0,xw(3,j,i))
   43 continue
   42 continue
      dbz0=dbz0+0.5d0*rdist
      endif
      endif
c
      nmols=0
      do 34 i=1,nsa
      if(iropt.eq.3) then
      nmols=max(nmols,ifra(i))
      else
      nmols=max(nmols,imol(i))
      endif
   34 continue
c
      if(lrgpdb.eq.0) then
      write(lfnpdb,1000) 1.0d1*dble(nxrep)*drep(1),
     + 1.0d1*dble(nyrep)*drep(2),1.0d1*dble(nzrep)*drep(3),
     + 90.0,90.0,90.0
      else
      write(lfnpdb,2000) 1.0d1*dble(nxrep)*drep(1),
     + 1.0d1*dble(nyrep)*drep(2),1.0d1*dble(nzrep)*drep(3),
     + 90.0,90.0,90.0
      endif
c
      ioff=0
      inum=0
      if(nmerge.gt.0) then
      do 6061 k=1,nmerge
      ioff=ioff+argos_prepare_merge(lfnmrg,xmerge(1,k),filmrg(k),
     + lfnout,lfnpdb,slvnam,inum,lrgpdb)
 6061 continue
      ioff=inum
c      if(filmrg(1:1).ne.' ') then
c      ioff=argos_prepare_merge(lfnmrg,filmrg,lfnout,lfnpdb,slvnam,inum)
      if(ioff.lt.0) call md_abort('Error 2 in argos_prepare_merge',0)
      endif
c
      if(iropt.le.1) then
      ilang=0
      do 4 izrep=1,iabs(nzrep)
      if(nzrep.gt.0) then
      dbz=(dble(izrep)-0.5d0*dble(nzrep+1))*drep(3)
      else
      dbz=dbz0
      endif
      do 5 iyrep=1,nyrep
      dby=(dble(iyrep)-0.5d0*dble(nyrep+1))*drep(2)
      do 6 ixrep=1,nxrep
      dbx=(dble(ixrep)-0.5d0*dble(nxrep+1))*drep(1)
c
      do 66 k=1,nskip
      if(ixrep.eq.iskip(1,k).and.iyrep.eq.iskip(2,k).and.
     + izrep.eq.iskip(3,k)) goto 6
   66 continue
c
      angle=8.0d0*atan(1.0d0)*util_random(0)
      ilang=ilang+1
      lang(ilang)=angle
c
      if(iropt.eq.0) then
      icharr=ichar(' ')
      else
      isoff=0
      endif
c
      do 111 k=1,nsb
      i=idsb(1,k)
      j=idsb(2,k)
      if(isgm(i).ne.isgm(j)) then
      dist=sqrt((xs(1,j)-xs(1,i))*(xs(1,j)-xs(1,i))+
     + (xs(2,j)-xs(2,i))*(xs(2,j)-xs(2,i))+
     + (xs(3,j)-xs(3,i))*(xs(3,j)-xs(3,i)))
      if(dist.gt.argos_prepare_atsiz(num(isat(j)))+
     + argos_prepare_atsiz(num(isat(i)))) then
      write(lfnpdb,1101) csa(i)(11:14),csa(i)(1:3),char(icharr),
     + isgm(i)+ioff+isoff,csa(j)(11:14),csa(j)(1:3),isgm(j)+ioff+isoff
      endif
      endif
  111 continue
c
      nums=0
c
      do 7 i=1,nsa
c
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      if(ips(i).ge.0.or.ihop(i).eq.0) then
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      else
      x(1)=xs(1,i)-xs(1,i+ihop(i))
      x(2)=xs(2,i)-xs(2,i+ihop(i))
      x(3)=xs(3,i)-xs(3,i+ihop(i))
      dist=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      dist=dumdst/dist
      x(1)=xs(1,i+ihop(i))+dist*x(1)
      x(2)=xs(2,i+ihop(i))+dist*x(2)
      x(3)=xs(3,i+ihop(i))+dist*x(3)
      endif
      y(1)=x(1)
      y(2)=x(2)
      y(3)=x(3)
      if(irrand.eq.1.or.irrand.eq.4) call rotate(o,px,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
      if(irrand.eq.2.or.irrand.eq.4) call rotate(o,py,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
      if(irrand.eq.3.or.irrand.eq.4) call rotate(o,pz,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
c
      temp=4.009103873d1*amass(isat(i))*
     + (vs(1,i)*vs(1,i)+vs(2,i)*vs(2,i)+vs(3,i)*vs(3,i))
      if(temp.gt.999.0) temp=999.989
      if(lrgpdb.eq.0) then
      if(nzrep.gt.0) then
      write(lfnpdb,1001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)+dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,1001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,1001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,-1.0d1*(y(1)+dbx),
     + -1.0d1*(y(2)+dby),-1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      endif
      endif
      else
      if(char(icharr).ne.' ') call md_abort('Error in LRGPDB ',icharr)
      if(nzrep.gt.0) then
      write(lfnpdb,2001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)+dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,2001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,2001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,-1.0d1*(y(1)+dbx),
     + -1.0d1*(y(2)+dby),-1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      endif
      endif
      endif
      if(isgm(i).gt.nums) nums=isgm(i)
    7 continue
      isoff=isoff+nums
      iaoff=iaoff+nsa
      icharr=icharr+1
    6 continue
    5 continue
    4 continue
c
      else
      natom=0
      nsegm=0
      isegm=0
      isoff=0
      do 20 l=1,nmols
      ilang=0
      do 14 izrep=1,iabs(nzrep)
      if(nzrep.gt.0) then
      dbz=(dble(izrep)-0.5d0*dble(nzrep+1))*drep(3)
      else
      dbz=dbz0
      endif
      do 15 iyrep=1,nyrep
      dby=(dble(iyrep)-0.5d0*dble(nyrep+1))*drep(2)
      do 16 ixrep=1,nxrep
      dbx=(dble(ixrep)-0.5d0*dble(nxrep+1))*drep(1)
c
      do 67 k=1,nskip
      if(ixrep.eq.iskip(1,k).and.iyrep.eq.iskip(2,k).and.
     + izrep.eq.iskip(3,k)) goto 16
   67 continue
c
      ilang=ilang+1
c
      if(l.eq.1) then
      angle=8.0d0*atan(1.0d0)*util_random(0)
      lang(ilang)=angle
      else
      angle=lang(ilang)
      endif
c
      icharr=ichar(' ')
c
      do 212 i=1,nsa
      if(iropt.eq.3) then
      if(ifra(i).eq.l) then
      isoff=nsegm-isgm(i)+1
      goto 213
      endif
      else
      if(imol(i).eq.l) then
      isoff=nsegm-isgm(i)+1
      goto 213
      endif
      endif
  212 continue
  213 continue
c
      do 211 k=1,nsb
      i=idsb(1,k)
      j=idsb(2,k)
      if(iropt.eq.3) then
      if(ifra(i).eq.l.and.ifra(j).eq.l.and.isgm(i).ne.isgm(j)) then
      dist=sqrt((xs(1,j)-xs(1,i))*(xs(1,j)-xs(1,i))+
     + (xs(2,j)-xs(2,i))*(xs(2,j)-xs(2,i))+
     + (xs(3,j)-xs(3,i))*(xs(3,j)-xs(3,i)))
      if(dist.gt.argos_prepare_atsiz(num(isat(j)))+
     + argos_prepare_atsiz(num(isat(i)))) then
      write(lfnpdb,1101) csa(i)(11:14),csa(i)(1:3),char(icharr),
     + isgm(i)+ioff+isoff,csa(j)(11:14),csa(j)(1:3),isgm(j)+ioff+isoff
      endif
      endif
      else
      if(imol(i).eq.l.and.imol(j).eq.l.and.isgm(i).ne.isgm(j)) then
      dist=sqrt((xs(1,j)-xs(1,i))*(xs(1,j)-xs(1,i))+
     + (xs(2,j)-xs(2,i))*(xs(2,j)-xs(2,i))+
     + (xs(3,j)-xs(3,i))*(xs(3,j)-xs(3,i)))
      if(dist.gt.argos_prepare_atsiz(num(isat(j)))+
     + argos_prepare_atsiz(num(isat(i)))) then
      write(lfnpdb,1101) csa(i)(11:14),csa(i)(1:3),char(icharr),
     + isgm(i)+ioff+isoff,csa(j)(11:14),csa(j)(1:3),isgm(j)+ioff+isoff
      endif
      endif
      endif
  211 continue
c
      nums=0
c
      do 17 i=1,nsa
c
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      if(ips(i).ge.0.or.ihop(i).eq.0) then
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      else
      x(1)=xs(1,i)-xs(1,i+ihop(i))
      x(2)=xs(2,i)-xs(2,i+ihop(i))
      x(3)=xs(3,i)-xs(3,i+ihop(i))
      dist=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      dist=dumdst/dist
      x(1)=xs(1,i+ihop(i))+dist*x(1)
      x(2)=xs(2,i+ihop(i))+dist*x(2)
      x(3)=xs(3,i+ihop(i))+dist*x(3)
      endif
      y(1)=x(1)
      y(2)=x(2)
      y(3)=x(3)
      if(irrand.eq.1.or.irrand.eq.4) call rotate(o,px,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
      if(irrand.eq.2.or.irrand.eq.4) call rotate(o,py,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
      if(irrand.eq.3.or.irrand.eq.4) call rotate(o,pz,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
c
      if(iropt.eq.3) then
      if(ifra(i).eq.l) then
      temp=4.009103873d1*amass(isat(i))*
     + (vs(1,i)*vs(1,i)+vs(2,i)*vs(2,i)+vs(3,i)*vs(3,i))
      if(temp.gt.999.0) temp=999.989
      natom=natom+1
      if(lrgpdb.eq.0) then
      if(nzrep.gt.0) then
      write(lfnpdb,1001) natom,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)+dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,1001) natom,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,1001) natom,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,-1.0d1*(y(1)+dbx),
     + -1.0d1*(y(2)+dby),-1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      endif
      endif
      else
      if(char(icharr).ne.' ') call md_abort('Error in LRGPDB ',icharr)
      if(nzrep.gt.0) then
      write(lfnpdb,2001) natom,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)+dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,2001) natom,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,2001) natom,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,-1.0d1*(y(1)+dbx),
     + -1.0d1*(y(2)+dby),-1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      endif
      endif
      endif
      nsegm=isgm(i)+isoff
      endif
      else
      if(imol(i).eq.l) then
      temp=4.009103873d1*amass(isat(i))*
     + (vs(1,i)*vs(1,i)+vs(2,i)*vs(2,i)+vs(3,i)*vs(3,i))
      if(temp.gt.999.0) temp=999.989
      natom=natom+1
      if(lrgpdb.eq.0) then
      if(nzrep.gt.0) then
      write(lfnpdb,1001) natom,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)+dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,1001) natom,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,1001) natom,csa(i)(11:14),csa(i)(1:3),
     + char(icharr),isgm(i)+ioff+isoff,-1.0d1*(y(1)+dbx),
     + -1.0d1*(y(2)+dby),-1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      endif
      endif
      else
      if(char(icharr).ne.' ') call md_abort('Error in LRGPDB ',icharr)
      if(nzrep.gt.0) then
      write(lfnpdb,2001) natom,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)+dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,2001) natom,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,1.0d1*(y(1)+dbx),
     + 1.0d1*(y(2)+dby),1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      else
      write(lfnpdb,2001) natom,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+ioff+isoff,-1.0d1*(y(1)+dbx),
     + -1.0d1*(y(2)+dby),-1.0d1*(y(3)-dbz),
     + temp,argos_prepare_atnam(num(isat(i)))
      endif
      endif
      endif
      nsegm=isgm(i)+isoff
      endif
      endif
   17 continue
   16 continue
   15 continue
   14 continue
   20 continue
      endif
c
      endif
c
      write(lfnpdb,1002)
 1002 format('TER')
c
      numslv=iopt
      if(iopt.lt.0) numslv=nwmc
      if(iopt.gt.nwm+nwmc) numslv=nwm+nwmc
c
      if(numslv.gt.0) then
c
      if(nxrep.le.1.and.nyrep.le.1.and.iabs(nzrep).le.1) then
c
      nats=nsa
      do 2 i=1,numslv
      do 3 j=1,nwa
      nats=nats+1
      if(i.le.nwmc) then
      temp=4.009103873d1*amass(iwat(j))*(vwc(1,j,i)*vwc(1,j,i)+
     + vwc(2,j,i)*vwc(2,j,i)+vwc(3,j,i)*vwc(3,j,i))
      if(temp.gt.999.0) temp=999.989
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,(1.0d1*xwc(k,j,i),k=1,3),temp,
     + argos_prepare_atnam(num(iwat(j)))
      else
      temp=4.009103873d1*amass(iwat(j))*(vw(1,j,i-nwmc)*vw(1,j,i-nwmc)+
     + vw(2,j,i-nwmc)*vw(2,j,i-nwmc)+vw(3,j,i-nwmc)*vw(3,j,i-nwmc))
      if(temp.gt.999.0) temp=999.989
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,(1.0d1*xw(k,j,i-nwmc),k=1,3),temp,
     + argos_prepare_atnam(num(iwat(j)))
      endif
 1003 format('ATOM',i7,1x,a4,1x,a3,i6,4x,3f8.3,6x,f6.2,4x,a2)
    3 continue
    2 continue
c
      else
c
      nats=nxrep*nyrep*iabs(nzrep)*nsa
      ilang=0
      do 44 izrep=1,iabs(nzrep)
      if(nzrep.gt.0) then
      dbz=(dble(izrep)-0.5d0*dble(nzrep+1))*drep(3)
      else
      dbz=dbz0
      endif
      do 45 iyrep=1,nyrep
      dby=(dble(iyrep)-0.5d0*dble(nyrep+1))*drep(2)
      do 46 ixrep=1,nxrep
      dbx=(dble(ixrep)-0.5d0*dble(nxrep+1))*drep(1)
c
      ilang=ilang+1
      angle=lang(ilang)
c
      nums=isoff
      do 47 i=1,numslv
      do 48 j=1,nwa
      nats=nats+1
c
      if(i.le.nwmc) then
      x(1)=xwc(1,j,i)
      x(2)=xwc(2,j,i)
      x(3)=xwc(3,j,i)
      y(1)=xwc(1,j,i)
      y(2)=xwc(2,j,i)
      y(3)=xwc(3,j,i)
      v(1)=vwc(1,j,i)
      v(2)=vwc(2,j,i)
      v(3)=vwc(3,j,i)
      else
      x(1)=xw(1,j,i-nwmc)
      x(2)=xw(2,j,i-nwmc)
      x(3)=xw(3,j,i-nwmc)
      y(1)=xw(1,j,i-nwmc)
      y(2)=xw(2,j,i-nwmc)
      y(3)=xw(3,j,i-nwmc)
      v(1)=vw(1,j,i-nwmc)
      v(2)=vw(2,j,i-nwmc)
      v(3)=vw(3,j,i-nwmc)
      endif
      if(irrand.eq.1.or.irrand.eq.4) call rotate(o,px,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
      if(irrand.eq.2.or.irrand.eq.4) call rotate(o,py,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
      if(irrand.eq.3.or.irrand.eq.4) call rotate(o,pz,angle,x,y)
      x(1)=y(1)
      x(2)=y(2)
      x(3)=y(3)
c
      if(nzrep.gt.0) then
      temp=4.009103873d1*amass(iwat(j))*(v(1)*v(1)+
     + v(2)*v(2)+v(3)*v(3))
      if(temp.gt.999.0) temp=999.989
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(y(1)+dbx),1.0d1*(y(2)+dby),
     + 1.0d1*(y(3)+dbz),temp,argos_prepare_atnam(num(iwat(j)))
      else
      if(izrep.eq.1) then
      temp=4.009103873d1*amass(iwat(j))*(v(1)*v(1)+
     + v(2)*v(2)+v(3)*v(3))
      if(temp.gt.999.0) temp=999.989
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(y(1)+dbx),1.0d1*(y(2)+dby),
     + 1.0d1*(y(3)-dbz),temp,argos_prepare_atnam(num(iwat(j)))
      else
      temp=4.009103873d1*amass(iwat(j))*(v(1)*v(1)+
     + v(2)*v(2)+v(3)*v(3))
      if(temp.gt.999.0) temp=999.989
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(y(1)+dbx),1.0d1*(y(2)+dby),
     + 1.0d1*(y(3)-dbz),temp,argos_prepare_atnam(num(iwat(j)))
      endif
      endif
   48 continue
   47 continue
      nums=nums+numslv
   46 continue
   45 continue
   44 continue
c
      endif
c
      endif
c
      write(lfnpdb,1004)
 1004 format('END')
c
c
      close(unit=lfnpdb)
c
      if(util_print('files',print_default)) then
      write(lfnout,3000) filpdb(1:length)
 3000 format(' Created pdb',t40,a,/)
      endif
c
      argos_prepare_wrtpdb=.true.
      return
c
 9999 continue
      argos_prepare_wrtpdb=.false.
      return
      end
