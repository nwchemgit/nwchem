      logical function pre_wrtpdb(lfnout,lfnpdb,filpdb,iopt,box,
     + num,amass,mat,
     + csa,isat,isgm,xs,vs,msa,nsa,cwa,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc,slvnam,nxrep,nyrep,nzrep,drep,msb,nsb,idsb,
     + rdist)
c
c $Id: pre_wrtpdb.F,v 1.17 2000-05-03 18:03:58 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
      character*2 pre_atnam
      real*8 pre_atrad
      external pre_atnam,pre_atrad
c
      integer lfnout,lfnpdb,iopt,msa,nsa,mwm,nwm,mwa,nwa,mat,mwmc,nwmc
      character*255 filpdb
      integer num(mat),isat(msa),isgm(msa)
      character*16 cwa(mwa),csa(msa)
      character*3 slvnam
      real*8 amass(mat),xs(3,msa),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),box(3)
      integer nxrep,nyrep,nzrep
      real*8 rdist,drep(3)
      integer msb,nsb
      integer idsb(2,msb)
c
      integer length
c
      integer i,j,k,nats,nums,numslv,isoff,iaoff
      integer ixrep,iyrep,izrep
      real*8 temp,dbx,dby,dbz,dist,dbz0
c
      isoff=0
      iaoff=0
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
      write(lfnpdb,1000) (1.0d1*box(i),i=1,3),90.0,90.0,90.0
 1000 format(
     + 'HEADER',/,
     + 'TITLE     ',/,
     + 'TITLE    2',/,
     + 'TITLE    3',/,
     + 'REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.1, 25-OCT-1996',/,
     + 'CRYST1',3f9.3,3f7.2)
      nums=0
      do 11 k=1,nsb
      i=idsb(1,k)
      j=idsb(2,k)
      if(isgm(i).ne.isgm(j)) then
      dist=sqrt((xs(1,j)-xs(1,i))*(xs(1,j)-xs(1,i))+
     + (xs(2,j)-xs(2,i))*(xs(2,j)-xs(2,i))+
     + (xs(3,j)-xs(3,i))*(xs(3,j)-xs(3,i)))
      if(dist.gt.pre_atrad(num(isat(j)))+pre_atrad(num(isat(i)))) then
      write(lfnpdb,1101) csa(i)(11:14),csa(i)(1:3),isgm(i),
     + csa(j)(11:14),csa(j)(1:3),isgm(j)
 1101 format('LINK',t13,a4,1x,a3,2x,i4,16x,a4,1x,a3,2x,i4)
      endif
      endif
   11 continue
      do 1 i=1,nsa
      temp=4.009103873d1*amass(isat(i))*
     + (vs(1,i)*vs(1,i)+vs(2,i)*vs(2,i)+vs(3,i)*vs(3,i))
      if(temp.gt.999.0) temp=999.989
      write(lfnpdb,1001) i,csa(i)(11:14),csa(i)(1:3),
     + isgm(i),(1.0d1*xs(k,i),k=1,3),temp,pre_atnam(num(isat(i)))
 1001 format('ATOM',i7,1x,a4,1x,a3,2x,i4,4x,3f8.3,6x,f6.2,4x,a2)
      if(isgm(i).gt.nums) nums=isgm(i)
    1 continue
      else
      dbz0=0.0d0
      if(nzrep.lt.0) then
      if(nsa.gt.0.and.nwm.eq.0) dbz0=xs(3,1)
      if(nsa.eq.0.and.nwm.gt.0) dbz0=xw(3,1,1)
      if(nsa.gt.0.and.nwm.gt.0) dbz=min(xs(3,1),xw(3,1,1))
      do 31 i=1,nsa
      dbz0=min(dbz0,xs(3,i))
   31 continue
      do 32 i=1,nwm
      do 33 j=1,nwa
      dbz0=min(dbz0,xw(3,j,i))
   33 continue
   32 continue
      dbz0=dbz0-0.5d0*rdist
      endif
      write(lfnpdb,1000) 1.0d1*dble(nxrep)*drep(1),
     + 1.0d1*dble(nyrep)*drep(2),1.0d1*dble(nzrep)*drep(3),
     + 90.0,90.0,90.0
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
      do 111 k=1,nsb
      i=idsb(1,k)
      j=idsb(2,k)
      if(isgm(i).ne.isgm(j)) then
      dist=sqrt((xs(1,j)-xs(1,i))*(xs(1,j)-xs(1,i))+
     + (xs(2,j)-xs(2,i))*(xs(2,j)-xs(2,i))+
     + (xs(3,j)-xs(3,i))*(xs(3,j)-xs(3,i)))
      if(dist.gt.pre_atrad(num(isat(j)))+pre_atrad(num(isat(i)))) then
      write(lfnpdb,1101) csa(i)(11:14),csa(i)(1:3),isgm(i)+isoff,
     + csa(j)(11:14),csa(j)(1:3),isgm(j)+isoff
      endif
      endif
  111 continue
c
      nums=0
c
      do 7 i=1,nsa
      temp=4.009103873d1*amass(isat(i))*
     + (vs(1,i)*vs(1,i)+vs(2,i)*vs(2,i)+vs(3,i)*vs(3,i))
      if(temp.gt.999.0) temp=999.989
      if(nzrep.gt.0) then
      write(lfnpdb,1001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+isoff,1.0d1*(xs(1,i)+dbx),
     + 1.0d1*(xs(2,i)+dby),1.0d1*(xs(3,i)+dbz),
     + temp,pre_atnam(num(isat(i)))
      else
      if(izrep.eq.1) then
      write(lfnpdb,1001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+isoff,1.0d1*(xs(1,i)+dbx),
     + 1.0d1*(xs(2,i)+dby),1.0d1*(xs(3,i)-dbz),
     + temp,pre_atnam(num(isat(i)))
      else
      write(lfnpdb,1001) i+iaoff,csa(i)(11:14),csa(i)(1:3),
     + isgm(i)+isoff,-1.0d1*(xs(1,i)+dbx),
     + -1.0d1*(xs(2,i)+dby),-1.0d1*(xs(3,i)-dbz),
     + temp,pre_atnam(num(isat(i)))
      endif
      endif
      if(isgm(i).gt.nums) nums=isgm(i)
    7 continue
      isoff=isoff+nums
      iaoff=iaoff+nsa
    6 continue
    5 continue
    4 continue
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
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,(1.0d1*xwc(k,j,i),k=1,3)
      else
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,(1.0d1*xw(k,j,i-nwmc),k=1,3)
      endif
 1003 format('ATOM',i7,1x,a4,1x,a3,i6,4x,3f8.3)
    3 continue
    2 continue
c
      else
c
      nats=nxrep*nyrep*iabs(nzrep)*nsa
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
      do 47 i=1,numslv
      do 48 j=1,nwa
      nats=nats+1
      if(nzrep.gt.0) then
      if(i.le.nwmc) then
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(xwc(1,j,i)+dbx),1.0d1*(xwc(2,j,i)+dby),
     + 1.0d1*(xwc(3,j,i)+dbz)
      else
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(xw(1,j,i-nwmc)+dbx),1.0d1*(xw(2,j,i-nwmc)+dby),
     + 1.0d1*(xw(3,j,i-nwmc)+dbz)
      endif
      else
      if(izrep.eq.1) then
      if(i.le.nwmc) then
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(xwc(1,j,i)+dbx),1.0d1*(xwc(2,j,i)+dby),
     + 1.0d1*(xwc(3,j,i)-dbz)
      else
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(xw(1,j,i-nwmc)+dbx),1.0d1*(xw(2,j,i-nwmc)+dby),
     + 1.0d1*(xw(3,j,i-nwmc)-dbz)
      endif
      else
      if(i.le.nwmc) then
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,1.0d1*(xwc(1,j,i)+dbx),1.0d1*(xwc(2,j,i)+dby),
     + 1.0d1*(xwc(3,j,i)-dbz)
      else
      write(lfnpdb,1003) nats,cwa(j)(11:14),slvnam,
     + i+nums,-1.0d1*(xw(1,j,i-nwmc)+dbx),-1.0d1*(xw(2,j,i-nwmc)+dby),
     + -1.0d1*(xw(3,j,i-nwmc)-dbz)
      endif
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
      write(lfnout,2000) filpdb(1:length)
 2000 format(' Created pdb',t40,a,/)
      endif
c
      pre_wrtpdb=.true.
      return
c
 9999 continue
      pre_wrtpdb=.false.
      return
      end
