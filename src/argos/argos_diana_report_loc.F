      subroutine argos_diana_report_loc(lfn,numi,numj,ndxi,ndxj,itxt,
     + jtxt,tim)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      integer lfn,numi,numj
      integer ndxi(numi,3),ndxj(numj,4)
      character*12 itxt(numi),jtxt(numj)
      real*8 tim(numj)
c
      integer i,j,k,nskip,it,jt(10)
      real*8 ctime,ptime,dtime,timi,timf
c
      rewind(lfn)
      read(lfn,1000) j
      read(lfn,1000) j
 1000 format(i7)
      do 1 i=1,numi
      read(lfn,1001) (ndxi(i,j),j=1,3),itxt(i)
 1001 format(i7,i5,i6,1x,a12)
    1 continue
      read(lfn,1000) j
      do 2 i=1,numj
      read(lfn,1001) (ndxj(i,j),j=1,3),jtxt(i)
    2 continue
      nskip=numi+numj+3
c
      do 3 i=1,numi
      rewind(lfn)
      do 4 j=1,nskip
      read(lfn,1000,end=3) k
    4 continue
      do 6 j=1,numj
      ndxj(j,4)=0
      tim(j)=0.0d0
    6 continue
      ptime=-1.0d0
      timi=0.0d0
    5 continue
      read(lfn,1002,end=11) ctime,it,(jt(k),k=1,10)
 1002 format(f12.6,i6,10i6)
      if(ptime.le.0.0d0) timi=ctime
      if(timf.lt.ctime) timf=ctime
      if(it.ne.ndxi(i,1)) goto 5
c
      if(ptime.ge.0.0d0) then
      dtime=ctime-ptime
      do 7 j=1,numj
      if(ndxj(j,4).ne.0) tim(j)=tim(j)+dtime
    7 continue
      endif
c
      do 8 j=1,numj
      ndxj(j,4)=0
    8 continue
      ptime=ctime
c
      do 9 k=1,10
      if(jt(k).eq.0) goto 9
      do 10 j=1,numj
      if(ndxj(j,1).eq.jt(k)) ndxj(j,4)=1
   10 continue
    9 continue
c
      goto 5
   11 continue
c
      dtime=ctime-ptime
      do 12 j=1,numj
      if(ndxj(j,4).ne.0) tim(j)=tim(j)+dtime
   12 continue
c
      do 13 j=1,numj
      if(tim(j).gt.0.0d0) then
      write(6,1003) ndxi(i,3),itxt(i),
     + ndxj(j,3),jtxt(j),tim(j),100.0d0*tim(j)/(timf-timi)
 1003 format(i7,1x,a,i7,1x,a,f12.6,' ps = ',f6.2,' %')
      endif
   13 continue
      write(6,1004)
 1004 format(' ')
c
    3 continue
c
      return
      end
c $Id$
