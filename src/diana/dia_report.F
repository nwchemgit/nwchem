      subroutine dia_report_file(opt,fil,lfn)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "mafdecls.fh"
c
      character*8 opt
      character*255 fil
      integer lfn
c
      integer i,j,num,numi,numj
      integer i_ii,l_ii,i_ic,l_ic
      integer i_ji,l_ji,i_jc,l_jc,i_jr,l_jr
c
      if(me.ne.0) return
c
      if(opt.eq.'local   '.and.index(fil,'.loc').gt.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1))
      rewind(lfn)
      read(lfn,1000) num
 1000 format(i7)
      read(lfn,1000) numi
      do 1 i=1,numi
      read(lfn,1000) j
    1 continue
      read(lfn,1000) numj
      if(.not.ma_push_get(mt_int,3*numi,'ii',l_ii,i_ii))
     + call md_abort('Could not allocate ii',0)
      if(.not.ma_push_get(mt_byte,12*numi,'ic',l_ic,i_ic))
     + call md_abort('Could not allocate ic',0)
      if(.not.ma_push_get(mt_int,4*numj,'ji',l_ji,i_ji))
     + call md_abort('Could not allocate ji',0)
      if(.not.ma_push_get(mt_byte,12*numj,'jc',l_jc,i_jc))
     + call md_abort('Could not allocate jc',0)
      if(.not.ma_push_get(mt_dbl,numj,'jr',l_jr,i_jr))
     + call md_abort('Could not allocate jr',0)
      call dia_report_loc(lfn,numi,numj,int_mb(i_ii),int_mb(i_ji),
     + byte_mb(i_ic),byte_mb(i_jc),dbl_mb(i_jr))
      if(.not.ma_pop_stack(l_jr))
     + call md_abort('Failed to deallocate jr',0)
      if(.not.ma_pop_stack(l_jc))
     + call md_abort('Failed to deallocate jc',0)
      if(.not.ma_pop_stack(l_ji))
     + call md_abort('Failed to deallocate ji',0)
      if(.not.ma_pop_stack(l_ic))
     + call md_abort('Failed to deallocate ic',0)
      if(.not.ma_pop_stack(l_ii))
     + call md_abort('Failed to deallocate ii',0)
      close(unit=lfn)
      endif
c
      return
      end
      subroutine dia_report_loc(lfn,numi,numj,ndxi,ndxj,itxt,jtxt,tim)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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
