      logical function pre_wrtseq(lfnseq,filseq,
     + cseq,lseq,mseq,nseq,lssb,mssb,nssb,llnk,clnk,mlnk,nlnk,nlnkf,
     + slvmdl,icount,scount)
c
c $Id: pre_wrtseq.F,v 1.9 2000-02-20 16:05:15 d3j191 Exp $
c
c     function to write the sequence file
c
      implicit none
c
      integer lfnseq,mseq,nseq,mssb,nssb,mlnk,nlnk,nlnkf
      integer lseq(4,mseq),lssb(2,mssb),llnk(4,mlnk)
      character*4 clnk(2,mlnk)
      character*5 cseq(2,mseq)
      character*10 slvmdl
      character*255 filseq
      integer i,length,mol,icount,ll
      character*4 scount
c
      length=index(filseq,' ')-1
c
      open(unit=lfnseq,file=filseq(1:length),form='formatted',
     + status='new',err=9999)
c
      ll=0
      mol=0
      do 1 i=1,nseq
      if(i.gt.1.and.lseq(4,i).ne.mol) then
      write(lfnseq,1001)
 1001 format('molecule')
      endif
      mol=lseq(4,i)
      length=index(cseq(2,i),'_')-1
      if(length.le.0) length=5
c      write(lfnseq,1000) lseq(1,i),cseq(2,i)(1:length)
      write(lfnseq,1000) lseq(1,i),cseq(2,i)
 1000 format(i5,a)
      ll=lseq(1,i)
    1 continue
c
      if(icount.gt.0) then
      write(lfnseq,1005)
 1005 format('fraction')
      do 3 i=1,icount
      write(lfnseq,1006) ll+i,scount
 1006 format(i5,a)
      if(i.lt.icount) write(lfnseq,1007)
 1007 format('molecule')
    3 continue
      endif
c
c      do 2 i=1,nssb
c      write(lfnseq,1002) lssb(1,i),3,lssb(2,i),3
c    2 continue
c
      do 4 i=1,nlnk
      if(i.le.nlnkf) then
      write(lfnseq,1002) lseq(1,llnk(1,i)),clnk(1,i),
     + lseq(1,llnk(2,i)),clnk(2,i)
 1002 format('LINK ',i5,a4,i5,a4)
      else
      write(lfnseq,1003) lseq(1,llnk(1,i)),clnk(1,i),
     + lseq(1,llnk(2,i)),clnk(2,i)
 1003 format('link ',i5,a4,i5,a4)
      endif
    4 continue
c
c     solvent model
c     -------------
c
      if(slvmdl(1:1).ne.' ') then
      write(lfnseq,1004) 0,slvmdl
 1004 format('solvent',/,i5,a)
      endif
c
      write(lfnseq,1009)
 1009 format('stop')
c
      close(unit=lfnseq)
      pre_wrtseq=.true.
      return
 9999 continue
      pre_wrtseq=.false.
      return
      end
