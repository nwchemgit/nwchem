      logical function pre_term(lfnout,lfnsgm,
     + lseq,cseq,mseq,nseq,llnk,clnk,mlnk,nlnk)
c
c $Id$
c
c     function to determine the terminal segments
c
      implicit none
c
#include "util.fh"
c
      logical pre_sgmfil,pre_atmscn
      external pre_sgmfil,pre_atmscn
c
      integer lfnout,lfnsgm,mseq,nseq,mlnk,nlnk
      integer lseq(6,mseq),llnk(3,mlnk)
      character*10 cseq(mseq)
      character*4 clnk(2,mlnk)
c
      character*255 filsgm
      integer i,isgm,j,jsgm,num,link
c
c     process the explicit link list
c     ------------------------------
c
      do 1 i=1,nlnk
c
      isgm=0
      jsgm=0
      do 2 j=1,nseq
      if(lseq(1,j).eq.llnk(1,i)) isgm=j
      if(lseq(1,j).eq.llnk(2,i)) jsgm=j
    2 continue
      if(isgm.eq.0.or.jsgm.eq.0) then
      write(lfnout,1003) llnk(1,i),llnk(2,i)
 1003 format(' Error in link ',2i5)
      call md_abort('error in pre_term',9999)
      endif
c
      if(util_print('where',print_high)) then
      write(lfnout,1000) llnk(1,i),clnk(1,i),llnk(2,i),clnk(2,i)
 1000 format(' Link ',i5,':',a4,4x,i5,':',a4)
      endif
c
      if(.not.pre_sgmfil(cseq(isgm),lfnsgm,filsgm)) then
      write(lfnout,1001) filsgm(1:index(filsgm,' ')-1)
 1001 format('Segment ',a,' not found')
      call md_abort('pre_sgmfil failed',9999)
      endif
c
      if(.not.pre_atmscn(lfnsgm,filsgm,clnk(1,i),num,link)) then
      write(lfnout,1002) clnk(1,i),filsgm(1:index(filsgm,' ')-1)
 1002 format('Atom ',a4,'not found in segment ',a)
      call md_abort('pre_atmscn failed',9999)
      endif
c
      if(link.eq.1) then
cx      if(lseq(2,isgm).ge.0)  call md_abort('multiple links',9999)
      if(lseq(2,isgm).eq.-1) lseq(2,isgm)=0
      if(lseq(2,isgm).eq.-2) lseq(2,isgm)=1
      endif
      if(link.eq.2) then
cx      if(lseq(2,isgm).eq.-1) call md_abort('multiple links',9999)
cx      if(lseq(2,isgm).eq.0)  call md_abort('multiple links',9999)
      if(lseq(2,isgm).eq.-2) lseq(2,isgm)=-1
      if(lseq(2,isgm).gt.0) lseq(2,isgm)=0
      endif
c
      if(.not.pre_sgmfil(cseq(jsgm),lfnsgm,filsgm)) then
      write(lfnout,1001) filsgm(1:index(filsgm,' ')-1)
      call md_abort('pre_sgmfil failed',9999)
      endif
c
      if(.not.pre_atmscn(lfnsgm,filsgm,clnk(2,i),num,link)) then
      write(lfnout,1001) filsgm(1:index(filsgm,' ')-1)
      write(lfnout,1002) clnk(2,i),filsgm(1:index(filsgm,' ')-1)
      call md_abort('pre_atmscn failed',9999)
      endif
c
      if(link.eq.1) then
cx      if(lseq(2,jsgm).ge.0)  call md_abort('multiple links',9999)
      if(lseq(2,jsgm).eq.-1) lseq(2,jsgm)=0
      if(lseq(2,jsgm).eq.-2) lseq(2,jsgm)=1
      endif
      if(link.eq.2) then
cx      if(lseq(2,jsgm).eq.-1) call md_abort('multiple links',9999)
cx      if(lseq(2,jsgm).eq.0)  call md_abort('multiple links',9999)
      if(lseq(2,jsgm).eq.-2) lseq(2,jsgm)=-1
      if(lseq(2,jsgm).gt.0) lseq(2,jsgm)=0
      endif
c
    1 continue
c
      pre_term=.true.
      return
c
 9999 continue
      pre_term=.false.
      return
      end
