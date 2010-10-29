      logical function argos_prep_wrtseq(lfnseq,filseq,lfnout,
     + cseq,lseq,mseq,nseq,lssb,mssb,nssb,llnk,clnk,mlnk,nlnk,nlnkf,
     + slvmdl,mcount,ncount,icount,scount,mfract,nfract,ifract)
c
c $Id$
c
c     function to write the sequence file
c
      implicit none
c
#include "util.fh"
c
      integer lfnseq,mseq,nseq,mssb,nssb,mlnk,nlnk,nlnkf,lfnout
      integer lseq(6,mseq),lssb(2,mssb),llnk(4,mlnk)
      character*4 clnk(2,mlnk)
      character*5 cseq(2,mseq)
      character*10 slvmdl
      character*255 filseq
      integer mcount,ncount
      integer i,j,l,length,mol,molec,icount(mcount),ll
      character*4 scount(mcount)
      integer mfract,nfract,ifract(mfract)
c
      length=index(filseq,' ')-1
c
      open(unit=lfnseq,file=filseq(1:length),form='formatted',
     + status='new',err=9999)
c
      ll=0
      mol=0
      molec=1
      do 1 i=1,nseq
      if(i.gt.1.and.lseq(4,i).ne.mol) then
      if(nfract.eq.0) then
      write(lfnseq,1001)
 1001 format('molecule')
      else
      do 11 j=1,nfract
      if(ifract(j).eq.molec) then
      write(lfnseq,1005)
      goto 111
      endif
   11 continue
      write(lfnseq,1001)
  111 continue
      endif
      molec=molec+1
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
      if(ncount.gt.0) then
      write(lfnseq,1005)
 1005 format('fraction')
      l=0
      do 3 i=1,ncount
      do 33 j=1,icount(i)
      l=l+1
      write(lfnseq,1006) ll+l,scount(i)
 1006 format(i5,a)
      if(i.lt.ncount.or.j.lt.icount(i)) write(lfnseq,1007)
 1007 format('molecule')
   33 continue
      if(util_print('topology',print_medium)) then
      write(lfnout,1060) icount(i),scount(i)
 1060 format(' Appended ',i5,' counter ions ',a,' to topology')
      endif
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
      argos_prep_wrtseq=.true.
      return
 9999 continue
      argos_prep_wrtseq=.false.
      return
      end
