      logical function pre_wrtseq(lfnseq,filseq,
     + cseq,lseq,mseq,nseq,lssb,mssb,nssb,llnk,clnk,mlnk,nlnk)
c
c     function to write the sequence file
c
      implicit none
c
      integer lfnseq,mseq,nseq,mssb,nssb,mlnk,nlnk
      integer lseq(4,mseq),lssb(2,mssb),llnk(2,mlnk)
      character*4 clnk(2,mlnk)
      character*5 cseq(2,mseq)
      character*255 filseq
      integer i,j,length,mol
c
      length=index(filseq,' ')-1
c
      open(unit=lfnseq,file=filseq(1:length),form='formatted',
     + status='new',err=9999)
c
      mol=0
      do 1 i=1,nseq
      if(i.gt.1.and.lseq(4,i).ne.mol) then
      write(lfnseq,1001)
 1001 format('molecule')
      endif
      mol=lseq(4,i)
      length=index(cseq(2,i),'_')-1
      if(length.le.0) length=5
      write(lfnseq,1000) lseq(1,i),cseq(2,i)(1:length)
 1000 format(i5,a)
    1 continue
c
c      do 2 i=1,nssb
c      write(lfnseq,1002) lssb(1,i),3,lssb(2,i),3
c    2 continue
c
      do 3 i=1,nlnk
      write(lfnseq,1002) lseq(1,llnk(1,i)),clnk(1,i),
     + lseq(1,llnk(2,i)),clnk(2,i)
 1002 format('link ',i5,a4,i5,a4)
    3 continue
c
      write(lfnseq,1003)
 1003 format('stop')
c
      close(unit=lfnseq)
      pre_wrtseq=.true.
      return
 9999 continue
      pre_wrtseq=.false.
      return
      end
