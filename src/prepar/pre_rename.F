      logical function pre_rename(lfnout,lprint,ffield,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c $Id: pre_rename.F,v 1.1 1998-06-03 14:36:37 d3j191 Exp $
c
c     function to comply with ffield segment names
c
      implicit none
c
      integer lfnout
      integer lprint
      character*80 ffield
      integer mseq,nseq,matm,natm
      integer lseq(4,mseq),latm(5,matm)
      character*5 cseq(2,mseq)
      character*6 catm(3,matm)
c
      integer i
c
      do 1 i=1,nseq
      cseq(2,i)=cseq(1,i)
    1 continue
c
      if(ffield(1:5).eq.'amber') then
      endif
c
      pre_rename=.true.
      return
c
 9999 continue
      pre_rename=.false.
      return
      end
