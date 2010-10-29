      logical function argos_prepare_sgmnum(lfnout,
     + lseq,cseq,mseq,nseq,nsgm)
c
c $Id$
c
c     function to scan database for unique segments in sequence
c
c     in  : filnam    = dbase file name
c           lfnout    = output file logical file number
c           lseq(1,*) = sequence segment numbers
c           cseq(2,*) = sequence segment names
c           mseq      = dimension of the sequence list
c           nseq      = length of sequence list
c
c     out : lseq(2,*) = index to segment list
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,mseq,nseq,nsgm
      integer lseq(6,mseq)
      character*5 cseq(2,mseq)
c
      integer i,j
c
      if(util_print('sequence',print_high)) then
      write(lfnout,2000)
 2000 format(/,' Unique segments  ',/)
      endif
c
      nsgm=0
      do 1 i=1,nseq
      do 2 j=1,i-1
      if(cseq(2,i).eq.cseq(2,j)) then
      lseq(2,i)=lseq(2,j)
      goto 1
      endif
    2 continue
      nsgm=nsgm+1
      lseq(2,i)=nsgm
      if(util_print('sequence',print_high)) then
      write(lfnout,2001) lseq(2,i),cseq(2,i)
 2001 format(i5,2x,a)
      endif
    1 continue
c
      argos_prepare_sgmnum=.true.
      return
      end
