      logical function pre_sgmnum(lfnout,
     + lseq,cseq,mseq,nseq,nsgm)
c
c $Id: pre_sgmnum.F,v 1.5 1998-06-08 23:29:38 d3j191 Exp $
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
      integer lseq(4,mseq)
      character*5 cseq(2,mseq)
c
      integer i,j
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000)
 2000 format(/,'pre_sgmnum  ')
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
    1 continue
c
      pre_sgmnum=.true.
      return
      end
