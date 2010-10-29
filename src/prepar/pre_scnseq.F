      logical function pre_scnseq(iunit,dirnam,lfnout,
     + kseq,cseq,mseq,ksgm,csgm,msgm)
c
c $Id$
c
c     function to scan database for segments in sequence
c
c     return value:
c          .false. : scan failed
c          .true.  : scan completed successfully
c
c     in  : iunit     = dbase logical file number
c           filnam    = dbase file name
c           lfnout    = output file logical file number
c           mseq      = dimension of the sequence list
c     out : kseq(1,*) = segment numbers
c           kseq(2,*) = number of atoms
c           kseq(3,*) = index to unique segment
c           cseq      = segment names
c
      implicit none
c
#include "util.fh"
c
      integer iunit,lfnout,mseq,msgm
      integer kseq(3,mseq),ksgm(3,msgm)
      character*5 cseq(2,mseq),csgm(msgm)
      character*255 dirnam,filnam
c
      character*80 card
      integer length,lseq,isgm,i
      logical found
c
      lseq=mseq
c
      length=index(dirnam,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2006) dirnam(1:length)
 2006 format(/,'pre_scnseq  ',a)
      endif
c
c     initialize the lists
c
      do 9 i=1,msgm
      ksgm(1,i)=0
    9 continue
c
      do 10 i=1,lseq
c      write(lfnout,1000) i,kseq(1,i),cseq(1,i),kseq(2,i),kseq(3,i)
c 1000 format(2i5,2x,a5,2i5)
      ksgm(1,kseq(3,i))=ksgm(1,kseq(3,i))+1
      if(ksgm(1,kseq(3,i)).eq.1) csgm(kseq(3,i))=cseq(1,i)
   10 continue
c
      do 11 i=1,msgm
      filnam=dirnam(1:length)//csgm(i)(1:index(csgm(i),' ')-1)//'.frg'
      found=.false.
      open(iunit,file=filnam,form='formatted',status='old',err=19)
      found=.true.
      close(iunit)
   19 continue
c      if(found) then
c      write(lfnout,1001) i,csgm(i),ksgm(1,i),
c     + filnam(1:index(filnam,' ')-1)
c 1001 format(i5,2x,a5,i5,2x,a)
c      else
c      write(lfnout,1002) i,csgm(i),ksgm(1,i)
c 1002 format(i5,2x,a5,i5)
c      endif
   11 continue
c
      pre_scnseq=.true.
      return
      end
