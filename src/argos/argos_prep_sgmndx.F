      logical function argos_prep_sgmndx(iunit,
     + lfnout,lseq,cseq,mseq,nseq,lsgm,csgm,msgm,nsgm,mato)
c
c $Id$
c
c     function to scan database directories for segments in sequence
c
c     in  : iunit        = dbase logical file number
c           lfnout       = output file logical file number
c           lseq(1,mseq) = segment number
c           lseq(2,mseq) = index to segment list
c           cseq(2,mseq) = segment name on topology
c           mseq         = dimension of the sequence list
c           nseq         = length of sequence list
c           msgm         = dimension of segment list
c           msgm         = length of segment list
c
c     out : lseq(1,mseq) = segment numbers
c                2       = number of atoms
c                3       = index to unique segment
c           csgm(msgm)   = unique segment names
c           lsgm(1,msgm) = number of segments of type i
c                2       = source: 0=not found; sgm: 1=s; 2=x; 3=u; 4=t;
c                                               frg:-1=s;-2=x;-3=u;-4=t;
c                3       = number of atoms in segment
c
      implicit none
c
#include "util.fh"
#include "argos_prep_common.fh"
c
      integer iunit,lfnout,mseq,nseq,msgm,nsgm,mato
      integer lseq(6,mseq),lsgm(3,msgm)
      character*5 cseq(2,mseq),csgm(msgm)
      character*255 filnam
c
      integer length
      integer i,j
c
      if(util_print('sequence',print_high)) then
      write(lfnout,2000)
 2000 format(/,' Segment definition files ',/)
      endif
c
      argos_prep_sgmndx=.true.
c
c     initialize  number of segments to zero
c     ----------  source of segments to unknown
c
      do 1 i=1,nsgm
      lsgm(1,i)=0
      lsgm(2,i)=0
      lsgm(3,i)=0
    1 continue
c
c     find all unique segments in sequence csgm(1:nsgm) = segment names
c     ------------------------------------ lsgm(1,1:nsgm)
c
      mato=0
      do 2 i=1,nseq
      mato=max(mato,lseq(3,i+1)-lseq(3,i))
      csgm(lseq(2,i))=cseq(2,i)
      lsgm(1,lseq(2,i))=lsgm(1,lseq(2,i))+1
      if(lsgm(1,lseq(2,i)).eq.1) then
      csgm(lseq(2,i))=cseq(2,i)
      lsgm(3,lseq(2,i))=lseq(3,i+1)-lseq(3,i)
      endif
    2 continue
c
c     find segment files for the segments in the list
c     -----------------------------------------------
c
      do 3 i=1,nsgm
      length=index(csgm(i),' ')-1
      if(length.le.0) length=5
      lsgm(2,i)=0
c
c     check the temporary dbase directory
c
      do 4 j=mdirpar,1,-1
      if(nfilpar(j).gt.0) then
      filnam=
     + dirpar(j)(1:index(dirpar(j),' ')-1)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=5)
      lsgm(2,i)=j
      close(iunit)
      goto 6
    5 continue
      filnam=
     + dirpar(j)(1:index(dirpar(j),' ')-1)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=4)
      lsgm(2,i)=-j
      close(iunit)
      goto 6
      endif
    4 continue
    6 continue
c
      if(lsgm(2,i).eq.0) argos_prep_sgmndx=.false.
c
      if(util_print('sequence',print_high)) then
      if(lsgm(2,i).lt.0)
     + write(lfnout,1001) i,csgm(i),lsgm(1,i),lsgm(2,i),
     + dirpar(-lsgm(2,i))(1:index(dirpar(-lsgm(2,i)),' ')-1)//
     + csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).gt.0)
     + write(lfnout,1001) i,csgm(i),lsgm(1,i),lsgm(2,i),
     + dirpar(lsgm(2,i))(1:index(dirpar(lsgm(2,i)),' ')-1)//
     + csgm(i)(1:length)//'.sgm '
 1001 format(i5,2x,a5,i5,2x,i5,1x,a)
      if(lsgm(2,i).eq.0)write(lfnout,1002) i,csgm(i),lsgm(1,i),
     + lsgm(2,i)
 1002 format(i5,2x,a5,i5,2x,i5,1x,' no definition file found')
      endif
c
    3 continue
c
      return
      end
