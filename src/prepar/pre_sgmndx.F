      logical function pre_sgmndx(iunit,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + lfnout,lseq,cseq,mseq,nseq,lsgm,csgm,msgm,nsgm,mato)
c
c $Id: pre_sgmndx.F,v 1.11 2001-11-06 06:15:21 d3j191 Exp $
c
c     function to scan database directories for segments in sequence
c
c     in  : iunit        = dbase logical file number
c           dir_s        = standard dbase directory
c           dir_x        = extensions dbase directory
c           dir_u        = user preferences dbase directory
c           dir_t        = temporary dbase directory
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
c
      integer iunit,lfnout,mseq,nseq,msgm,nsgm,mato
      integer lseq(4,mseq),lsgm(3,msgm)
      character*5 cseq(2,mseq),csgm(msgm)
      character*255 dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,filnam
c
      integer len_s,len_x,len_q,len_u,len_t,len_c,length
      integer i
c
      len_s=index(dir_s,' ')-1
      len_x=index(dir_x,' ')-1
      len_q=index(dir_q,' ')-1
      len_u=index(dir_u,' ')-1
      len_t=index(dir_t,' ')-1
      len_c=index(dir_c,' ')-1
c
      if(util_print('sequence',print_high)) then
      write(lfnout,2000)
 2000 format(/,' Segment definition files ',/)
      endif
c
      pre_sgmndx=.true.
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
      if(len_c.gt.0) then
      filnam=dir_c(1:len_c)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=4)
      lsgm(2,i)=6
      close(iunit)
      goto 15
    4 continue
      filnam=dir_c(1:len_c)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=5)
      lsgm(2,i)=-6
      close(iunit)
      goto 15
      endif
    5 continue
c
      if(len_t.gt.0) then
      filnam=dir_t(1:len_t)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=6)
      lsgm(2,i)=5
      close(iunit)
      goto 15
    6 continue
      filnam=dir_t(1:len_t)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=7)
      lsgm(2,i)=-5
      close(iunit)
      goto 15
      endif
    7 continue
c
c     check the user preferences dbase directory
c
      if(len_u.gt.0) then
      filnam=dir_u(1:len_u)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=8)
      lsgm(2,i)=4
      close(iunit)
      goto 15
    8 continue
      filnam=dir_u(1:len_u)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=9)
      lsgm(2,i)=-4
      close(iunit)
      goto 15
      endif
    9 continue
c
c
      if(len_q.gt.0) then
      filnam=dir_q(1:len_q)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=10)
      lsgm(2,i)=3
      close(iunit)
      goto 15
   10 continue
      filnam=dir_u(1:len_u)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=11)
      lsgm(2,i)=-3
      close(iunit)
      goto 15
      endif
   11 continue
c
c     check the extensions dbase directory
c
      if(len_x.gt.0) then
      filnam=dir_x(1:len_x)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=12)
      lsgm(2,i)=2
      close(iunit)
      goto 15
   12 continue
      filnam=dir_x(1:len_x)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=13)
      lsgm(2,i)=-2
      close(iunit)
      goto 15
      endif
   13 continue
c
c     check the standard dbase directory
c
      if(len_s.gt.0) then
      filnam=dir_s(1:len_s)//csgm(i)(1:length)//'.sgm '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=14)
      lsgm(2,i)=1
      close(iunit)
      goto 15
   14 continue
      filnam=dir_s(1:len_s)//csgm(i)(1:length)//'.frg '
      open(iunit,file=filnam(1:index(filnam,' ')-1),form='formatted',
     + status='old',err=15)
      lsgm(2,i)=-1
      close(iunit)
      goto 15
      endif
c
   15 continue
c
      if(lsgm(2,i).eq.0) pre_sgmndx=.false.
c
      if(util_print('sequence',print_high)) then
      if(lsgm(2,i).eq.-1) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_s(1:len_s)//csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).eq.1) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_s(1:len_s)//csgm(i)(1:length)//'.sgm '
      if(lsgm(2,i).eq.-2) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_x(1:len_x)//csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).eq.2) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_x(1:len_x)//csgm(i)(1:length)//'.sgm '
      if(lsgm(2,i).eq.-3) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_q(1:len_q)//csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).eq.3) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_q(1:len_q)//csgm(i)(1:length)//'.sgm '
      if(lsgm(2,i).eq.-4) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_u(1:len_u)//csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).eq.4) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_u(1:len_u)//csgm(i)(1:length)//'.sgm '
      if(lsgm(2,i).eq.-5) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_t(1:len_t)//csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).eq.5) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_t(1:len_t)//csgm(i)(1:length)//'.sgm '
      if(lsgm(2,i).eq.-6) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_c(1:len_c)//csgm(i)(1:length)//'.frg '
      if(lsgm(2,i).eq.5) write(lfnout,1001) i,csgm(i),lsgm(1,i),
     + lsgm(2,i),dir_c(1:len_c)//csgm(i)(1:length)//'.sgm '
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
