      logical function pre_gensgm(lfnout,lprint,lfnfrg,lfnsgm,
     + ffield,dir_s,dir_x,dir_u,dir_t,lsgm,csgm,msgm,nsgm)
c
c     function to generate segment files
c
c     in  : lfnfrg    = fragment file logical file number
c           lfnsgm    = segment file logical file number
c           dir_s     = dbase directory name
c           dir_x     = dbase directory name
c           dir_u     = dbase directory name
c           dir_t     = dbase directory name
c           lfnout    = output file logical file number
c           lprint    = integer for print to output
c           mseq      = dimension of the sequence list
c           lsgm(1,i) = number of segments of type i
c           lsgm(2,i) = source: 0=not found; 1=s; 2=x; 3=u; 4=t;
c           csgm(i)   = segment name
c
      implicit none
c
#include "mafdecls.fh"
c
      logical pre_mksgm
      external pre_mksgm
c
      integer lfnout,lfnfrg,lfnsgm
      integer lprint
      integer msgm,nsgm
      integer lsgm(3,msgm)
      character*5 csgm(msgm)
      character*80 ffield,dir_s,dir_x,dir_u,dir_t
c
      character*80 filfrg,filsgm
c
      integer len_s,len_x,len_u,len_t,length
      integer i
c
      if(lprint.ge.5) then
      write(lfnout,2000)
 2000 format(/,'pre_gensgm ')
      endif
c
      len_s=index(dir_s,' ')-1
      len_x=index(dir_x,' ')-1
      len_u=index(dir_u,' ')-1
      len_t=index(dir_t,' ')-1
c
      do 1 i=1,nsgm
      if(lsgm(2,i).lt.0) then
c
      length=index(csgm(i),' ')-1
      if(length.le.0) length=5
c
      if(lsgm(2,i).eq.-1) then
      filfrg=dir_s(1:len_s)//csgm(i)(1:length)//'.frg '
      filsgm=dir_t(1:len_t)//csgm(i)(1:length)//'.sgm '
      else
      if(lsgm(2,i).eq.-2) then
      filfrg=dir_x(1:len_x)//csgm(i)(1:length)//'.frg '
      filsgm=dir_t(1:len_t)//csgm(i)(1:length)//'.sgm '
      else
      if(lsgm(2,i).eq.-3) then
      filfrg=dir_u(1:len_u)//csgm(i)(1:length)//'.frg '
      filsgm=dir_t(1:len_t)//csgm(i)(1:length)//'.sgm '
      else
      if(lsgm(2,i).eq.-4) then
      filfrg=dir_t(1:len_t)//csgm(i)(1:length)//'.frg '
      filsgm=dir_t(1:len_t)//csgm(i)(1:length)//'.sgm '
      else
      goto 9999
      endif
      endif
      endif
      endif
c
      if(.not.pre_mksgm(lfnout,lprint,ffield,
     + lfnfrg,filfrg,lfnsgm,filsgm))
     + call errquit('pre_mksgm failed',9999)
c
      endif
    1 continue
c
      pre_gensgm=.true.
      return
c
 9999 continue
      pre_gensgm=.false.
      return
      end
