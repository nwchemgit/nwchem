      subroutine argos_diana_siztop()
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "argos_diana_common.fh"
c
      character*1 cdummy
      real*8 rdummy
c
      integer i,j,k,nat,naq,nparm,nseq
      integer naw,nbw,nhw,ndw,now,ntw,nnw
      integer nas,nbs,num
      character*5 sname,aname
      character*80 string
      integer ibroad(10)
c
      ltop=.false.
      nsgm=0
c
      if(me.eq.0) then
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=999)
c
      read(lfntop,1000,end=999,err=999) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
 1000 format(a1)
 5000 format(a)
      read(lfntop,1001) nparm
      read(lfntop,1001) nat
      read(lfntop,1001) naq
      read(lfntop,1001) nseq
 1001 format(i5)
      read(lfntop,1000) cdummy
      do 1 i=1,nat*nparm
      read(lfntop,1000) cdummy
    1 continue
      do 2 i=1,nat
      do 3 j=i,nat
      do 4 k=1,nparm
      read(lfntop,1000) cdummy
    4 continue
    3 continue
    2 continue
      do 5 i=1,naq*nparm
      read(lfntop,1002) rdummy
 1002 format(5x,f12.6)
    5 continue
      do 6 i=1,nseq
      read(lfntop,1000) cdummy
    6 continue
      read(lfntop,1003) naw,nbw,nhw,ndw,now,ntw,nnw
 1003 format(5i7,2i10)
      mwb=nbw
      read(lfntop,1003) nas,nbs
      msb=nbs
      do 7 i=1,naw
      read(lfntop,1005) sname,aname,num,j
 2004 format(2a5,i6)
    7 continue
      nat=(nparm+1)*(nbw+nhw+ndw+now)
      do 8 i=1,nat
      read(lfntop,1000) cdummy
    8 continue
      if(ntw.gt.0) then
      read(lfntop,1004) (j,i=1,ntw)
      read(lfntop,1004) (j,i=1,ntw)
 1004 format(11i7)
      endif
      if(nnw.gt.0) then
      read(lfntop,1004) (j,i=1,nnw)
      read(lfntop,1004) (j,i=1,nnw)
      endif
      read(lfntop,1000) cdummy
      do 9 i=1,nparm
      read(lfntop,1000) cdummy
    9 continue
      do 10  i=1,nas
      read(lfntop,1005) sname,aname,num,j
 1005 format(a5,5x,a5,11x,i7,14x,i5)
c 1005 format(a5,5x,a5,11x,i5,15x,i5)
      nsgm=max(nsgm,num)
   10 continue
      do 11 i=1,nbs
      read(lfntop,1006) j,num
      read(lfntop,1000) cdummy
 1006 format(2i7)
 2005 format(2i8)
   11 continue
      close(unit=lfntop)
c
      msgm=nsgm
      ltop=.true.
      endif
  999 continue
c
      if(np.gt.1) then
      ibroad(1)=nsgm
      ibroad(2)=msgm
      call ga_brdcst(mag_d02,ibroad,2*ma_sizeof(mt_int,1,mt_byte),0)
      nsgm=ibroad(1)
      msgm=ibroad(2)
      endif
c
      return
      end
c $Id$
