      subroutine ana_rdtop(qs)
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "ana_common.fh"
c
      real*8 qs(msa)
c
      character*1 cdummy
      integer i,j,nat,naq,i_tmp,l_tmp
      integer naw,nbw,nhw,ndw,now,ntw,nnw
      integer nas
c
      if(me.eq.0) then
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old')
c
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
 1000 format(a1)
      read(lfntop,1001) nat
      read(lfntop,1001) naq
 1001 format(i5)
      do 1 i=1,nat
      read(lfntop,1000) cdummy
    1 continue
      do 2 i=1,nat
      do 3 j=i,nat
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    3 continue
    2 continue
      if(.not.ma_push_get(mt_dbl,naq,'tmp',l_tmp,i_tmp))
     + call errquit('Failed to allocate tmp',0)
      do 4 i=1,naq
      read(lfntop,1002) dbl_mb(i_tmp-1+i)
 1002 format(f12.6)
    4 continue
      read(lfntop,1003) naw,nbw,nhw,ndw,now,ntw,nnw
 1003 format(5i7,2i10)
      read(lfntop,1003) nas
      nat=naw+2*(nbw+nhw+ndw+now)
      do 5 i=1,nat
      read(lfntop,1000) cdummy
    5 continue
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
      do 6 i=1,nas
      read(lfntop,1005) j
 1005 format(46x,i5)
      qs(i)=dbl_mb(i_tmp-1+j)
    6 continue
      close(unit=lfntop)
      if(.not.ma_pop_stack(l_tmp))
     + call errquit('Failed to deallocate tmp',0)
c
      endif
c
c     broadcast charges
c
      if(np.gt.1) then
      call ga_brdcst(mag_d01,qs,nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      return
      end
