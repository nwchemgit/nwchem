      subroutine dia_rdtop(sgmnam,qs,qw,iram,imol,iatmn,ibnd)
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "dia_common.fh"
c
      character*16 sgmnam(msa)
      real*8 qs(msa),qw(mwa)
      integer iram(msgm,7),imol(msa),iatmn(msa),ibnd(msb,2)
c
      character*1 cdummy
      integer i,j,k,nat,naq,nparm,nseq,i_tmp,l_tmp,i_itmp,l_itmp
      integer naw,nbw,nhw,ndw,now,ntw,nnw
      integer nas,nbs,num
      character*5 sname,aname
      integer i_atm,l_atm
c
      ltop=.false.
c
      do 1 i=1,nsgm
      iram(i,1)=0
      iram(i,2)=0
      iram(i,3)=0
      iram(i,4)=0
      iram(i,5)=0
      iram(i,6)=0
      iram(i,7)=0
    1 continue
c
      if(me.eq.0) then
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=999)
c
c      print*,'TOPOLOGY FILE ',filtop(1:index(filtop,' ')-1)
      read(lfntop,1000,end=999,err=999) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
 1000 format(a1)
      read(lfntop,1001) nparm
      read(lfntop,1001) nat
      read(lfntop,1001) naq
      read(lfntop,1001) nseq
      read(lfntop,1000) cdummy
 1001 format(i5)
      if(.not.ma_push_get(mt_int,nat,'atm',l_atm,i_atm))
     + call md_abort('Could not allocate atm',0)
      do 2 i=1,nat
      read(lfntop,1222) int_mb(i_atm-1+i)
 1222 format(5x,i5)
      if(nparm.gt.1) then
      do 222 j=1,nparm-1
      read(lfntop,1000) cdummy
  222 continue
      endif
    2 continue
      do 3 i=1,nat
      do 4 j=i,nat
      do 5 k=1,nparm
      read(lfntop,1000) cdummy
    5 continue
    4 continue
    3 continue
      if(.not.ma_push_get(mt_dbl,naq,'tmp',l_tmp,i_tmp))
     + call md_abort('Failed to allocate tmp',0)
      do 6 i=1,naq
      read(lfntop,1002) dbl_mb(i_tmp-1+i)
 1002 format(5x,f12.6)
      do 7 j=1,nparm-1
      read(lfntop,1000) cdummy
    7 continue
    6 continue
      do 8 i=1,nseq
      read(lfntop,1000) cdummy
    8 continue
      read(lfntop,1003) naw,nbw,nhw,ndw,now,ntw,nnw
 1003 format(5i7,2i10)
      read(lfntop,1003) nas,nbs
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      write(44,2003) naw,nas,nbs
 2003 format(3i10)
      do 9 i=1,naw
      read(lfntop,1005) sname,aname,num,num,j
      write(44,2004) sname,aname
 2004 format(2a5,i6)
      qw(i)=dbl_mb(i_tmp-1+j)
    9 continue
      do 10 i=1,nbw
      read(lfntop,'(2i7)') ibndw(i,1),ibndw(i,2)
      do 11 j=1,nparm
      read(lfntop,1000) cdummy
   11 continue
   10 continue
      nat=4*(nhw+ndw+now)
      do 12 i=1,nat
      read(lfntop,1000) cdummy
   12 continue
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
      do 13 i=1,nparm
      read(lfntop,1000) cdummy
   13 continue
      if(.not.ma_push_get(mt_int,nas,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',0)
      do 14 i=1,nas
      read(lfntop,1005) sname,aname,imol(i),num,k,j
 1005 format(a5,5x,a5,4x,2i7,14x,2i5)
c 1005 format(a5,5x,a5,6x,2i5,15x,i5)
c      write(*,1005) sname,aname,imol(i),num,j
      write(44,2004) sname,aname,num
      write(sgmnam(i),2004) sname,aname,num
      qs(i)=dbl_mb(i_tmp-1+j)
      iatmn(i)=int_mb(i_atm-1+k)
      int_mb(i_itmp-1+i)=num
      if(aname.eq.' N   ') iram(num,2)=i
      if(aname.eq.' CA  ') iram(num,3)=i
      if(aname.eq.' C   ') iram(num,4)=i
      if(aname.eq.' H   ') iram(num,6)=i
      if(aname.eq.' O   ') iram(num,7)=i
   14 continue
c      write(*,'(10i5)') (iatmn(i),i=1,nas)
      do 15 i=1,nbs
      read(lfntop,1006) j,num
      ibnd(i,1)=j
      ibnd(i,2)=num
      do 16 k=1,nparm
      read(lfntop,1000) cdummy
   16 continue
 1006 format(2i7)
      write(44,2005) j,num
 2005 format(2i8)
      if(int_mb(i_itmp-1+j).ne.int_mb(i_itmp-1+num)) then
      if(sgmnam(j)(6:10).eq.' C   '.and.
     + sgmnam(num)(6:10).eq.' N   ') then
      iram(int_mb(i_itmp-1+j),5)=num
      iram(int_mb(i_itmp-1+num),1)=j
      endif
      if(sgmnam(num)(6:10).eq.' C   '.and.
     + sgmnam(j)(6:10).eq.' N   ') then
      iram(int_mb(i_itmp-1+j),1)=num
      iram(int_mb(i_itmp-1+num),5)=j
      endif
      endif
   15 continue
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to deallocate itmp',0)
      close(unit=44)
      close(unit=lfntop)
      if(.not.ma_pop_stack(l_tmp))
     + call md_abort('Failed to deallocate tmp',0)
      if(.not.ma_pop_stack(l_atm))
     + call md_abort('Could not deallocate atm',0)
c
      ltop=.true.
      endif
  999 continue      
c
c     broadcast charges
c
      if(np.gt.1) then
      call ga_brdcst(mag_d01,ltop,ma_sizeof(mt_log,1,mt_byte),0)
      if(ltop) then
      call ga_brdcst(mag_d02,qs,nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d03,iatmn,nsa*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d04,iram,7*msgm*ma_sizeof(mt_int,1,mt_byte),0)
      endif
      endif
c
c      print*,'DIA_RDTOP DONE'
      return
      end
      subroutine dia_siztop()
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "dia_common.fh"
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
