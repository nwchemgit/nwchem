      subroutine ana_select(sgmnam)
c
c $Id: ana_select.F,v 1.1 1999-01-05 23:46:39 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
c
      call ana_isel(sgmnam,int_mb(i_isel),dbl_mb(i_wt))
c
      return
      end
      subroutine ana_isel(sgmnam,isel,wt)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa)
c
      integer i,j
      character*255 select,fpref
      integer is,js,it
      character*5 cs,ct
c
c     select
c
      do 1 i=1,nsa
      isel(i)=1
    1 continue
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,fpref))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      select=fpref(1:index(fpref,' ')-1)//'.sel '
      open(unit=69,file=select(1:index(select,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(unit=69)
c
      do 2 i=1,nsa
      isel(i)=0
    2 continue
c
  100 continue
      read(69,6901,end=9998) is,js,cs
 6901 format(2i7,a5)
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 3 i=1,nsa
      read(sgmnam(i),'(5x,a5,i6)') ct,it
      do 4 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
    4 continue
      if((is.eq.0.or.(is.le.it.and.js.ge.it)).and.
     + (cs.eq.'     '.or.cs.eq.ct)) isel(i)=1
    3 continue
      goto 100
c
 9998 continue
      close(unit=69)
c
 9999 continue
      do 5 i=1,nsa
      if(isel(i).eq.0) wt(i)=0.0d0
    5 continue
      write(*,'(100i1)') (isel(i),i=1,nsa)
      return
      end
