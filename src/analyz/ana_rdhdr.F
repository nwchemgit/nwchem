      subroutine ana_rdhdr(sgmnam)
c
c $Id: ana_rdhdr.F,v 1.5 1999-12-10 23:54:47 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
c
      character*16 wnam,sgmnam(msa)
      character*80 card
      character*255 fname
      integer i,na,lq
      character*3 cnum
      character*1 cdum
c
      integer ib,jb,nsb
c
      if(me.eq.0) then
c
      iscof=ifrst
      write(cnum,'(i3.3)') iscof
c
      fname=filtrj
      lq=index(filtrj,'?')
      if(lq.gt.0) then
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
      endif
c
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9999)
    1 continue
      read(lfntrj,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 1
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(lfntrj,1003) nwa,nsa,nsb
      write(44,1003) nwa,nsa,nsb
 1003 format(3i10)
      if(nsa.gt.msa) call errquit('Error in trj',0)
      if(nwa.gt.0) then
      do 2 i=1,nwa
      read(lfntrj,1004) wnam
      write(44,1004) wnam
 1004 format(a16)
    2 continue
      endif
      if(nsa.gt.0) then
      do 3 i=1,nsa
      read(lfntrj,1005) sgmnam(i)
      write(44,1005) sgmnam(i)
 1005 format(a16)
    3 continue
      endif
      if(nsb.gt.0) then
      do 4 i=1,nsb
      read(lfntrj,1001) ib,jb
      write(44,1001) ib,jb
 1001 format(2i8)
    4 continue
      endif
      close(unit=44)
c
      endif
c
      return
c
 9999 continue
      call errquit('Failed to open trajectory file',0)
      return
      end
