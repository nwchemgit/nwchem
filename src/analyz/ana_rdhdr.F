      subroutine ana_rdhdr(sgmnam)
c
c $Id: ana_rdhdr.F,v 1.2 1999-07-13 14:31:31 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
c
      character*16 sgmnam(msa)
      character*80 card
      character*255 fname
      integer i,na,lq
      character*3 cnum
c
      integer ib,jb,nsb
c
      if(me.eq.0) then
c
      iscof=ifrst
      write(cnum,'(i3.3)') iscof
c
      fname=filsco
      lq=index(filsco,'?')
      if(lq.gt.0) then
      fname=filsco(1:lq-1)//cnum//filsco(lq+1:index(filsco,' ')-1)
      endif
c
      open(unit=lfnsco,file=fname(1:index(fname,' ')-1),
     + status='old',err=9999)
      write(*,'(a,a)') 'Opening sco file ',fname(1:index(fname,' ')-1)
    1 continue
      read(lfnsco,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:17).ne.'num_solute_atoms:') goto 1
      read(card(18:25),'(i8)') na
      if(na.ne.nsa) call errquit('Inconsistent files',0)
    2 continue
      read(lfnsco,1000,err=9999,end=9999) card
      if(card(1:11).ne.'atom_names:') goto 2
      do 3 i=1,nsa
      read(lfnsco,1000,err=9999,end=9999) card
      if(card(16:16).eq.' ') then
      sgmnam(i)=card(5:9)//card(1:4)//' '//card(10:15)
      else
      sgmnam(i)=card(5:9)//card(1:4)//' '//card(11:16)
      endif
    3 continue
      read(lfnsco,1000,err=9999,end=9999) card
      if(card(1:10).eq.'num_bonds:') then
      open(unit=44,file='bond.list',form='formatted',
     + status='unknown')
      rewind(44)
      read(card(12:19),'(i8)') nsb
      write(44,4401) card(1:20)
 4401 format(a)
      read(lfnsco,1000,err=9999,end=9999) card
      do 4 i=1,nsb
      read(lfnsco,1001) ib,jb
 1001 format(2i8)
      write(44,1001) ib,jb
    4 continue
      close(unit=44)
      endif
      endif
c
      return
c
 9999 continue
      call errquit('Failed to open trajectory file',0)
      return
      end
