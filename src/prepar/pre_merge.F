      integer function  pre_merge(iunit,xmerge,filnam,
     + lfnout,lfnpdb,namslv,inum,lrgpdb)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer iunit,lfnout,natm,nseq,nssb,nchain,lfnpdb
      character*255 filnam
      integer numslv,inum
      character*3 namslv
      character*1 altloc,chain
c
      character*80 card
      integer length,mdl,nsgm,isgm,i
      logical lslvnt,lrgpdb
      character*5 cnum
      real*8 x(3),xmerge(3)
c
      logical lrgmrg
c
      natm=0
      nseq=0
      nssb=0
      nchain=0
      mdl=0
      isgm=0
      nsgm=0
      lslvnt=.false.
      numslv=0
      cnum='     '
      lrgmrg=.false.
c
      length=index(filnam,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2006) filnam(1:length)
 2006 format(/,'merge file ',a)
      endif
c
      open(iunit,file=filnam(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filnam(1:length)
 2000 format(' PDB OPEN  ',a)
      endif
c
    1 continue
c
      read(iunit,1000,end=999,err=9999) card
 1000 format(a)
c
      if(card(1:6).eq.'END   ') goto 9
c
c     HEADER
c
      if(card(1:6).eq.'HEADER') then
      if(util_print('where',print_debug)) then
      write(lfnout,2001) card
 2001 format(10x,a)
      endif
      goto 1
      endif
c
c     LRGPDB
c
      if(card(1:6).eq.'LRGPDB') then
      lrgmrg=.true.
      if(.not.lrgpdb) call md_abort('Specify large on write pdb',0)
      goto 1
      endif
c
c     FORMUL
c
      if(card(1:6).eq.'FORMUL'.and.card(19:19).eq.'*') then
      namslv=card(13:15)
      goto 1
      endif
c
c     ATOM
c
      if(card(1:4).eq.'ATOM') then
      if(card(18:20).eq.namslv) lslvnt=.true.
      if((card(17:17).ne.' ').and.
     +   (card(17:17).ne.altloc).and.(altloc.ne.'*')) goto 1
      if(.not.lslvnt) natm=natm+1
      if(lrgmrg) then
      read(card(22:26),1003) i
 1003 format(i5)
      else
      if((card(22:22).ne.' ').and.
     +   (card(22:22).ne.chain).and.(chain.ne.'*')) goto 1
      read(card(23:26),1013) i
 1013 format(i4)
      endif
      isgm=max(isgm,i)
      if(lslvnt) then
      numslv=numslv+1
      else
      nsgm=nsgm+1
      nchain=nchain+1
      if(inum.ge.0) then
      if(card(22:26).ne.cnum) then
      inum=inum+1
      cnum=card(22:26)
      endif
      write(card(22:26),'(i5)') inum
      endif
      read(card(31:54),'(3f8.3)') x
      x(1)=x(1)+xmerge(1)
      x(2)=x(2)+xmerge(2)
      x(3)=x(3)+xmerge(3)
      write(card(31:54),'(3f8.3)') x
      write(lfnpdb,1000) card
      endif
      goto 1
      endif
c
c     HETATM
c
      if(card(1:6).eq.'HETATM') then
      if(card(18:20).eq.namslv) lslvnt=.true.
      if(card(17:17).ne.' '.and.
     + card(17:17).ne.altloc.and.altloc.ne.'*') goto 1
      if(.not.lslvnt) natm=natm+1
      if(lrgmrg) then
      read(card(22:26),1004) i
 1004 format(i5)
      else
      if(card(22:22).ne.' '.and.
     + card(22:22).ne.chain.and.chain.ne.'*') goto 1
      read(card(23:26),1014) i
 1014 format(i4)
      endif
      isgm=max(isgm,i)
      if(lslvnt) then
      numslv=numslv+1
      else
      nsgm=nsgm+1
      nchain=nchain+1
      if(inum.ge.0) then
      if(card(22:26).ne.cnum) then
      inum=inum+1
      cnum=card(22:26)
      endif
      write(card(22:26),'(i5)') inum
      endif
      read(card(31:54),'(3f8.3)') x
      x(1)=x(1)+xmerge(1)
      x(2)=x(2)+xmerge(2)
      x(3)=x(3)+xmerge(3)
      write(card(31:54),'(3f8.3)') x
      write(lfnpdb,1000) card
      endif
      goto 1
      endif
c
      goto 1
c
  999 continue
c
    9 continue
c
      close(iunit)
c
c      inum=isgm
c      pre_merge=isgm
      pre_merge=inum
      return
c
 9999 continue
c
      if(util_print('where',print_none)) then
      write(lfnout,2003) filnam(1:length)
 2003 format('PDB ERROR ',a)
      endif
c
      pre_merge=-9999
      return
      end
