      logical function pre_modify(lfnout,lfnmod,filmod,latm,catm,qatm,
     + matm,natm)
c
c     $Id: pre_modify.F,v 1.1 1998-06-11 18:51:53 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfnmod,matm,natm
      character*255 filmod
      integer latm(8,matm)
      character*6 catm(4,matm)
      real*8 qatm(6,matm)
c
      character*80 card
      character*12 target
      integer isgm
      character*6 aname
      integer ndx,ndy,mset
      real*8 charge,polar
      character*5 type
c
      integer i,j
c
      open(unit=lfnmod,file=filmod(1:index(filmod,' ')-1),
     + form='formatted',status='old',err=999)
      rewind(unit=lfnmod)
    1 continue
      read(lfnmod,1000,end=99) card
 1000 format(a)
c
c     atom modifications
c     ------------------
c
      if(card(1:4).eq.'atom') then
      target=card(6:17)
      isgm=0
      aname='      '
      ndx=index(target,':')
      ndy=index(target,' ')
      if(ndx.gt.0) then
      read(target(1:ndx-1),*) isgm
      if(ndy-ndx.gt.5)
     + call errquit('pre_modify: atom name problem',9999)
      aname(1:ndy-ndx)=target(ndx+1:ndy)
      else
      if(ndy.gt.5) call errquit('pre_modify: atom name problem',9999)
      aname(1:ndy-1)=target(1:ndy-1)
      endif
      if(aname(1:1).eq.'_') aname(1:1)=' '
c
      read(card(19:19),'(i1)') mset
c
      if(card(21:23).eq.'chg') then
      read(card(25:36),'(f12.6)') charge
      if(util_print('where',print_default)) then
      write(lfnout,2000) isgm,aname,mset,charge
 2000 format('modify atom ',i5,':',a6,' set',i2,' charge ',f12.6)
      endif
      do 2 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.aname) then
      do 3 j=1,3
      if(mset.eq.0.or.mset.eq.j) qatm(2*j-1,i)=charge
    3 continue
      endif
    2 continue
      endif
      if(card(21:23).eq.'pol') then
      read(card(25:36),'(f12.9)') polar
      if(util_print('where',print_default)) then
      write(lfnout,2001) isgm,aname,mset,polar
 2001 format('modify atom ',i5,':',a6,' set',i2,' polar ',f12.9)
      endif
      do 4 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.aname) then
      do 5 j=1,3
      if(mset.eq.0.or.mset.eq.j) qatm(2*j,i)=polar
    5 continue
      endif
    4 continue
      endif
      if(card(21:23).eq.'typ') then
      type=card(25:29)
      if(util_print('where',print_default)) then
      write(lfnout,2002) isgm,aname,mset,type
 2002 format('modify atom ',i5,':',a6,' set',i2,' type ',a5)
      endif
      do 6 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.aname) then
      do 7 j=1,3
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(1:5)=type
    7 continue
      endif
    6 continue
      endif
      if(card(21:23).eq.'dum') then
      if(util_print('where',print_default)) then
      write(lfnout,2003) isgm,aname,mset
 2003 format('modify atom ',i5,':',a6,' set',i2,' dummy ')
      endif
      do 8 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.aname) then
      do 9 j=1,3
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='D'
    9 continue
      endif
    8 continue
      endif
      if(card(21:23).eq.'ego') then
      if(util_print('where',print_default)) then
      write(lfnout,2004) isgm,aname,mset
 2004 format('modify atom ',i5,':',a6,' set',i2,' self ')
      endif
      do 10 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.aname) then
      do 11 j=1,3
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='S'
   11 continue
      endif
   10 continue
      endif
      if(card(21:23).eq.'qtm') then
      if(util_print('where',print_default)) then
      write(lfnout,2005) isgm,aname,mset
 2005 format('modify atom ',i5,':',a6,' set',i2,' quantum ')
      endif
      do 12 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.aname) then
      do 13 j=1,3
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='Q'
   13 continue
      endif
   12 continue
      endif
      endif
c
c     segment modifications
c     ---------------------
c
      if(card(1:7).eq.'segment') then
      read(card(9:13),'(i5)') isgm
      read(card(19:19),'(i1)') mset
      if(card(21:23).eq.'dum') then
      if(util_print('where',print_default)) then
      write(lfnout,2006) isgm,mset
 2006 format('modify segment ',i5,' set',i2,' dummy ')
      endif
      endif
      if(card(21:23).eq.'ego') then
      if(util_print('where',print_default)) then
      write(lfnout,2007) isgm,mset
 2007 format('modify segment ',i5,' set',i2,' self ')
      endif
      endif
      if(card(21:23).eq.'qtm') then
      if(util_print('where',print_default)) then
      write(lfnout,2008) isgm,mset
 2008 format('modify segment ',i5,' set',i2,' quantum ')
      endif
      endif
      endif
c
      goto 1
   99 continue
      close(unit=lfnmod)
  999 continue
      pre_modify=.true.
      return
      end
