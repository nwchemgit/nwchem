      logical function pre_modify(iwhich,lfnout,lfnmod,filmod,ffield,
     + latm,catm,qatm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,kdih,rdih,mdih,ndih,limp,rimp,mimp,nimp)
c
c     $Id: pre_modify.F,v 1.3 1998-09-17 16:38:44 d3j191 Exp $
c
      implicit none
c
#include "inp.fh"
#include "util.fh"
c
      logical str_replace
      external str_replace
c
      integer iwhich,lfnout,lfnmod,matm,natm
      character*80 ffield
      character*255 filmod
      integer latm(8,matm)
      character*6 catm(4,matm)
      real*8 qatm(6,matm)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(6,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(6,mang)
      integer mdih,ndih
      integer ldih(9,mdih),kdih(6,3,mdih)
      real*8 rdih(6,6,mdih)
      integer mimp,nimp
      integer limp(9,mimp)
      real*8 rimp(6,mimp)
c
      character*80 card
      character*12 target,atomi,atomj,atomk,atoml
      integer isgm,jsgm,ksgm,lsgm,ifr,ito
      integer iatom,jatom,katom,latom
      integer ibnd,iang,idih,iimp
      character*6 aname
      integer ndx,ndy,mset,multip
      real*8 charge,polar,value,forcon
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
      if(iwhich.eq.1) then
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
      if(aname(2:2).eq.'_') aname(2:2)=' '
      if(aname(3:3).eq.'_') aname(3:3)=' '
      if(aname(4:4).eq.'_') aname(4:4)=' '
      if(aname(5:5).eq.'_') aname(5:5)=' '
      if(aname(6:6).eq.'_') aname(6:6)=' '
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
      endif
c
      if(iwhich.eq.2) then
c
c     bond modifications
c     ------------------
c
      if(card(1:4).eq.'bond') then
      read(card(6:6),'(i1)') mset
      read(card(8:19),'(f12.6)') value
      read(card(20:31),'(e12.5)') forcon
      card=card(33:80)
      ifr=0
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call errquit('str_replace',9999)
      iatom=0
      jatom=0
      do 14 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(latm(5,i).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
   14 continue
      if(iatom.gt.0.and.jatom.gt.0) then
      ibnd=0
      do 15 i=1,nbnd
      if((lbnd(1,i).eq.iatom.and.lbnd(2,i).eq.jatom).or.
     + (lbnd(1,i).eq.jatom.and.lbnd(2,i).eq.iatom)) then
      ibnd=i
      goto 16
      endif
   15 continue
   16 continue
      if(ibnd.eq.0) call errquit('Bond could not be found',9999)
      do 17 i=1,3
      if(mset.eq.0.or.mset.eq.i) then
      rbnd(2*i-1,ibnd)=value
      rbnd(2*i,ibnd)=forcon
      endif
   17 continue
      endif
      endif
c
c     angle modifications
c     ------------------
c
      if(card(1:5).eq.'angle') then
      read(card(7:7),'(i1)') mset
      read(card(9:20),'(f12.6)') value
      read(card(21:32),'(e12.5)') forcon
      card=card(34:80)
      ifr=0
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) ksgm
      atomk=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomk,'_',' '))
     + call errquit('str_replace',9999)
      iatom=0
      jatom=0
      katom=0
      do 18 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(latm(5,i).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
      if(latm(5,i).eq.ksgm.and.catm(1,i).eq.atomk(1:6)) katom=i
   18 continue
      if(iatom.gt.0.and.jatom.gt.0.and.katom.gt.0) then
      iang=0
      do 19 i=1,nang
      if(lang(2,i).eq.jatom.and.
     + ((lang(1,i).eq.iatom.and.lang(3,i).eq.katom).or.
     + (lang(1,i).eq.katom.and.lang(3,i).eq.iatom))) then
      iang=i
      goto 20
      endif
   19 continue
   20 continue
      if(iang.eq.0) call errquit('Angle could not be found',9999)
      do 21 i=1,3
      if(mset.eq.0.or.mset.eq.i) then
      rang(2*i-1,iang)=value
      rang(2*i,iang)=forcon
      endif
   21 continue
      endif
      endif
c
c     torsion modifications
c     ---------------------
c
      if(card(1:7).eq.'torsion') then
      read(card(9:9),'(i1)') mset
      read(card(11:11),'(i1)') multip
      read(card(13:24),'(f12.6)') value
      read(card(25:36),'(e12.5)') forcon
      card=card(38:80)
      ifr=0
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) ksgm
      atomk=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomk,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) lsgm
      atoml=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atoml,'_',' '))
     + call errquit('str_replace',9999)
      iatom=0
      jatom=0
      katom=0
      latom=0
      do 22 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(latm(5,i).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
      if(latm(5,i).eq.ksgm.and.catm(1,i).eq.atomk(1:6)) katom=i
      if(latm(5,i).eq.lsgm.and.catm(1,i).eq.atoml(1:6)) latom=i
   22 continue
      if(iatom.gt.0.and.jatom.gt.0.and.katom.gt.0.and.latom.gt.0) then
      idih=0
      do 23 i=1,ndih
      if((ldih(1,i).eq.iatom.and.ldih(2,i).eq.jatom.and.
     + ldih(3,i).eq.katom.and.ldih(4,i).eq.latom).or.
     + (ldih(1,i).eq.latom.and.ldih(2,i).eq.katom.and.
     + ldih(3,i).eq.jatom.and.ldih(4,i).eq.iatom)) then
      idih=i
      goto 24
      endif
   23 continue
   24 continue
      if(idih.eq.0) call errquit('Torsion could not be found',9999)
      do 25 j=1,max(ldih(7,i),ldih(8,i),ldih(9,i))
      do 26 i=1,3
      if((mset.eq.0.or.mset.eq.i).and.
     + (multip.eq.0.or.multip.eq.iabs(kdih(j,i,idih)))) then
      rdih(j,2*i-1,idih)=value
      rdih(j,2*i,idih)=forcon
      endif
   26 continue
      if(multip.eq.0) goto 27
   25 continue
   27 continue
      endif
      endif
c
c     out-of-plane modifications
c     --------------------------
c
      if(card(1:5).eq.'plane') then
      read(card(7:7),'(i1)') mset
      read(card(9:20),'(f12.6)') value
      read(card(21:32),'(e12.5)') forcon
      card=card(34:80)
      ifr=0
      if(.not.inp_strtok(card,' ',ifr,ito)) then
      call errquit('Unable to find atom',9999)
      endif
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) ksgm
      atomk=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomk,'_',' '))
     + call errquit('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) lsgm
      atoml=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atoml,'_',' '))
     + call errquit('str_replace',9999)
      iatom=0
      jatom=0
      katom=0
      latom=0
      do 28 i=1,natm
      if(latm(5,i).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(latm(5,i).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
      if(latm(5,i).eq.ksgm.and.catm(1,i).eq.atomk(1:6)) katom=i
      if(latm(5,i).eq.lsgm.and.catm(1,i).eq.atoml(1:6)) latom=i
   28 continue
      if(iatom.gt.0.and.jatom.gt.0.and.katom.gt.0.and.latom.gt.0) then
      iimp=0
      do 29 i=1,nimp
      if(ffield(1:5).eq.'amber') then
      if(limp(1,i).eq.katom.and.
     + ((limp(2,i).eq.iatom.and.limp(3,i).eq.jatom.and.
     + limp(4,i).eq.latom).or.
     + (limp(2,i).eq.iatom.and.limp(3,i).eq.latom.and.
     + limp(4,i).eq.jatom).or.
     + (limp(2,i).eq.jatom.and.limp(3,i).eq.iatom.and.
     + limp(4,i).eq.latom).or.
     + (limp(2,i).eq.jatom.and.limp(3,i).eq.latom.and.
     + limp(4,i).eq.iatom).or.
     + (limp(2,i).eq.latom.and.limp(3,i).eq.iatom.and.
     + limp(4,i).eq.jatom).or.
     + (limp(2,i).eq.latom.and.limp(3,i).eq.jatom.and.
     + limp(4,i).eq.iatom))) then
      iimp=i
      goto 30
      endif
      else
      if(limp(3,i).eq.katom.and.
     + ((limp(1,i).eq.iatom.and.limp(2,i).eq.jatom.and.
     + limp(4,i).eq.latom).or.
     + (limp(1,i).eq.jatom.and.limp(2,i).eq.iatom.and.
     + limp(4,i).eq.latom).or.
     + (limp(1,i).eq.iatom.and.limp(2,i).eq.latom.and.
     + limp(4,i).eq.latom).or.
     + (limp(1,i).eq.jatom.and.limp(2,i).eq.latom.and.
     + limp(4,i).eq.iatom).or.
     + (limp(1,i).eq.latom.and.limp(2,i).eq.iatom.and.
     + limp(4,i).eq.jatom).or.
     + (limp(1,i).eq.latom.and.limp(2,i).eq.jatom.and.
     + limp(4,i).eq.iatom))) then
      iimp=i
      goto 30
      endif
      endif
   29 continue
   30 continue
      if(iimp.eq.0) call errquit('Out-of-plane could not be found',9999)
      do 31 i=1,3
      if(mset.eq.0.or.mset.eq.i) then
      rimp(2*i-1,iimp)=value
      rimp(2*i,iimp)=forcon
      endif
   31 continue
      endif
      endif
c
      endif
c
      goto 1
   99 continue
      close(unit=lfnmod)
  999 continue
      pre_modify=.true.
      return
      end
      logical function str_replace(string,charf,chart)
      implicit none
#include "inp.fh"
      character*(*) string
      character*1 charf,chart
      integer i
      do 1 i=1,inp_strlen(string)
      if(string(i:i).eq.charf) string(i:i)=chart
    1 continue
      str_replace=.true.
      return
      end
