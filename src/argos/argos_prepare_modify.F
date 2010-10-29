      logical function argos_prepare_modify(iwhich,lfnout,lfnmod,filmod,
     + nparms,mparms,latm,catm,qatm,matm,natm,lseq,cseq,mseq,nseq,
     + lbnd,rbnd,mbnd,nbnd,nbndt,lang,rang,mang,nang,nangt,
     + ldih,nfdih,kdih,rdih,mdih,ndih,ndiht,limp,kimp,rimp,mimp,nimp,
     + nimpt,lupdat)
c
c     $Id$
c
      implicit none
c
#include "inp.fh"
#include "util.fh"
#include "argos_prepare_common.fh"
c
      logical str_replace
      external str_replace
c
      integer iwhich,lfnout,lfnmod,matm,natm,mseq,nseq,nparms,mparms
      character*255 filmod
      integer latm(11,matm)
      character*6 catm(mparms,matm)
      integer lseq(6,mseq)
      character*10 cseq(mseq)
      real*8 qatm(nparms,2,matm)
      integer mbnd,nbnd,nbndt
      integer lbnd(4,mbnd)
      real*8 rbnd(nparms,2,mbnd)
      integer mang,nang,nangt
      integer lang(5,mang)
      real*8 rang(nparms,4,mang)
      integer mdih,ndih,ndiht
      integer ldih(6,mdih),nfdih(nparms,mdih),kdih(6,nparms,mdih)
      real*8 rdih(6,nparms,2,mdih)
      integer mimp,nimp,nimpt
      integer limp(6,mimp),kimp(nparms,mimp)
      real*8 rimp(nparms,2,mimp)
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
      logical foundit,lupdat
c
      integer i,j,k,l,m,n
c
      lupdat=.false.
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
     + call md_abort('argos_prepare_modify: atom name problem',9999)
      aname(1:ndy-ndx)=target(ndx+1:ndy)
      else
      if(ndy.gt.5) 
     + call md_abort('argos_prepare_modify: atom name problem',9999)
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
 2000 format(' modify atom ',i5,':',a6,' set',i2,' charge ',f12.6)
      endif
      do 2 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      do 3 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) qatm(j,1,i)=charge
    3 continue
      endif
    2 continue
      endif
      if(card(21:23).eq.'pol') then
      read(card(25:36),'(f12.9)') polar
      if(util_print('where',print_default)) then
      write(lfnout,2001) isgm,aname,mset,polar
 2001 format(' modify atom ',i5,':',a6,' set',i2,' polar ',f12.9)
      endif
      do 4 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      do 5 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) qatm(j,2,i)=polar
    5 continue
      endif
    4 continue
      endif
      if(card(21:23).eq.'typ') then
      type=card(25:29)
      if(util_print('where',print_default)) then
      write(lfnout,2002) isgm,aname,mset,type
 2002 format(' modify atom ',i5,':',a6,' set',i2,' type ',a5)
      endif
      do 6 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      do 7 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(1:5)=type
    7 continue
      endif
    6 continue
      endif
      if(card(21:23).eq.'dum') then
      if(util_print('where',print_default)) then
      write(lfnout,2003) isgm,aname,mset
 2003 format(' modify atom ',i5,':',a6,' set',i2,' dummy ')
      endif
      do 8 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      do 9 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='D'
    9 continue
      endif
    8 continue
      endif
      if(card(21:23).eq.'ego') then
      if(util_print('where',print_default)) then
      write(lfnout,2004) isgm,aname,mset
 2004 format(' modify atom ',i5,':',a6,' set',i2,' self ')
      endif
      do 10 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      do 11 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='S'
   11 continue
      endif
   10 continue
      endif
      if(card(21:23).eq.'qtm') then
      if(util_print('where',print_default)) then
      write(lfnout,2005) isgm,aname,mset
 2005 format(' modify atom ',i5,':',a6,' set',i2,' quantum ')
      endif
      foundit=.false.
      do 12 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      foundit=.true.
      do 13 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='Q'
   13 continue
      endif
   12 continue
      if(.not.foundit) call md_abort('Atom not found for quantum',0)
      endif
      if(card(21:23).eq.'QTM') then
      if(util_print('where',print_default)) then
      write(lfnout,2015) isgm,aname,mset
 2015 format(' modify atom ',i5,':',a6,' set',i2,' high level quantum ')
      endif
      foundit=.false.
      do 112 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.aname) then
      foundit=.true.
      do 113 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='H'
  113 continue
      endif
  112 continue
      if(.not.foundit) call md_abort('Atom not found for quantum',0)
      endif
      endif
c
c     segment modifications
c     ---------------------
c
      if(card(1:7).eq.'segment') then
c      print*,'CARD=',card
      read(card(9:13),'(i5)') isgm
c      print*,'ISGM=',isgm
      read(card(19:19),'(i1)') mset
c      print*,'MSET=',mset
      if(card(21:23).eq.'dum') then
      if(util_print('where',print_default)) then
      write(lfnout,2006) isgm,mset
 2006 format(' modify segment ',i5,' set',i2,' dummy ')
      endif
      do 14 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm) then
      do 15 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='D'
   15 continue
      endif
   14 continue
      endif
      if(card(21:23).eq.'ego') then
      if(util_print('where',print_default)) then
      write(lfnout,2007) isgm,mset
 2007 format(' modify segment ',i5,' set',i2,' self ')
      endif
      do 16 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm) then
      do 17 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='S'
   17 continue
      endif
   16 continue
      endif
      if(card(21:23).eq.'noq') then
      if(util_print('where',print_default)) then
      write(lfnout,2108) isgm,mset
 2108 format(' modify segment ',i5,' set',i2,' nocharge ')
      endif
      do 218 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm) then
      do 219 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) qatm(j,1,i)=0.0d0
  219 continue
      endif
  218 continue
      endif
      if(card(21:23).eq.'qtm') then
      if(util_print('where',print_default)) then
      write(lfnout,2008) isgm,mset
 2008 format(' modify segment ',i5,' set',i2,' quantum ')
      endif
      do 18 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm) then
      do 19 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='Q'
   19 continue
      endif
   18 continue
      endif
      if(card(21:23).eq.'QTM') then
      if(util_print('where',print_default)) then
      write(lfnout,2018) isgm,mset
 2018 format(' modify segment ',i5,' set',i2,' high level quantum ')
      endif
      do 118 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm) then
      do 119 j=1,nparms
      if(mset.eq.0.or.mset.eq.j) catm(j+1,i)(6:6)='H'
  119 continue
      endif
  118 continue
      endif
      endif
c
      if(card(1:6).eq.'update') then
      lupdat=.true.
      do 50 j=1,nparms
      n=0
      m=0
      do 51 i=1,nbndt
      if((catm(j+1,lbnd(1,i))(6:6).ne.'Q'.and.
     + catm(j+1,lbnd(1,i))(6:6).ne.'H').or.
     + (catm(j+1,lbnd(2,i))(6:6).ne.'Q'.and.
     + catm(j+1,lbnd(2,i))(6:6).ne.'H')) then
      n=n+1
      do 52 k=1,4
      lbnd(k,n)=lbnd(k,i)
   52 continue
      do 552 k=1,nparms
      rbnd(k,1,n)=rbnd(k,1,i)
      rbnd(k,2,n)=rbnd(k,2,i)
  552 continue
      if(i.le.nbnd) m=n
      endif
   51 continue
      nbnd=m
      nbndt=n
      n=0
      m=0
      do 53 i=1,nangt
      if((catm(j+1,lang(1,i))(6:6).ne.'Q'.and.
     + catm(j+1,lang(1,i))(6:6).ne.'H').or.
     + (catm(j+1,lang(2,i))(6:6).ne.'Q'.and.
     + catm(j+1,lang(2,i))(6:6).ne.'H').or.
     + (catm(j+1,lang(3,i))(6:6).ne.'Q'.and.
     + catm(j+1,lang(3,i))(6:6).ne.'H')) then
      n=n+1
      do 54 k=1,5
      lang(k,n)=lang(k,i)
   54 continue
      do 554 k=1,nparms
      rang(k,1,n)=rang(k,1,i)
      rang(k,2,n)=rang(k,2,i)
      if(ffield(1:6).eq.'charmm') then
      rang(k,3,n)=rang(k,3,i)
      rang(k,4,n)=rang(k,4,i)
      endif
  554 continue
      if(i.le.nang) m=n
      endif
   53 continue
      nang=m
      nangt=n
      n=0
      m=0
      do 55 i=1,ndiht
      if((catm(j+1,ldih(1,i))(6:6).ne.'Q'.and.
     + catm(j+1,ldih(1,i))(6:6).ne.'H').or.
     + (catm(j+1,ldih(2,i))(6:6).ne.'Q'.and.
     + catm(j+1,ldih(2,i))(6:6).ne.'H').or.
     + (catm(j+1,ldih(3,i))(6:6).ne.'Q'.and.
     + catm(j+1,ldih(3,i))(6:6).ne.'H').or.
     + (catm(j+1,ldih(4,i))(6:6).ne.'Q'.and.
     + catm(j+1,ldih(4,i))(6:6).ne.'H')) then
      n=n+1
      do 56 k=1,6
      ldih(k,n)=ldih(k,i)
   56 continue
      do 57 l=1,nparms
      nfdih(l,n)=nfdih(l,i)
      do 556 k=1,6
      rdih(k,l,1,n)=rdih(k,l,1,i)
      rdih(k,l,2,n)=rdih(k,l,2,i)
  556 continue
   57 continue
      if(i.le.ndih) m=n
      endif
   55 continue
      ndih=m
      ndiht=n
      n=0
      m=0
      do 58 i=1,nimpt
      if((catm(j+1,limp(1,i))(6:6).ne.'Q'.and.
     + catm(j+1,limp(1,i))(6:6).ne.'H').or.
     + (catm(j+1,limp(2,i))(6:6).ne.'Q'.and.
     + catm(j+1,limp(2,i))(6:6).ne.'H').or.
     + (catm(j+1,limp(3,i))(6:6).ne.'Q'.and.
     + catm(j+1,limp(3,i))(6:6).ne.'H').or.
     + (catm(j+1,limp(4,i))(6:6).ne.'Q'.and.
     + catm(j+1,limp(4,i))(6:6).ne.'H')) then
      n=n+1
      do 59 k=1,6
      limp(k,n)=limp(k,i)
   59 continue
      do 559 k=1,nparms
      kimp(k,n)=kimp(k,i)
      rimp(k,1,n)=rimp(k,1,i)
      rimp(k,2,n)=rimp(k,2,i)
  559 continue
      if(i.le.nimp) m=n
      endif
   58 continue
      nimp=m
      nimpt=n
   50 continue
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
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call md_abort('str_replace',9999)
      iatom=0
      jatom=0
      do 20 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(lseq(1,latm(5,i)).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
   20 continue
      if(iatom.gt.0.and.jatom.gt.0) then
      ibnd=0
      do 21 i=1,nbndt
      if((lbnd(1,i).eq.iatom.and.lbnd(2,i).eq.jatom).or.
     + (lbnd(1,i).eq.jatom.and.lbnd(2,i).eq.iatom)) then
      ibnd=i
      goto 22
      endif
   21 continue
   22 continue
      if(ibnd.eq.0) call md_abort('Bond could not be found',9999)
      do 23 i=1,3
      if(mset.eq.0.or.mset.eq.i) then
      rbnd(i,1,ibnd)=value
      rbnd(i,2,ibnd)=forcon
      endif
   23 continue
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
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) ksgm
      atomk=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomk,'_',' '))
     + call md_abort('str_replace',9999)
      iatom=0
      jatom=0
      katom=0
      do 24 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(lseq(1,latm(5,i)).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
      if(lseq(1,latm(5,i)).eq.ksgm.and.catm(1,i).eq.atomk(1:6)) katom=i
   24 continue
      if(iatom.gt.0.and.jatom.gt.0.and.katom.gt.0) then
      iang=0
      do 25 i=1,nangt
      if(lang(2,i).eq.jatom.and.
     + ((lang(1,i).eq.iatom.and.lang(3,i).eq.katom).or.
     + (lang(1,i).eq.katom.and.lang(3,i).eq.iatom))) then
      iang=i
      goto 26
      endif
   25 continue
   26 continue
      if(iang.eq.0) call md_abort('Angle could not be found',9999)
      do 27 i=1,3
      if(mset.eq.0.or.mset.eq.i) then
      rang(i,1,iang)=value
      rang(i,2,iang)=forcon
      endif
   27 continue
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
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) ksgm
      atomk=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomk,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) lsgm
      atoml=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atoml,'_',' '))
     + call md_abort('str_replace',9999)
      iatom=0
      jatom=0
      katom=0
      latom=0
      do 28 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(lseq(1,latm(5,i)).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
      if(lseq(1,latm(5,i)).eq.ksgm.and.catm(1,i).eq.atomk(1:6)) katom=i
      if(lseq(1,latm(5,i)).eq.lsgm.and.catm(1,i).eq.atoml(1:6)) latom=i
   28 continue
      if(iatom.gt.0.and.jatom.gt.0.and.katom.gt.0.and.latom.gt.0) then
      idih=0
      do 29 i=1,ndiht
      if((ldih(1,i).eq.iatom.and.ldih(2,i).eq.jatom.and.
     + ldih(3,i).eq.katom.and.ldih(4,i).eq.latom).or.
     + (ldih(1,i).eq.latom.and.ldih(2,i).eq.katom.and.
     + ldih(3,i).eq.jatom.and.ldih(4,i).eq.iatom)) then
      idih=i
      goto 30
      endif
   29 continue
   30 continue
      if(idih.eq.0) call md_abort('Torsion could not be found',9999)
      m=0
      do 531 j=1,nparms
      m=max(m,nfdih(j,i))
  531 continue
      do 31 j=1,m
      do 32 i=1,nparms
      if((mset.eq.0.or.mset.eq.i).and.
     + (multip.eq.0.or.multip.eq.iabs(kdih(j,i,idih)))) then
      rdih(j,i,1,idih)=value
      rdih(j,i,2,idih)=forcon
      endif
   32 continue
      if(multip.eq.0) goto 33
   31 continue
   33 continue
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
      call md_abort('Unable to find atom',9999)
      endif
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) ksgm
      atomk=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomk,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) lsgm
      atoml=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atoml,'_',' '))
     + call md_abort('str_replace',9999)
      iatom=0
      jatom=0
      katom=0
      latom=0
      do 34 i=1,natm
      if(lseq(1,latm(5,i)).eq.isgm.and.catm(1,i).eq.atomi(1:6)) iatom=i
      if(lseq(1,latm(5,i)).eq.jsgm.and.catm(1,i).eq.atomj(1:6)) jatom=i
      if(lseq(1,latm(5,i)).eq.ksgm.and.catm(1,i).eq.atomk(1:6)) katom=i
      if(lseq(1,latm(5,i)).eq.lsgm.and.catm(1,i).eq.atoml(1:6)) latom=i
   34 continue
      if(iatom.gt.0.and.jatom.gt.0.and.katom.gt.0.and.latom.gt.0) then
      iimp=0
      do 35 i=1,nimpt
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
      goto 36
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
      goto 36
      endif
      endif
   35 continue
   36 continue
      if(iimp.eq.0) call
     + md_abort('Out-of-plane could not be found',9999)
      do 37 i=1,3
      if(mset.eq.0.or.mset.eq.i) then
      rimp(i,1,iimp)=value
      rimp(i,2,iimp)=forcon
      endif
   37 continue
      endif
      endif
c
c
      endif
c
      goto 1
   99 continue
      close(unit=lfnmod)
  999 continue
      argos_prepare_modify=.true.
      return
      end
