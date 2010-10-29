      logical function argos_prep_rdsgm(lfnout,nparms,mparms,
     + imol,isgm,igrp,
     + ipgrp,lfnsgm,filsgm,lfnmat,filmat,
     + latm,catm,qatm,matm,natm,lbnd,rbnd,mbnd,nbnd,
     + lang,rang,mang,nang,ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,
     + wcorr,npar,ipardef,itopol)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "argos_prep_common.fh"
c
      external loc
      integer loc
c
      integer lfnout,lfnsgm,imol,isgm,igrp,ipgrp,lfnmat,nparms,mparms
      character*255 filsgm,filmat
      integer matm,natm,itopol
      integer latm(11,matm)
      character*6 catm(mparms,matm)
      real*8 qatm(nparms,2,matm)
      integer mbnd,nbnd
      integer lbnd(4,mbnd)
      real*8 rbnd(nparms,2,mbnd)
      integer mang,nang
      integer lang(5,mang)
      real*8 rang(nparms,4,mang)
      integer mdih,ndih
      integer ldih(6,mdih),nfdih(nparms,mdih),kdih(6,nparms,mdih)
      real*8 rdih(6,nparms,2,mdih)
      integer mimp,nimp
      integer limp(6,mimp),kimp(nparms,mimp)
      real*8 rimp(nparms,2,mimp)
      real*8 wcorr(10)
      real*8 version
c
      character*80 card
      integer i,j,k,l,length,na,jmol,nzmat,npar,ipardef,idhop
      integer nsatm,nsbnd,nsang,nsdih,nsimp,ld(6),md(10),izm(4)
      real*8 rd(10,2),zm(3)
c
      integer itemp
      real*8 dtemp
      character*6 ctemp
c
      jmol=imol
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,1110) filsgm(1:length)
 1110 format('READING SEGMENT FILE ',a)
      endif
c
    1 continue
      read(lfnsgm,1000) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
      read(card,1001,err=9999) version
 1001 format(f12.6)
      read(lfnsgm,1002,err=9999)
     + nsatm,nsbnd,nsang,nsdih,nsimp,nzmat,npar,ipardef
 1002 format(8i5)
      if(ipardef.eq.0) ipardef=1
      if(npar.gt.nparms) call md_abort('Error in number par sets',npar)
      do 2 i=1,npar
      read(lfnsgm,1003) wcorr(i)
 1003 format(f12.6)
    2 continue
      do 102 i=npar+1,nparms
      wcorr(i)=wcorr(npar)
  102 continue
c
      if(natm+nsatm.gt.matm) call md_abort('increase matm',9999)
      if(nbnd+nsbnd.gt.mbnd) call md_abort('increase mbnd',9999)
      if(nang+nsang.gt.mang) call md_abort('increase mang',9999)
      if(ndih+nsdih.gt.mdih) call md_abort('increase mdih',9998)
      if(nimp+nsimp.gt.mimp) call md_abort('increase mimp',9999)
c
c     read the atom list
c     ------------------
c
      na=natm
      do 3 i=1,nsatm
      natm=natm+1
      read(lfnsgm,1004) catm(1,natm),latm(3,natm),latm(4,natm),
     + latm(10,natm),latm(1,natm),latm(2,natm)
 1004 format(5x,a6,5i5)
c      write(*,1004) catm(1,natm),latm(3,natm),latm(4,natm),
c     + latm(10,natm),latm(1,natm),latm(2,natm)
      read(lfnsgm,1005) (catm(j+1,natm),qatm(j,1,natm),qatm(j,2,natm),
     + j=1,npar)
 1005 format(5x,a6,2f12.6)
c
      latm(5,natm)=isgm
      latm(6,natm)=jmol
      jmol=iabs(jmol)
      latm(1,natm)=latm(1,natm)+igrp
      latm(2,natm)=latm(2,natm)+ipgrp
c
      if(isgm.eq.0) then
      do 4 j=2,npar+1
      catm(j,natm)(6:6)='w'
    4 continue
      endif
c
      do 103 j=npar+1,nparms
      catm(j+1,natm)=catm(npar+1,natm)
      qatm(j,1,natm)=qatm(npar,1,natm)
      qatm(j,2,natm)=qatm(npar,2,natm)
  103 continue
c
      idhop=0
      do 1103 j=1,npar
      if(catm(j+1,natm)(6:6).ne.'D') idhop=idhop+2**(j-1)
 1103 continue
      latm(11,natm)=idhop
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1004) catm(1,natm),(latm(j,natm),j=1,4)
      write(lfnout,1005) (catm(j+1,natm),qatm(j,1,natm),qatm(j,2,natm),
     + j=1,nparms)
      endif
c
      if(itopol.eq.0.and.ipardef.gt.1) then
      ctemp=catm(2,natm)
      catm(2,natm)=catm(1+ipardef,natm)
      catm(1+ipardef,natm)=ctemp
      dtemp=qatm(1,1,natm)
      qatm(1,1,natm)=qatm(ipardef,1,natm)
      qatm(ipardef,1,natm)=dtemp
      dtemp=qatm(1,2,natm)
      qatm(1,2,natm)=qatm(ipardef,2,natm)
      qatm(ipardef,2,natm)=dtemp
      endif
    3 continue
c
      igrp=latm(1,natm)
      ipgrp=latm(2,natm)
c
c     read the bond list
c     ------------------
c
      do 5 i=1,nsbnd
      nbnd=nbnd+1
      read(lfnsgm,1006) (lbnd(j,nbnd),j=1,4)
 1006 format(5x,4i5)
      read(lfnsgm,1007) (rbnd(j,1,nbnd),rbnd(j,2,nbnd),j=1,npar)
 1007 format(f12.6,e12.5)
c
      lbnd(1,nbnd)=lbnd(1,nbnd)+na
      lbnd(2,nbnd)=lbnd(2,nbnd)+na
c
      do 105 j=npar+1,nparms
      rbnd(j,1,nbnd)=rbnd(npar,1,nbnd)
      rbnd(j,2,nbnd)=rbnd(npar,2,nbnd)
  105 continue
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1006) (lbnd(j,nbnd),j=1,4)
      write(lfnout,1007) (rbnd(j,1,nbnd),rbnd(j,2,nbnd),j=1,nparms)
      endif
c
      if(itopol.eq.0.and.ipardef.gt.1) then
      dtemp=rbnd(1,1,nbnd)
      rbnd(1,1,nbnd)=rbnd(ipardef,1,nbnd)
      rbnd(ipardef,1,nbnd)=dtemp
      dtemp=rbnd(1,2,nbnd)
      rbnd(1,2,nbnd)=rbnd(ipardef,2,nbnd)
      rbnd(ipardef,2,nbnd)=dtemp
      endif
c
    5 continue
c
c     read the angle list
c     -------------------
c
      do 6 i=1,nsang
      nang=nang+1
      read(lfnsgm,1008) (lang(j,nang),j=1,5)
 1008 format(5x,5i5)
      if(ffield(1:6).ne.'charmm') then
      read(lfnsgm,1009) (rang(j,1,nang),rang(j,2,nang),j=1,npar)
 1009 format(f10.6,e12.5)
      else
      read(lfnsgm,1019) (rang(j,1,nang),rang(j,2,nang),
     + rang(j,3,nang),rang(j,4,nang),j=1,npar)
 1019 format(2(f10.6,e12.5))
      endif
      lang(1,nang)=lang(1,nang)+na
      lang(2,nang)=lang(2,nang)+na
      lang(3,nang)=lang(3,nang)+na
c
      do 106 j=npar+1,nparms
      rang(j,1,nbnd)=rang(npar,1,nang)
      rang(j,2,nbnd)=rang(npar,2,nang)
      if(ffield(1:6).eq.'charmm') then
      rang(j,3,nbnd)=rang(npar,3,nang)
      rang(j,4,nbnd)=rang(npar,4,nang)
      endif
  106 continue
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1008) (lang(j,nang),j=1,5)
      if(ffield(1:6).ne.'charmm') then
      write(lfnout,1009) (rang(j,1,nang),rang(j,2,nang),j=1,nparms)
      else
      write(lfnout,1019) (rang(j,1,nang),rang(j,2,nang),
     + rang(j,3,nang),rang(j,4,nang),j=1,nparms)
      endif
      endif
c
      if(itopol.eq.0.and.ipardef.gt.1) then
      dtemp=rang(1,1,nang)
      rang(1,1,nang)=rang(ipardef,1,nang)
      rang(ipardef,1,nang)=dtemp
      dtemp=rang(1,2,nang)
      rang(1,2,nang)=rang(ipardef,2,nang)
      rang(ipardef,2,nang)=dtemp
      if(ffield(1:6).eq.'charmm') then
      dtemp=rang(1,1,nang)
      rang(1,3,nang)=rang(ipardef,3,nang)
      rang(ipardef,3,nang)=dtemp
      dtemp=rang(1,4,nang)
      rang(1,4,nang)=rang(ipardef,4,nang)
      rang(ipardef,4,nang)=dtemp
      endif
      endif
c
    6 continue
c
c     read the torsion list
c     ---------------------
c
      do 7 i=1,nsdih
      read(lfnsgm,1010) (ld(j),j=1,6)
 1010 format(5x,6i5)
      read(lfnsgm,1011) (md(j),rd(j,1),rd(j,2),j=1,npar)
 1011 format(i3,f10.6,e12.5)
c
      if(itopol.eq.0.and.ipardef.gt.1) then
      itemp=md(1)
      md(1)=md(ipardef)
      md(ipardef)=itemp
      dtemp=rd(1,1)
      rd(1,1)=rd(ipardef,1)
      rd(ipardef,1)=dtemp
      dtemp=rd(1,2)
      rd(1,2)=rd(ipardef,2)
      rd(ipardef,2)=dtemp
      endif
c
      l=0
      do 8 j=1,npar
      if(md(j).ge.0) l=l+1
    8 continue
      if(l.eq.npar) then
      ndih=ndih+1
      do 9 j=1,nparms
      nfdih(j,ndih)=1
    9 continue
      do 10 j=1,6
      ldih(j,ndih)=ld(j)
   10 continue
      else
      do 11 j=1,npar
      if(md(j).lt.0) nfdih(j,ndih)=nfdih(j,ndih)+1
   11 continue
      endif
c
      do 12 l=1,nparms
      do 13 j=1,nparms
      rdih(nfdih(l,ndih),j,1,ndih)=rd(j,1)
      rdih(nfdih(l,ndih),j,2,ndih)=rd(j,2)
      kdih(nfdih(l,ndih),j,ndih)=md(j)
   13 continue
   12 continue
c
      ldih(1,ndih)=ldih(1,ndih)+na
      ldih(2,ndih)=ldih(2,ndih)+na
      ldih(3,ndih)=ldih(3,ndih)+na
      ldih(4,ndih)=ldih(4,ndih)+na
c
      do 107 j=npar+1,nparms
      nfdih(j,ndih)=nfdih(npar,ndih)
      do 117 k=1,nfdih(j,ndih)
      kdih(k,j,ndih)=kdih(k,npar,ndih)
      rdih(k,j,1,ndih)=rdih(k,npar,1,ndih)
      rdih(k,j,2,ndih)=rdih(k,npar,2,ndih)
  117 continue
  107 continue
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1010) (ldih(j,ndih),j=1,6)
      write(lfnout,1011) ((kdih(l,j,ndih),rdih(l,j,1,ndih),
     + rdih(l,j,2,ndih),l=1,nfdih(j,ndih)),j=1,nparms)
      endif
c
    7 continue
c
c     read the improper torsion list
c     ------------------------------
c
      do 14 i=1,nsimp
      nimp=nimp+1
      if(ffield(1:5).eq.'amber') then
      read(lfnsgm,1012) limp(2,nimp),limp(3,nimp),limp(1,nimp),
     + (limp(j,nimp),j=4,6)
      else
      read(lfnsgm,1012) (limp(j,nimp),j=1,6)
      endif
 1012 format(5x,6i5)
      read(lfnsgm,1013) (kimp(j,nimp),rimp(j,1,nimp),
     + rimp(j,2,nimp),j=1,npar)
 1013 format(i3,f10.6,e12.5)
c
      limp(1,nimp)=limp(1,nimp)+na
      limp(2,nimp)=limp(2,nimp)+na
      limp(3,nimp)=limp(3,nimp)+na
      limp(4,nimp)=limp(4,nimp)+na
c
      do 114 j=npar+1,nparms
      kimp(j,nimp)=kimp(npar,nimp)
      rimp(j,1,nimp)=rimp(npar,1,nimp)
      rimp(j,2,nimp)=rimp(npar,2,nimp)
  114 continue
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1012) (limp(j,nimp),j=1,6)
      write(lfnout,1013) (kimp(j,nimp),rimp(j,1,nimp),
     + rimp(j,2,nimp),j=1,nparms)
      endif
c
      if(itopol.eq.0.and.ipardef.gt.1) then
      itemp=kimp(1,nimp)
      kimp(1,nimp)=kimp(ipardef,nimp)
      kimp(ipardef,nimp)=itemp
      dtemp=rimp(1,1,nimp)
      rimp(1,1,nimp)=rimp(ipardef,1,nimp)
      rimp(ipardef,1,nimp)=dtemp
      dtemp=rimp(1,2,nimp)
      rimp(1,2,nimp)=rimp(ipardef,2,nimp)
      rimp(ipardef,2,nimp)=dtemp
      endif
c
   14 continue
c
c     copy the z-matrix information
c     -----------------------------
c
      if(nzmat.gt.0) then
      open(unit=lfnmat,file=filmat(1:index(filmat,' ')-1),
     + form='formatted',status='unknown',err=9999)
      do 15 i=1,nzmat
      read(lfnsgm,1014) izm,zm
      write(lfnmat,1015) izm,zm
 1014 format(5x,4i5,3f12.6)
 1015 format(4i5,3f12.6)
   15 continue
      close(unit=lfnmat)
      endif
c
      close(unit=lfnsgm)
c
      if(util_print('where',print_debug)) then
      write(lfnout,1120) filsgm(1:length)
 1120 format('READING SEGMENT FILE ',a,' DONE')
      endif
c
      argos_prep_rdsgm=.true.
      return
 9999 continue
      argos_prep_rdsgm=.false.
      return
      end
 
