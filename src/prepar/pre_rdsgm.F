      logical function pre_rdsgm(lfnout,ffield,imol,isgm,igrp,
     + ipgrp,lfnsgm,filsgm,latm,catm,qatm,matm,natm,lbnd,rbnd,mbnd,nbnd,
     + lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,limp,rimp,mimp,nimp,
     + wcorr)
c
c $Id: pre_rdsgm.F,v 1.12 1999-05-10 22:38:06 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfnsgm,imol,isgm,igrp,ipgrp
      character*255 filsgm
      character*80 ffield
      integer matm,natm
      integer latm(9,matm)
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
      real*8 wcorr(3)
c
      character*80 card
      integer i,j,l,length,na,jmol
      integer nsatm,nsbnd,nsang,nsdih,nsimp,ld(6),md(3)
      real*8 rd(6)
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
      read(card,1001,err=9999) nsatm,nsbnd,nsang,nsdih,nsimp,wcorr
 1001 format(5i5,3f12.6)
      if(natm+nsatm.gt.matm) call errquit('increase matm',9999)
      if(nbnd+nsbnd.gt.mbnd) call errquit('increase mbnd',9999)
      if(nang+nsang.gt.mang) call errquit('increase mang',9999)
      if(ndih+nsdih.gt.mdih) call errquit('increase mdih',9998)
      if(nimp+nsimp.gt.mimp) call errquit('increase mimp',9999)
c
c     read the atom list
c     ------------------
c
      na=natm
      do 2 i=1,nsatm
      natm=natm+1
      read(lfnsgm,1002) (catm(j,natm),j=1,4),(latm(j,natm),j=1,4)
 1002 format(5x,4a6,4i5)
      read(lfnsgm,1003) (qatm(j,natm),j=1,6)
 1003 format(6f12.6)
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1002) (catm(j,natm),j=1,4),(latm(j,natm),j=1,4)
      write(lfnout,1003) (qatm(j,natm),j=1,6)
      endif
c
      latm(5,natm)=isgm
      latm(6,natm)=jmol
      jmol=iabs(jmol)
      latm(1,natm)=latm(1,natm)+igrp
      latm(2,natm)=latm(2,natm)+ipgrp
c
      if(isgm.eq.0) then
      catm(2,natm)(6:6)='w'
      catm(3,natm)(6:6)='w'
      catm(4,natm)(6:6)='w'
      endif
c
    2 continue
      igrp=latm(1,natm)
      ipgrp=latm(2,natm)
c
c     read the bond list
c     ------------------
c
      do 3 i=1,nsbnd
      nbnd=nbnd+1
      read(lfnsgm,1004) (lbnd(j,nbnd),j=1,4)
 1004 format(5x,4i5)
      read(lfnsgm,1005) (rbnd(j,nbnd),j=1,6)
 1005 format(3(f12.6,e12.5))
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1004) (lbnd(j,nbnd),j=1,4)
      write(lfnout,1005) (rbnd(j,nbnd),j=1,6) 
      endif
c
      lbnd(1,nbnd)=lbnd(1,nbnd)+na
      lbnd(2,nbnd)=lbnd(2,nbnd)+na
    3 continue
c
c     read the angle list
c     -------------------
c
      do 4 i=1,nsang
      nang=nang+1
      read(lfnsgm,1006) (lang(j,nang),j=1,5)
 1006 format(5x,5i5)
      read(lfnsgm,1007) (rang(j,nang),j=1,6)
 1007 format(3(f10.6,e12.5))
      lang(1,nang)=lang(1,nang)+na
      lang(2,nang)=lang(2,nang)+na
      lang(3,nang)=lang(3,nang)+na
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1006) (lang(j,nang),j=1,5)
      write(lfnout,1007) (rang(j,nang),j=1,6)
      endif
    4 continue
c
c     read the torsion list
c     ---------------------
c
      do 5 i=1,nsdih
      read(lfnsgm,1008) (ld(j),j=1,6)
 1008 format(5x,6i5)
      read(lfnsgm,1009) (md(j),rd(2*j-1),rd(2*j),j=1,3)
 1009 format(3(i3,f10.6,e12.5))
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1008) (ld(j),j=1,6)
      write(lfnout,1009) (md(j),rd(2*j-1),rd(2*j),j=1,3)
      endif
c
      if(md(1).ge.0.and.md(2).ge.0.and.md(3).ge.0) then
      ndih=ndih+1
      ldih(7,ndih)=1
      ldih(8,ndih)=1
      ldih(9,ndih)=1
      do 51 j=1,6
      ldih(j,ndih)=ld(j)
   51 continue
      else
      if(md(1).lt.0) ldih(7,ndih)=ldih(7,ndih)+1
      if(md(2).lt.0) ldih(8,ndih)=ldih(8,ndih)+1
      if(md(3).lt.0) ldih(9,ndih)=ldih(9,ndih)+1
      endif
c
      do 150 l=1,3
      do 52 j=1,6
      rdih(ldih(6+l,ndih),j,ndih)=rd(j)
   52 continue
      do 53 j=1,3
      kdih(ldih(6+l,ndih),j,ndih)=md(j)
   53 continue
  150 continue
c
      ldih(1,ndih)=ldih(1,ndih)+na
      ldih(2,ndih)=ldih(2,ndih)+na
      ldih(3,ndih)=ldih(3,ndih)+na
      ldih(4,ndih)=ldih(4,ndih)+na
    5 continue
c
c     read the improper torsion list
c     ------------------------------
c
      do 6 i=1,nsimp
      nimp=nimp+1
      if(ffield(1:5).eq.'amber') then
      read(lfnsgm,1010) limp(2,nimp),limp(3,nimp),limp(1,nimp),
     + (limp(j,nimp),j=4,6)
      else
      read(lfnsgm,1010) (limp(j,nimp),j=1,6)
      endif
 1010 format(5x,6i5)
      read(lfnsgm,1011) (limp(6+j,nimp),rimp(2*j-1,nimp),
     + rimp(2*j,nimp),j=1,3)
 1011 format(3(i3,f10.6,e12.5))
c
      if(util_print('connectivity',print_debug)) then
      write(lfnout,1010) (limp(j,nimp),j=1,6)
      write(lfnout,1011) (limp(6+j,nimp),rimp(2*j-1,nimp),
     + rimp(2*j,nimp),j=1,3)
      endif
c
      limp(1,nimp)=limp(1,nimp)+na
      limp(2,nimp)=limp(2,nimp)+na
      limp(3,nimp)=limp(3,nimp)+na
      limp(4,nimp)=limp(4,nimp)+na
c
    6 continue
c
      close(unit=lfnsgm)
      pre_rdsgm=.true.
      return
 9999 continue
      pre_rdsgm=.false.
      return
      end

