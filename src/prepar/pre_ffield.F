      logical function pre_ffield(idir,lfnpar,filpar,lfnout,
     + q14fac,releps,ffield,
     + icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,ldih,kdih,rdih,mdih,ndih,
     + limp,rimp,mimp,nimp)
c
c $Id: pre_ffield.F,v 1.17 2000-04-12 17:21:25 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer idir,lfnpar,lfnout
      character*255 filpar
      real*8 q14fac,releps
      character*80 ffield
      integer matt,natt,icomb
      integer latt(3,matt)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(9,matm)
      character*6 catm(4,matm)
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
      integer i,j,length,nt,ng,latnum,mult,ifnd,nd
      character*80 card
      character*6 atyp1,atyp2,atyp3,atyp4,btyp1,btyp2,btyp3,btyp4
      character*1 qtm
      real*8 rmass,p1,p2,p3,p4,pt1,pt2,pt3,pt4
c
      length=index(filpar,' ')-1
      open(unit=lfnpar,file=filpar(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('files',print_default)) then
      write(lfnout,1111) filpar(1:length)
 1111 format(' Parameter file',t40,a)
      endif
c
      icomb=0
c
    1 continue
      read(lfnpar,1000,err=9999,end=999) card
 1000 format(a)
      if(card(1:1).eq.'#') goto 1
    2 continue
c
c     1-4 electrostatic scaling
c     -------------------------
c
      if(card(1:32).eq.'Electrostatic 1-4 scaling factor') then
      read(card(34:45),'(f12.6)') q14fac
      goto 1
      endif
c
c     relative dielectric constant
c     ----------------------------
c
      if(card(1:28).eq.'Relative dielectric constant') then
      read(card(30:41),'(f12.6)') releps
      goto 1
      endif
c
c     non-bonded parameter format
c     ---------------------------
c
      if(card(1:17).eq.'Parameters C6 C12') then
      if(icomb.ne.0.and.icomb.ne.1)
     + call errquit('Combination rules should not be mixed',9999)
      icomb=1
      goto 1
      endif
      if(card(1:24).eq.'Parameters epsilon sigma') then
      if(icomb.ne.0.and.icomb.ne.2)
     + call errquit('Combination rules should not be mixed',9999)
      icomb=2
      goto 1
      endif
      if(card(1:21).eq.'Parameters epsilon R*') then
      if(icomb.ne.0.and.icomb.ne.3)
     + call errquit('Combination rules should not be mixed',9999)
      icomb=3
      goto 1
      endif
c
c     atom types
c     ----------
c
      if(card(1:5).eq.'Atoms') then
      atyp1(6:6)=' '
      atyp2(6:6)=' '
    3 continue
      read(lfnpar,1000,end=999) card
      if(card(1:1).eq.'#') goto 3
      if(card(1:5).eq.'Cross') goto 2
      if(card(1:5).eq.'Bonds') goto 2
      if(card(1:6).eq.'Angles') goto 2
      if(card(1:16).eq.'Proper dihedrals') goto 2
      if(card(1:18).eq.'Improper dihedrals') goto 2
      if(card(1:10).eq.'Atom types') goto 2
      read(card,1001,err=2) atyp1(1:5),rmass,p1,p2,p3,p4,nt,ng
 1001 format(a5,f10.5,4e12.5,i5,1x,i10)
      if(atyp1.eq.'      ') goto 1
      read(lfnpar,1002,err=9999) atyp2(1:5),qtm,latnum,pt1,pt2,pt3,pt4
 1002 format(a5,1x,a1,3x,i5,4e12.5)
      if(atyp2.eq.'      ') atyp2=atyp1
      do 4 i=1,natt
      if(catt(1,i)(1:5).eq.atyp1(1:5).and.
     + (qtm.eq.' '.or.qtm.eq.catt(1,i)(6:6))) then
      catt(2,i)=atyp2
      ratt(i)=rmass
      latt(1,i)=latnum
      latt(2,i)=nt
      latt(3,i)=ng
      patt(1,1,i,i)=p1
      patt(2,1,i,i)=p2
      patt(3,1,i,i)=p3
      patt(4,1,i,i)=p4
      patt(1,2,i,i)=pt1
      patt(2,2,i,i)=pt2
      patt(3,2,i,i)=pt3
      patt(4,2,i,i)=pt4
      do 5 j=1,natt
      if(i.ne.j) then
      if(icomb.eq.1) then
      endif
      if(icomb.eq.2) then
      patt(1,1,i,j)=sqrt(p1*patt(1,1,j,j))
      patt(1,1,j,i)=patt(1,1,i,j)
      patt(2,1,i,j)=p2+patt(2,1,j,j)
      patt(2,1,j,i)=patt(2,1,i,j)
      patt(1,2,i,j)=sqrt(pt1*patt(1,2,j,j))
      patt(1,2,j,i)=patt(1,2,i,j)
      patt(2,2,i,j)=pt2+patt(2,2,j,j)
      patt(2,2,j,i)=patt(2,2,i,j)
      endif
      if(icomb.eq.3) then
      patt(1,1,i,j)=sqrt(p1*patt(1,1,j,j))
      patt(1,1,j,i)=patt(1,1,i,j)
      patt(2,1,i,j)=p2+patt(2,1,j,j)
      patt(2,1,j,i)=patt(2,1,i,j)
      patt(1,2,i,j)=sqrt(pt1*patt(1,2,j,j))
      patt(1,2,j,i)=patt(1,2,i,j)
      patt(2,2,i,j)=pt2+patt(2,2,j,j)
      patt(2,2,j,i)=patt(2,2,i,j)
      endif
      endif
    5 continue
      endif
    4 continue
      goto 3
      endif
c
c     cross terms
c     -----------
c
      if(card(1:5).eq.'Cross') then
      atyp1(6:6)=' '
      atyp2(6:6)=' '
    6 continue
      read(lfnpar,1000,end=999) card
      if(card(1:1).eq.'#') goto 6
      if(card(1:5).eq.'Atoms') goto 2
      if(card(1:5).eq.'Bonds') goto 2
      if(card(1:6).eq.'Angles') goto 2
      if(card(1:16).eq.'Proper dihedrals') goto 2
      if(card(1:18).eq.'Improper dihedrals') goto 2
      if(card(1:10).eq.'Atom types') goto 2
      read(card,1003,err=2) atyp1(1:5),atyp2(1:5),p1,p2
 1003 format(a5,1x,a5,2x,2e12.5)
      if(atyp1.eq.'      ') goto 1
      if(atyp2.eq.'      ') goto 1
      read(lfnpar,1004,err=9999) pt1,pt2
 1004 format(13x,2e12.5)
      do 7 i=1,natt
      if(catt(1,i)(1:5).eq.atyp1) then
      do 8 j=1,natt
      if(catt(1,i)(1:5).eq.atyp2) then
      if(icomb.eq.1) then
      patt(1,1,i,j)=p1
      patt(2,1,i,j)=p2
      patt(1,2,i,j)=pt1
      patt(2,2,i,j)=pt2
      endif
      if(icomb.eq.2) then
      endif
      if(icomb.eq.3) then
      patt(1,1,i,j)=p1*p1
      patt(2,1,i,j)=p2+p2
      patt(1,2,i,j)=pt1*pt1
      patt(2,2,i,j)=pt2+pt2
      endif
      patt(1,1,i,j)=patt(1,1,j,i)
      patt(2,1,i,j)=patt(2,1,j,i)
      patt(1,2,i,j)=patt(1,2,j,i)
      patt(2,2,i,j)=patt(2,2,j,i)
      endif
    8 continue
      endif
    7 continue
      goto 6
      endif
c
c     bonds
c     -----
c
      if(card(1:5).eq.'Bonds') then
      atyp1(6:6)=' '
      atyp2(6:6)=' '
      btyp1(6:6)=' '
      btyp2(6:6)=' '
    9 continue
      read(lfnpar,1000,end=999) card
      if(card(1:1).eq.'#') goto 9
      if(card(1:5).eq.'Atoms') goto 2
      if(card(1:5).eq.'Cross') goto 2
      if(card(1:6).eq.'Angles') goto 2
      if(card(1:16).eq.'Proper dihedrals') goto 2
      if(card(1:18).eq.'Improper dihedrals') goto 2
      if(card(1:10).eq.'Atom types') goto 2
      read(card,1005,err=2) atyp1(1:5),atyp2(1:5),p1,p2
 1005 format(a5,1x,a5,f10.5,e12.5)
      if(atyp1.eq.'      ') goto 1
      if(atyp2.eq.'      ') goto 1
      do 10 i=1,nbnd
      if(lbnd(3,i).le.0) then
      do 11 j=1,3
      btyp1(1:5)=catm(1+j,lbnd(1,i))(1:5)
      btyp2(1:5)=catm(1+j,lbnd(2,i))(1:5)
      if((atyp1.eq.btyp1.and.atyp2.eq.btyp2).or.
     + (atyp2.eq.btyp1.and.atyp1.eq.btyp2)) then
      rbnd(j*2-1,i)=p1
      rbnd(j*2,i)=p2
      lbnd(3,i)=-(ior(iabs(lbnd(3,i)),2**(j-1)))
      endif
   11 continue
      endif
   10 continue
      goto 9
      endif
c
c     angles
c     ------
c
      if(card(1:6).eq.'Angles') then
      atyp1(6:6)=' '
      atyp2(6:6)=' '
      atyp3(6:6)=' '
      btyp1(6:6)=' '
      btyp2(6:6)=' '
      btyp3(6:6)=' '
   12 continue
      read(lfnpar,1000,end=999) card
      if(card(1:1).eq.'#') goto 12
      if(card(1:5).eq.'Atoms') goto 2
      if(card(1:5).eq.'Cross') goto 2
      if(card(1:5).eq.'Bonds') goto 2
      if(card(1:16).eq.'Proper dihedrals') goto 2
      if(card(1:18).eq.'Improper dihedrals') goto 2
      if(card(1:10).eq.'Atom types') goto 2
      read(card,1006,err=2) atyp1(1:5),atyp2(1:5),atyp3(1:5),p1,p2
 1006 format(a5,1x,a5,1x,a5,f10.5,e12.5)
      if(atyp1.eq.'      ') goto 1
      if(atyp2.eq.'      ') goto 1
      if(atyp3.eq.'      ') goto 1
      do 13 i=1,nang
      if(lang(4,i).le.0) then
      do 14 j=1,3
      btyp1(1:5)=catm(1+j,lang(1,i))(1:5)
      btyp2(1:5)=catm(1+j,lang(2,i))(1:5)
      btyp3(1:5)=catm(1+j,lang(3,i))(1:5)
      if(atyp2.eq.btyp2.and.((atyp1.eq.btyp1.and.atyp3.eq.btyp3).or.
     + (atyp3.eq.btyp1.and.atyp1.eq.btyp3))) then
      rang(j*2-1,i)=p1
      rang(j*2,i)=p2
      lang(4,i)=-(ior(iabs(lang(4,i)),2**(j-1)))
      endif
   14 continue
      endif
   13 continue
      goto 12
      endif
c
c     torsions
c     --------
c
      if(card(1:16).eq.'Proper dihedrals') then
      atyp1(6:6)=' '
      atyp2(6:6)=' '
      atyp3(6:6)=' '
      atyp4(6:6)=' '
      btyp1(6:6)=' '
      btyp2(6:6)=' '
      btyp3(6:6)=' '
      btyp4(6:6)=' '
      nd=0
   15 continue
      read(lfnpar,1000,end=999) card
      if(card(1:1).eq.'#') goto 15
      if(card(1:5).eq.'Atoms') goto 2
      if(card(1:5).eq.'Cross') goto 2
      if(card(1:5).eq.'Bonds') goto 2
      if(card(1:6).eq.'Angles') goto 2
      if(card(1:18).eq.'Improper dihedrals') goto 2
      if(card(1:10).eq.'Atom types') goto 2
      read(card,1007,err=2) atyp1(1:5),atyp2(1:5),atyp3(1:5),atyp4(1:5),
     + p1,p2,mult
 1007 format(a5,1x,a5,1x,a5,1x,a5,f10.5,e12.5,i5)
      nd=nd+1
      if(atyp2.eq.'      ') goto 1
      if(atyp3.eq.'      ') goto 1
      do 16 i=1,ndih
      if(ldih(5,i).le.0) then
      do 17 j=1,3
      btyp1(1:5)=catm(1+j,ldih(1,i))(1:5)
      btyp2(1:5)=catm(1+j,ldih(2,i))(1:5)
      btyp3(1:5)=catm(1+j,ldih(3,i))(1:5)
      btyp4(1:5)=catm(1+j,ldih(4,i))(1:5)
c
      if(((atyp1.eq.'      '.or.atyp1.eq.btyp1).and.
     + atyp2.eq.btyp2.and.atyp3.eq.btyp3.and.
     + (atyp4.eq.'      '.or.atyp4.eq.btyp4)).or.
     + ((atyp1.eq.'      '.or.atyp1.eq.btyp4).and.
     + atyp2.eq.btyp3.and.atyp3.eq.btyp2.and.
     + (atyp4.eq.'      '.or.atyp4.eq.btyp1))) then
      ldih(6+j,i)=nd
      kdih(ldih(6+j,i),j,i)=mult
      rdih(ldih(6+j,i),2*j-1,i)=p1 
      rdih(ldih(6+j,i),2*j,i)=p2 
      ldih(5,i)=-(ior(iabs(ldih(5,i)),2**(j-1)))
      endif
   17 continue
      endif
   16 continue
      if(mult.gt.0) nd=0
      goto 15
      endif
c
c     impropers
c     ---------
c
      if(card(1:18).eq.'Improper dihedrals') then
      atyp1(6:6)=' '
      atyp2(6:6)=' '
      atyp3(6:6)=' '
      atyp4(6:6)=' '
      btyp1(6:6)=' '
      btyp2(6:6)=' '
      btyp3(6:6)=' '
      btyp4(6:6)=' '
   18 continue
      read(lfnpar,1000,end=999) card
      if(card(1:1).eq.'#') goto 18
      if(card(1:5).eq.'Atoms') goto 2
      if(card(1:5).eq.'Cross') goto 2
      if(card(1:5).eq.'Bonds') goto 2
      if(card(1:6).eq.'Angles') goto 2
      if(card(1:16).eq.'Proper dihedrals') goto 2
      if(card(1:10).eq.'Atom types') goto 2
      if(ffield(1:5).eq.'amber') then
      read(card,1008,err=2) atyp2(1:5),atyp3(1:5),atyp1(1:5),atyp4(1:5),
     + p1,p2,mult
 1008 format(a5,1x,a5,1x,a5,1x,a5,f10.5,e12.5,i5)
      else
      read(card,1008,err=2) atyp1(1:5),atyp2(1:5),atyp3(1:5),atyp4(1:5),
     + p1,p2,mult
      endif
      if(atyp1.eq.'      ') goto 1
      if(atyp4.eq.'      ') goto 1
      do 19 i=1,nimp
      if(limp(5,i).le.0) then
      do 20 j=1,3
      btyp1(1:5)=catm(1+j,limp(1,i))(1:5)
      btyp2(1:5)=catm(1+j,limp(2,i))(1:5)
      btyp3(1:5)=catm(1+j,limp(3,i))(1:5)
      btyp4(1:5)=catm(1+j,limp(4,i))(1:5)
c
      ifnd=0
      if(atyp1.eq.btyp1) then
      if(atyp2.eq.'      ') then
      if(atyp3.eq.'      ') then
      if(atyp4.eq.btyp4) ifnd=1
      if(atyp4.eq.btyp3) ifnd=1
      if(atyp4.eq.btyp2) ifnd=1
      else
      if(atyp3.eq.btyp3.and.atyp4.eq.btyp4) ifnd=1
      if(atyp3.eq.btyp4.and.atyp4.eq.btyp3) ifnd=1
      if(atyp3.eq.btyp2.and.atyp4.eq.btyp4) ifnd=1
      if(atyp3.eq.btyp4.and.atyp4.eq.btyp2) ifnd=1
      if(atyp3.eq.btyp2.and.atyp4.eq.btyp3) ifnd=1
      if(atyp3.eq.btyp3.and.atyp4.eq.btyp2) ifnd=1
      endif
      else
      if(atyp3.eq.'      ') then
      if(atyp2.eq.btyp3.and.atyp4.eq.btyp4) ifnd=1
      if(atyp2.eq.btyp4.and.atyp4.eq.btyp3) ifnd=1
      if(atyp2.eq.btyp2.and.atyp4.eq.btyp4) ifnd=1
      if(atyp2.eq.btyp4.and.atyp4.eq.btyp2) ifnd=1
      if(atyp2.eq.btyp2.and.atyp4.eq.btyp3) ifnd=1
      if(atyp2.eq.btyp3.and.atyp4.eq.btyp2) ifnd=1
      else
      if(atyp2.eq.btyp2.and.atyp3.eq.btyp3.and.atyp4.eq.btyp4) ifnd=1
      if(atyp2.eq.btyp3.and.atyp3.eq.btyp4.and.atyp4.eq.btyp2) ifnd=1
      if(atyp2.eq.btyp4.and.atyp3.eq.btyp2.and.atyp4.eq.btyp3) ifnd=1
      if(atyp2.eq.btyp3.and.atyp3.eq.btyp2.and.atyp4.eq.btyp4) ifnd=-1
      if(atyp2.eq.btyp2.and.atyp3.eq.btyp4.and.atyp4.eq.btyp3) ifnd=-1
      if(atyp2.eq.btyp4.and.atyp3.eq.btyp3.and.atyp4.eq.btyp2) ifnd=-1
      endif
      endif
      endif
c
      if(ffield(1:5).eq.'charm'.and.ifnd.eq.0) then
      if(atyp4.eq.btyp1) then
      if(atyp3.eq.'      ') then
      if(atyp2.eq.'      ') then
      if(atyp1.eq.btyp4) ifnd=1
      if(atyp1.eq.btyp3) ifnd=1
      if(atyp1.eq.btyp2) ifnd=1
      else
      if(atyp2.eq.btyp3.and.atyp1.eq.btyp4) ifnd=1
      if(atyp2.eq.btyp4.and.atyp1.eq.btyp3) ifnd=1
      if(atyp2.eq.btyp2.and.atyp1.eq.btyp4) ifnd=1
      if(atyp2.eq.btyp4.and.atyp1.eq.btyp2) ifnd=1
      if(atyp2.eq.btyp2.and.atyp1.eq.btyp3) ifnd=1
      if(atyp2.eq.btyp3.and.atyp1.eq.btyp2) ifnd=1
      endif
      else
      if(atyp2.eq.'      ') then
      if(atyp3.eq.btyp3.and.atyp1.eq.btyp4) ifnd=1
      if(atyp3.eq.btyp4.and.atyp1.eq.btyp3) ifnd=1
      if(atyp3.eq.btyp2.and.atyp1.eq.btyp4) ifnd=1
      if(atyp3.eq.btyp4.and.atyp1.eq.btyp2) ifnd=1
      if(atyp3.eq.btyp2.and.atyp1.eq.btyp3) ifnd=1
      if(atyp3.eq.btyp3.and.atyp1.eq.btyp2) ifnd=1
      else
      if(atyp3.eq.btyp2.and.atyp2.eq.btyp3.and.atyp1.eq.btyp4) ifnd=1
      if(atyp3.eq.btyp3.and.atyp2.eq.btyp4.and.atyp1.eq.btyp2) ifnd=1
      if(atyp3.eq.btyp4.and.atyp2.eq.btyp2.and.atyp1.eq.btyp3) ifnd=1
      if(atyp3.eq.btyp3.and.atyp2.eq.btyp2.and.atyp1.eq.btyp4) ifnd=-1
      if(atyp3.eq.btyp2.and.atyp2.eq.btyp4.and.atyp1.eq.btyp3) ifnd=-1
      if(atyp3.eq.btyp4.and.atyp2.eq.btyp3.and.atyp1.eq.btyp2) ifnd=-1
      endif
      endif
      endif
      endif
c
      if(ifnd.ne.0) then
      limp(6+j,i)=mult
      rimp(2*j-1,i)=ifnd*p1 
      rimp(2*j,i)=p2 
      limp(5,i)=-(ior(iabs(limp(5,i)),2**(j-1)))
      endif
c
   20 continue
      endif
   19 continue
      goto 18
      endif
c
      goto 1
  999 continue
      close(unit=lfnpar)
c
      pre_ffield=.true.
      return
 9999 continue
      pre_ffield=.false.
      return
      end

