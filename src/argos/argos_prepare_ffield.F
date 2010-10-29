      logical function argos_prepare_ffield(idir,lfnpar,filep,lfnout,
     + q14fac,releps,
     + icomb,latt,catt,patt,ratt,matt,natt,latm,catm,matm,natm,
     + lbnd,rbnd,mbnd,nbnd,lang,rang,mang,nang,
     + ldih,nfdih,kdih,rdih,mdih,ndih,
     + limp,kimp,rimp,mimp,nimp,nparms,mparms)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "argos_prepare_common.fh"
c
c
      integer idir,lfnpar,lfnout,nparms,mparms
      character*255 filep
      real*8 q14fac,releps
      integer matt,natt,icomb
      integer latt(3,matt)
      character*6 catt(2,matt)
      real*8 patt(4,2,matt,matt),ratt(matt)
      integer matm,natm
      integer latm(11,matm)
      character*6 catm(mparms,matm)
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
      integer limp(6,mimp),kimp(nparms,mdih)
      real*8 rimp(nparms,2,mimp)
c
      integer i,j,length,nt,ng,latnum,mult,ifnd,nd,jcomb
      character*80 card
      character*6 atyp1,atyp2,atyp3,atyp4,btyp1,btyp2,btyp3,btyp4
      character*1 qtm
      real*8 rmass,p1,p2,p3,p4,pt1,pt2,pt3,pt4
c
      integer mask(3)
c
      mask(1)=4
      mask(2)=2
      mask(3)=1
c
      length=index(filep,' ')-1
      open(unit=lfnpar,file=filep(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('files',print_default)) then
      write(lfnout,1111) filep(1:length)
 1111 format(' Parameter file',t40,a)
      endif
c
      if(icomb.eq.0) icomb=3
c
    1 continue
      read(lfnpar,1000,err=9999,end=999) card
 1000 format(a)
      if(card(1:1).eq.'#') goto 1
      if(card(1:3).eq.'End') goto 999
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
c      if(icomb.ne.0.and.icomb.ne.1)
c     + call md_abort('Combination rules should not be mixed',9999)
      jcomb=1
      if(util_print('files',print_default).and.icomb.ne.jcomb) then
      if(icomb.eq.2) write(lfnout,1102)
      if(icomb.eq.3) write(lfnout,1103)
 1101 format(/,' Parameters converted to C6 and C12')
 1102 format(/,' Parameters converted to epsilon and sigma')
 1103 format(/,' Parameters converted to epsilon and R*')
      endif
      goto 1
      endif
      if(card(1:24).eq.'Parameters epsilon sigma') then
c      if(icomb.ne.0.and.icomb.ne.2)
c     + call md_abort('Combination rules should not be mixed',9999)
      jcomb=2
      if(util_print('files',print_default).and.icomb.ne.jcomb) then
      if(icomb.eq.1) write(lfnout,1101)
      if(icomb.eq.3) write(lfnout,1103)
      endif
      goto 1
      endif
      if(card(1:21).eq.'Parameters epsilon R*') then
c      if(icomb.ne.0.and.icomb.ne.3)
c     + call md_abort('Combination rules should not be mixed',9999)
      jcomb=3
      if(util_print('files',print_default).and.icomb.ne.jcomb) then
      if(icomb.eq.1) write(lfnout,1101)
      if(icomb.eq.2) write(lfnout,1102)
      endif
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
      if(card(1:3).eq.'End') goto 999
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
c
      if(icomb.ne.jcomb) then
      if(icomb.eq.1) then
      if(jcomb.eq.2) then
      call md_abort('Parameter conversion not implemented',12)
      elseif(jcomb.eq.3) then
      call md_abort('Parameter conversion not implemented',13)
      elseif(jcomb.eq.4) then
      call md_abort('Parameter conversion not implemented',14)
      endif
      elseif(icomb.eq.2) then
      if(jcomb.eq.1) then
      call md_abort('Parameter conversion not implemented',21)
      elseif(jcomb.eq.3) then
      call md_abort('Parameter conversion not implemented',23)
      elseif(jcomb.eq.4) then
      call md_abort('Parameter conversion not implemented',24)
      endif
      elseif(icomb.eq.3) then
      if(jcomb.eq.1) then
      patt(1,1,i,i)=0.25d0*(p1*p1/p2)
      patt(2,1,i,i)=0.5d0*((2.0d0*p2/p1)**(1.0d0/6.0d0))
      patt(3,1,i,i)=p3
      patt(4,1,i,i)=p4
      patt(1,2,i,i)=0.25d0*(pt1*pt1/pt2)
      patt(2,2,i,i)=0.5d0*((2.0d0*pt2/pt1)**(1.0d0/6.0d0))
      patt(3,2,i,i)=pt3
      patt(4,2,i,i)=pt4
      elseif(jcomb.eq.2) then
      call md_abort('Parameter conversion not implemented',32)
      elseif(jcomb.eq.4) then
      call md_abort('Parameter conversion not implemented',34)
      endif
      elseif(icomb.eq.4) then
      if(jcomb.eq.1) then
      call md_abort('Parameter conversion not implemented',41)
      elseif(jcomb.eq.2) then
      call md_abort('Parameter conversion not implemented',42)
      elseif(jcomb.eq.3) then
      call md_abort('Parameter conversion not implemented',43)
      endif
      endif
      else
      patt(1,1,i,i)=p1
      patt(2,1,i,i)=p2
      patt(3,1,i,i)=p3
      patt(4,1,i,i)=p4
      patt(1,2,i,i)=pt1
      patt(2,2,i,i)=pt2
      patt(3,2,i,i)=pt3
      patt(4,2,i,i)=pt4
      endif
c
      do 5 j=1,natt
      if(i.ne.j) then
      if(icomb.eq.1) then
      endif
      if(icomb.eq.2) then
      patt(1,1,i,j)=sqrt(patt(1,1,i,i)*patt(1,1,j,j))
      patt(1,1,j,i)=patt(1,1,i,j)
      patt(2,1,i,j)=patt(2,1,i,i)+patt(2,1,j,j)
      patt(2,1,j,i)=patt(2,1,i,j)
      patt(1,2,i,j)=sqrt(patt(1,2,i,i)*patt(1,2,j,j))
      patt(1,2,j,i)=patt(1,2,i,j)
      patt(2,2,i,j)=patt(2,2,i,i)+patt(2,2,j,j)
      patt(2,2,j,i)=patt(2,2,i,j)
      endif
      if(icomb.eq.3) then
      patt(1,1,i,j)=sqrt(patt(1,1,i,i)*patt(1,1,j,j))
      patt(1,1,j,i)=patt(1,1,i,j)
      patt(2,1,i,j)=patt(2,1,i,i)+patt(2,1,j,j)
      patt(2,1,j,i)=patt(2,1,i,j)
      patt(1,2,i,j)=sqrt(patt(1,2,i,i)*patt(1,2,j,j))
      patt(1,2,j,i)=patt(1,2,i,j)
      patt(2,2,i,j)=patt(2,2,i,i)+patt(2,2,j,j)
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
      if(card(1:3).eq.'End') goto 999
      read(card,1003,err=2) atyp1(1:5),atyp2(1:5),p1,p2
 1003 format(a5,1x,a5,2x,2e12.5)
      if(atyp1.eq.'      ') goto 1
      if(atyp2.eq.'      ') goto 1
      read(lfnpar,1004,err=9999) pt1,pt2
 1004 format(13x,2e12.5)
c
      do 7 i=1,natt
      if(catt(1,i)(1:5).eq.atyp1) then
      do 8 j=1,natt
      if(catt(1,i)(1:5).eq.atyp2) then
c
      if(icomb.ne.jcomb) then
      if(icomb.eq.1) then
      if(jcomb.eq.2) then
      call md_abort('Parameter conversion not implemented',12)
      elseif(jcomb.eq.3) then
      call md_abort('Parameter conversion not implemented',13)
      elseif(jcomb.eq.4) then
      call md_abort('Parameter conversion not implemented',14)
      endif
      elseif(icomb.eq.2) then
      if(jcomb.eq.1) then
      call md_abort('Parameter conversion not implemented',21)
      elseif(jcomb.eq.3) then
      call md_abort('Parameter conversion not implemented',23)
      elseif(jcomb.eq.4) then
      call md_abort('Parameter conversion not implemented',24)
      endif
      elseif(icomb.eq.3) then
      if(jcomb.eq.1) then
      patt(1,1,i,j)=0.25d0*(p1*p1/p2)
      patt(2,1,i,j)=0.5d0*((2.0d0*p2/p1)**(1.0d0/6.0d0))
      patt(3,1,i,j)=p3
      patt(4,1,i,j)=p4
      patt(1,2,i,j)=0.25d0*(pt1*pt1/pt2)
      patt(2,2,i,j)=0.5d0*((2.0d0*pt2/pt1)**(1.0d0/6.0d0))
      patt(3,2,i,j)=pt3
      patt(4,2,i,j)=pt4
      elseif(jcomb.eq.2) then
      call md_abort('Parameter conversion not implemented',32)
      elseif(jcomb.eq.4) then
      call md_abort('Parameter conversion not implemented',34)
      endif
      elseif(icomb.eq.4) then
      if(jcomb.eq.1) then
      call md_abort('Parameter conversion not implemented',41)
      elseif(jcomb.eq.2) then
      call md_abort('Parameter conversion not implemented',42)
      elseif(jcomb.eq.3) then
      call md_abort('Parameter conversion not implemented',43)
      endif
      endif
      else
      if(icomb.eq.1) then
      patt(1,1,i,j)=p1
      patt(2,1,i,j)=p2
      patt(3,1,i,j)=p3
      patt(4,1,i,j)=p4
      patt(1,2,i,j)=pt1
      patt(2,2,i,j)=pt2
      patt(3,2,i,j)=pt3
      patt(4,2,i,j)=pt4
      elseif(icomb.eq.2) then
      patt(1,1,i,j)=p1
      patt(2,1,i,j)=p2+p2
      patt(3,1,i,j)=p3
      patt(4,1,i,j)=p4
      patt(1,2,i,j)=pt1
      patt(2,2,i,j)=pt2+pt2
      patt(3,2,i,j)=pt3
      patt(4,2,i,j)=pt4
      elseif(icomb.eq.3) then
      patt(1,1,i,j)=p1
      patt(2,1,i,j)=p2+p2
      patt(3,1,i,j)=p3
      patt(4,1,i,j)=p4
      patt(1,2,i,j)=pt1
      patt(2,2,i,j)=pt2+pt2
      patt(3,2,i,j)=pt3
      patt(4,2,i,j)=pt4
      elseif(icomb.eq.4) then
      patt(1,1,i,j)=p1
      patt(2,1,i,j)=p2
      patt(3,1,i,j)=p3
      patt(4,1,i,j)=p4
      patt(1,2,i,j)=pt1
      patt(2,2,i,j)=pt2
      patt(3,2,i,j)=pt3
      patt(4,2,i,j)=pt4
      endif
      endif
c
c      if(icomb.eq.1) then
c      patt(1,1,i,j)=p1
c      patt(2,1,i,j)=p2
c      patt(1,2,i,j)=pt1
c      patt(2,2,i,j)=pt2
c      endif
c      if(icomb.eq.2) then
c      endif
c      if(icomb.eq.3) then
c      patt(1,1,i,j)=sqrt(p1*p1)
c      patt(2,1,i,j)=p2+p2
c      patt(1,2,i,j)=sqrt(pt1*pt1)
c      patt(2,2,i,j)=pt2+pt2
c      endif
c
      patt(1,1,j,i)=patt(1,1,i,j)
      patt(2,1,j,i)=patt(2,1,i,j)
      patt(1,2,j,i)=patt(1,2,i,j)
      patt(2,2,j,i)=patt(2,2,i,j)
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
      if(card(1:3).eq.'End') goto 999
      read(card,1005,err=2) atyp1(1:5),atyp2(1:5),p1,p2
 1005 format(a5,1x,a5,f10.5,e12.5)
      if(atyp1.eq.'      ') goto 1
      if(atyp2.eq.'      ') goto 1
      do 10 i=1,nbnd
c      if(lbnd(3,i).le.0.and.lbnd(4,i).eq.0) then
      if(lbnd(3,i).le.0) then
      do 11 j=1,nparms
      if(iand(mask(j),lbnd(4,i)).eq.0) then
      btyp1(1:5)=catm(1+j,lbnd(1,i))(1:5)
      btyp2(1:5)=catm(1+j,lbnd(2,i))(1:5)
      if((atyp1.eq.btyp1.and.atyp2.eq.btyp2).or.
     + (atyp2.eq.btyp1.and.atyp1.eq.btyp2)) then
      rbnd(j,1,i)=p1
      rbnd(j,2,i)=p2
      lbnd(3,i)=-(ior(iabs(lbnd(3,i)),2**(j-1)))
      endif
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
      if(card(1:3).eq.'End') goto 999
      if(ffield(1:6).ne.'charmm') then
      read(card,1006,err=2) atyp1(1:5),atyp2(1:5),atyp3(1:5),p1,p2
      else
      read(card,1006,err=2) atyp1(1:5),atyp2(1:5),atyp3(1:5),p1,p2,p3,p4
      endif
 1006 format(a5,1x,a5,1x,a5,2(f10.5,e12.5))
      if(atyp1.eq.'      ') goto 1
      if(atyp2.eq.'      ') goto 1
      if(atyp3.eq.'      ') goto 1
      do 13 i=1,nang
c      if(lang(4,i).le.0.and.lang(5,i).eq.0) then
      if(lang(4,i).le.0) then
      do 14 j=1,nparms
      if(iand(mask(j),lang(5,i)).eq.0) then
      btyp1(1:5)=catm(1+j,lang(1,i))(1:5)
      btyp2(1:5)=catm(1+j,lang(2,i))(1:5)
      btyp3(1:5)=catm(1+j,lang(3,i))(1:5)
      if(atyp2.eq.btyp2.and.((atyp1.eq.btyp1.and.atyp3.eq.btyp3).or.
     + (atyp3.eq.btyp1.and.atyp1.eq.btyp3))) then
      rang(j,1,i)=p1
      rang(j,2,i)=p2
      if(ffield(1:6).eq.'charmm') then
      rang(j,3,i)=p3
      rang(j,4,i)=p4
      endif
      lang(4,i)=-(ior(iabs(lang(4,i)),2**(j-1)))
c      write(*,'(6i5,6a6)') lang(1,i),lang(2,i),lang(3,i),lang(4,i),
c     + lang(5,i),j,atyp1,atyp2,atyp3,btyp1,btyp2,btyp3
c      write(*,'(4i5,3a6)') (lang(j,i),j=1,4),atyp1,atyp2,atyp3
      endif
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
      if(card(1:3).eq.'End') goto 999
      read(card,1007,err=2) atyp1(1:5),atyp2(1:5),atyp3(1:5),atyp4(1:5),
     + p1,p2,mult
 1007 format(a5,1x,a5,1x,a5,1x,a5,f10.5,e12.5,i5)
      nd=nd+1
      if(atyp2.eq.'      ') goto 1
      if(atyp3.eq.'      ') goto 1
      do 16 i=1,ndih
c      if(ldih(5,i).le.0.and.ldih(6,i).eq.0) then
      if(ldih(5,i).le.0) then
      do 17 j=1,nparms
      if(iand(mask(j),ldih(6,i)).eq.0) then
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
      nfdih(j,i)=nd
      kdih(nfdih(j,i),j,i)=mult
      rdih(nfdih(j,i),j,1,i)=p1
      rdih(nfdih(j,i),j,2,i)=p2
      ldih(5,i)=-(ior(iabs(ldih(5,i)),2**(j-1)))
      endif
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
      if(card(1:3).eq.'End') goto 999
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
c      if(limp(5,i).le.0.and.limp(6,i).eq.0) then
      if(limp(5,i).le.0) then
      do 20 j=1,nparms
      if(iand(mask(j),limp(6,i)).eq.0) then
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
      kimp(j,i)=mult
      rimp(j,1,i)=ifnd*p1
      rimp(j,2,i)=p2
      limp(5,i)=-(ior(iabs(limp(5,i)),2**(j-1)))
      endif
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
      argos_prepare_ffield=.true.
      if(util_print('where',print_debug)) then
      write(lfnout,2006)
 2006 format(/,'argos_prepare_ffield done')
      endif
      return
 9999 continue
      argos_prepare_ffield=.false.
      return
      end
 
