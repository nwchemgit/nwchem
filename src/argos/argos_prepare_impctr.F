      logical function argos_prepare_impctr(lfnout,latm,matm,
     + jlo,ilo,ihi,jhi,lbnd,mbnd,nbnd,lang,mang,nang,ldih,mdih,ndih,
     + lring,aring,mring,nring3,nring4,nring5,nring6)
c
c $Id$
c
c     function to determine the center type
c
c     in  : latm(2,matm)    = atomic numbers
c                3          = number of bonds to atom
c           matm            = dimension of atom arrays
c           ilo:ihi         = atoms for which center is determined
c           lbnd(1:3,mbnd)  = bond list
c           mbnd            = dimension of bond list
c           nbnd            = length of bond list
c           lang(1:3,mang)  = angle list
c           mang            = dimension of angle list
c           nang            = length of angle list
c           ldih(1:3,mdih)  = dihedral list
c           mdih            = dimension of dihedral list
c           ndih            = length of dihedral list
c
c     out : latm(4,matm)    = atom center type
c                             0 = undetermined
c                             1 = planar
c                             2 = tetrahedral
c                             4 = aromatic ring
c
      implicit none
c
#include "util.fh"
c
      logical argos_prepare_torang
      external argos_prepare_torang
c
      integer matm,ilo,ihi,jlo,jhi,mbnd,nbnd,mang,nang,mdih,ndih
      integer lfnout
      integer latm(5,matm),lbnd(2,mbnd),lang(3,mang),ldih(4,mdih)
c
      integer mring
      integer lring(6,mring)
      logical aring(mring)
      integer nring3,nring4,nring5,nring6
      logical valid
c
      integer i,j,k,l,li,lj,nh,ic
c
      argos_prepare_impctr=.false.
c
      do 53 i=jlo,jhi
      latm(4,i)=0
   53 continue
c
c     3 membered rings
c     ----------------
c
      nring3=0
      do 29 li=1,nang
      if(latm(2,lang(1,li)).eq.1) goto 29
      if(latm(2,lang(2,li)).eq.1) goto 29
      if(latm(2,lang(3,li)).eq.1) goto 29
      do 30 lj=1,nbnd
      if((lang(1,li).eq.lbnd(1,lj).and.lang(3,li).eq.lbnd(2,lj)).or.
     + (lang(1,li).eq.lbnd(2,lj).and.lang(3,li).eq.lbnd(1,lj))) then
      nring3=nring3+1
      if(nring3.gt.mring)
     +  call md_abort('increase mring in argos_prepare_impctr',9999)
      lring(1,nring3)=lang(1,li)
      lring(2,nring3)=lang(2,li)
      lring(3,nring3)=lang(3,li)
      do 31 i=1,2
      do 32 j=i+1,3
      if(lring(i,nring3).gt.lring(j,nring3)) then
      l=lring(i,nring3)
      lring(i,nring3)=lring(j,nring3)
      lring(j,nring3)=l
      endif
   32 continue
   31 continue
      do 33 i=1,nring3-1
      do 34 j=1,3
      if(lring(j,i).ne.lring(j,nring3)) goto 33
   34 continue
      nring3=nring3-1
      goto 29
   33 continue
      if(util_print('connectivity',print_high)) then
      write(lfnout,1000) nring3,(lring(j,nring3),j=1,3)
 1000 format('Ring',i5,': atoms ',6i5)
      endif
      endif
   30 continue
   29 continue
c
c     4 membered rings
c     ----------------
c
      nring4=nring3
      do 35 li=1,nang-1
      if(latm(2,lang(1,li)).eq.1) goto 35
      if(latm(2,lang(2,li)).eq.1) goto 35
      if(latm(2,lang(3,li)).eq.1) goto 35
      do 36 lj=li+1,nang
      if((lang(1,li).eq.lang(1,lj).and.lang(3,li).eq.lang(3,lj)).or.
     + (lang(1,li).eq.lang(3,lj).and.lang(3,li).eq.lang(1,lj))) then
      nring4=nring4+1
      if(nring4.gt.mring)
     +  call md_abort('increase mring in argos_prepare_impctr',9999)
      lring(1,nring4)=lang(1,li)
      lring(2,nring4)=lang(2,li)
      lring(3,nring4)=lang(3,li)
      lring(4,nring4)=lang(2,lj)
      do 37 i=1,3
      do 38 j=i+1,4
      if(lring(i,nring4).gt.lring(j,nring4)) then
      l=lring(i,nring4)
      lring(i,nring4)=lring(j,nring4)
      lring(j,nring4)=l
      endif
   38 continue
   37 continue
      do 39 i=nring3+1,nring4-1
      do 40 j=1,4
      if(lring(j,i).ne.lring(j,nring4)) goto 39
   40 continue
      nring4=nring4-1
      goto 35
   39 continue
      if(util_print('connectivity',print_high)) then
      write(lfnout,1000) nring4,(lring(j,nring4),j=1,4)
      endif
      endif
   36 continue
   35 continue
c
c     5 membered rings
c     ----------------
c
      nring5=nring4
      do 4 li=1,ndih
      if(latm(2,ldih(1,li)).eq.1) goto 4
      if(latm(2,ldih(2,li)).eq.1) goto 4
      if(latm(2,ldih(3,li)).eq.1) goto 4
      if(latm(2,ldih(4,li)).eq.1) goto 4
      do 5 lj=1,nang
      if((ldih(1,li).eq.lang(1,lj).and.ldih(4,li).eq.lang(3,lj)).or.
     + (ldih(1,li).eq.lang(3,lj).and.ldih(4,li).eq.lang(1,lj))) then
      nring5=nring5+1
      if(nring5.gt.mring)
     +  call md_abort('increase mring in argos_prepare_impctr',9999)
      do 6 i=1,4
      lring(i,nring5)=ldih(i,li)
    6 continue
      lring(5,nring5)=lang(2,lj)
      do 7 i=1,4
      do 8 j=i+1,5
      if(lring(i,nring5).gt.lring(j,nring5)) then
      l=lring(i,nring5)
      lring(i,nring5)=lring(j,nring5)
      lring(j,nring5)=l
      endif
    8 continue
    7 continue
      do 9 i=nring4+1,nring5-1
      do 10 j=1,5
      if(lring(j,i).ne.lring(j,nring5)) goto 9
   10 continue
      nring5=nring5-1
      goto 4
    9 continue
      if(util_print('connectivity',print_high)) then
      write(lfnout,1000) nring5,(lring(j,nring5),j=1,5)
      endif
      endif
    5 continue
    4 continue
c
c     6 membered rings
c     ----------------
c
      nring6=nring5
      do 11 li=1,ndih-1
      if(latm(2,ldih(1,li)).eq.1) goto 11
      if(latm(2,ldih(2,li)).eq.1) goto 11
      if(latm(2,ldih(3,li)).eq.1) goto 11
      if(latm(2,ldih(4,li)).eq.1) goto 11
      do 12 lj=li+1,ndih
      if((ldih(1,li).eq.ldih(1,lj).and.ldih(4,li).eq.ldih(4,lj)).or.
     + (ldih(1,li).eq.ldih(4,lj).and.ldih(4,li).eq.ldih(1,lj))) then
      nring6=nring6+1
      if(nring6.gt.mring)
     +  call md_abort('increase mring in argos_prepare_impctr',9999)
      do 13 i=1,4
      lring(i,nring6)=ldih(i,li)
   13 continue
      lring(5,nring6)=ldih(2,lj)
      lring(6,nring6)=ldih(3,lj)
      do 14 i=1,5
      do 15 j=i+1,6
      if(lring(i,nring6).gt.lring(j,nring6)) then
      l=lring(i,nring6)
      lring(i,nring6)=lring(j,nring6)
      lring(j,nring6)=l
      endif
   15 continue
   14 continue
      do 16 i=nring5+1,nring6-1
      do 17 j=1,6
      if(lring(j,i).ne.lring(j,nring6)) goto 16
   17 continue
      nring6=nring6-1
      goto 11
   16 continue
      if(util_print('connectivity',print_high)) then
      write(lfnout,1000) nring6,(lring(j,nring6),j=1,6)
      endif
      endif
   12 continue
   11 continue
c
c     select aromatic 3-membered rings
c     --------------------------------
c
      do 41 i=1,nring3
      aring(i)=.true.
      do 42 j=1,3
      valid=.false.
      l=lring(j,i)
c      if(latm(2,l).eq.6.and.latm(3,l).eq.3) valid=.true.
c      if(latm(2,l).eq.7.and.latm(3,l).eq.3) valid=.true.
c      if(latm(2,l).eq.7.and.latm(3,l).eq.2) valid=.true.
      if(.not.valid) then
      aring(i)=.false.
      goto 41
      endif
   42 continue
   41 continue
c
c     select aromatic 4-membered rings
c     --------------------------------
c
      do 43 i=nring3+1,nring4
      aring(i)=.true.
      do 44 j=1,4
      valid=.false.
      l=lring(j,i)
c      if(latm(2,l).eq.6.and.latm(3,l).eq.3) valid=.true.
c      if(latm(2,l).eq.7.and.latm(3,l).eq.3) valid=.true.
c      if(latm(2,l).eq.7.and.latm(3,l).eq.2) valid=.true.
      if(.not.valid) then
      aring(i)=.false.
      goto 43
      endif
   44 continue
   43 continue
c
c     select aromatic 5-membered rings
c     --------------------------------
c
      do 18 i=nring4+1,nring5
      aring(i)=.true.
      do 19 j=1,5
      valid=.false.
      l=lring(j,i)
      if(latm(2,l).eq.6.and.latm(3,l).eq.3) valid=.true.
      if(latm(2,l).eq.7.and.latm(3,l).eq.3) valid=.true.
      if(latm(2,l).eq.7.and.latm(3,l).eq.2) valid=.true.
      if(.not.valid) then
      aring(i)=.false.
      goto 18
      endif
   19 continue
   18 continue
c
c     select aromatic 6-membered rings
c     --------------------------------
c
      do 20 i=nring5+1,nring6
      aring(i)=.true.
      do 21 j=1,6
      valid=.false.
      l=lring(j,i)
      if(latm(2,l).eq.6.and.latm(3,l).eq.3) valid=.true.
      if(latm(2,l).eq.7.and.latm(3,l).eq.3) valid=.true.
      if(latm(2,l).eq.7.and.latm(3,l).eq.2) valid=.true.
      if(.not.valid) then
      aring(i)=.false.
      goto 20
      endif
   21 continue
   20 continue
c
c     set ring types
c     --------------
c
      do 22 i=1,nring5
      if(aring(i)) then
      do 23 j=1,5
      if(latm(4,lring(j,i)).lt.5) latm(4,lring(j,i))=0
      latm(4,lring(j,i))=latm(4,lring(j,i))+5
   23 continue
      endif
   22 continue
c
      do 24 i=nring5+1,nring6
      if(aring(i)) then
      do 25 j=1,6
      if(latm(4,lring(j,i)).lt.5) latm(4,lring(j,i))=0
      latm(4,lring(j,i))=latm(4,lring(j,i))+6
   25 continue
      endif
   24 continue
c
c     peptide bonds
c     -------------
c
      do 1 i=jlo,jhi
      if(latm(2,i).eq.7.and.latm(3,i).eq.3) then
      do 2 li=1,nbnd
      j=0
      if(lbnd(1,li).eq.i) j=lbnd(2,li)
      if(lbnd(2,li).eq.i) j=lbnd(1,li)
      if(j.ge.jlo.and.j.le.jhi) then
      if(latm(2,j).eq.6.and.latm(3,j).eq.3) then
      do 3 lj=1,nbnd
      k=0
      if(lbnd(1,lj).eq.j.and.lbnd(2,lj).ne.i) k=lbnd(2,lj)
      if(lbnd(2,lj).eq.j.and.lbnd(1,lj).ne.i) k=lbnd(1,lj)
      if(k.ge.jlo.and.k.le.jhi) then
      if(latm(2,k).eq.8.and.latm(3,k).eq.1) then
      if(latm(4,i).eq.0) latm(4,i)=1
      if(latm(4,j).eq.0) latm(4,j)=1
      goto 1
      endif
      endif
    3 continue
      endif
      endif
    2 continue
      endif
    1 continue
c
c     planar C=O
c     ----------
c
      do 26 i=jlo,jhi
      if(latm(2,i).eq.6.and.latm(3,i).eq.3.and.latm(4,i).eq.0) then
      do 27 li=1,nbnd
      j=0
      if(lbnd(1,li).eq.i) j=lbnd(2,li)
      if(lbnd(2,li).eq.i) j=lbnd(1,li)
      if(j.ge.jlo.and.j.le.jhi) then
      if(latm(2,j).eq.8.and.latm(3,j).eq.1.and.latm(4,i).eq.0)
     + latm(4,i)=1
      endif
   27 continue
      endif
   26 continue
c
c     planar NH2 in C(sp2)-NH2
c     ------------------------
c
      do 45 i=jlo,jhi
      if(latm(2,i).eq.7.and.latm(3,i).eq.3.and.latm(4,i).eq.0) then
      nh=0
      ic=0
      do 46 li=1,nbnd
      j=0
      if(lbnd(1,li).eq.i) j=lbnd(2,li)
      if(lbnd(2,li).eq.i) j=lbnd(1,li)
      if(j.ge.jlo.and.j.le.jhi) then
      if(latm(2,j).eq.1) nh=nh+1
      if(latm(2,j).eq.6.and.latm(3,j).eq.3) ic=j
      endif
   46 continue
      if((nh.eq.1.or.nh.eq.2).and.ic.gt.0.and.latm(4,i).eq.0)
     + latm(4,i)=1
      endif
   45 continue
c
c     sp2 carbon
c     ----------
c
      do 47 i=jlo,jhi
      if(latm(2,i).eq.6.and.latm(3,i).eq.3.and.latm(4,i).eq.0) then
      if(latm(4,i).eq.0) latm(4,i)=1
      endif
   47 continue
c
      argos_prepare_impctr=.true.
      return
      end
 
