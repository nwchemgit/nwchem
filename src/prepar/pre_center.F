      logical function pre_center(xatm,latm,matm,jlo,ilo,ihi,jhi,
     + lbnd,mbnd,nbnd,lang,mang,nang,ldih,mdih,ndih,lring,aring,
     + mring,nring5,nring6)
c
c     function to determine the center type
c
c     in  : xatm(1:3,matm)  = atomic coordinates
c           latm(2,matm)    = atomic numbers
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
      logical pre_torang
      external pre_torang
c
      integer matm,ilo,ihi,jlo,jhi,mbnd,nbnd,mang,nang,mdih,ndih
      real*8 xatm(3,matm)
      integer latm(5,matm),lbnd(2,mbnd),lang(3,mang),ldih(4,mdih)
c
      integer mring
      integer lring(6,mring)
      logical aring(mring)
      integer nring5,nring6
      logical valid
c
      integer i,j,k,l,li,lj
      real*8 angle
c
      pre_center=.false.
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
      latm(4,i)=1
      latm(4,j)=1
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
c     5 membered rings
c     ----------------
c
      nring5=0
      do 4 li=1,ndih
      do 5 lj=1,nang
      if((ldih(1,li).eq.lang(1,lj).and.ldih(4,li).eq.lang(3,lj)).or.
     + (ldih(1,li).eq.lang(3,lj).and.ldih(4,li).eq.lang(1,lj))) then
      nring5=nring5+1
      if(nring5.gt.mring)
     +  call errquit('increase mring in pre_center',9999)
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
      do 9 i=1,nring5-1
      do 10 j=1,5
      if(lring(j,i).ne.lring(j,nring5)) goto 9
   10 continue
      nring5=nring5-1
      goto 4
    9 continue
      endif
    5 continue
    4 continue
c
c     6 membered rings
c     ----------------
c
      nring6=nring5
      do 11 li=1,ndih-1
      do 12 lj=li+1,ndih
      if((ldih(1,li).eq.ldih(1,lj).and.ldih(4,li).eq.ldih(4,lj)).or.
     + (ldih(1,li).eq.ldih(4,lj).and.ldih(4,li).eq.ldih(1,lj))) then
      nring6=nring6+1
      if(nring6.gt.mring)
     +  call errquit('increase mring in pre_center',9999)
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
      endif
   12 continue
   11 continue
c
c     select aromatic 5-membered rings
c     --------------------------------
c
      do 18 i=1,nring5
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
      pre_center=.true.
      return
      end
      
