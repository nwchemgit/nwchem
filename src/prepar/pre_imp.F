      logical function pre_imp(ffield,latm,matm,natm,lbnd,mbnd,nbnd,
     + limp,mimp,nimp)
c
      implicit none
c
      integer matm,natm,mbnd,nbnd,mimp,nimp
      integer latm(5,matm),lbnd(2,mbnd),limp(4,mimp)
      character*80 ffield
c
      integer i,j,k,l,ltemp,na
      integer ia(10),nia(10)
c
      nimp=0
c
c     loop over all atoms
c
      do 1 i=1,natm
      if(latm(2,i).gt.0) then
c
c     count the number of bonds to this atom
c
      na=0
      do 2 j=1,nbnd
      if(lbnd(1,j).eq.i) then
      na=na+1
      ia(na)=lbnd(2,j)
      endif
      if(lbnd(2,j).eq.i) then
      na=na+1
      ia(na)=lbnd(1,j)
      endif
    2 continue
c
      nia(1)=0
      nia(2)=0
      nia(3)=0
      do 54 k=1,nbnd
      if(lbnd(1,k).eq.ia(1)) nia(1)=nia(1)+1
      if(lbnd(2,k).eq.ia(1)) nia(1)=nia(1)+1
      if(lbnd(1,k).eq.ia(2)) nia(2)=nia(2)+1
      if(lbnd(2,k).eq.ia(2)) nia(2)=nia(2)+1
      if(lbnd(1,k).eq.ia(3)) nia(3)=nia(3)+1
      if(lbnd(2,k).eq.ia(3)) nia(3)=nia(3)+1
   54 continue
      do 55 k=1,2
      do 56 l=2,3
      if(nia(k).lt.nia(l)) then
      ltemp=nia(k)
      nia(k)=nia(l)
      nia(l)=ltemp
      ltemp=ia(k)
      ia(k)=ia(l)
      ia(l)=ltemp
      endif
   56 continue
   55 continue
c
c     planar
c
      if(latm(2,i).eq.1.and.na.eq.3) then
      nimp=nimp+1
      limp(1,nimp)=i
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
c      write(*,'(6i5)') nimp,1,(limp(j,nimp),j=1,4)
      endif
c
c     chiral
c
      if((latm(2,i).eq.2.or.latm(2,i).eq.3).and.na.eq.3) then
      nimp=nimp+1
      limp(1,nimp)=i
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      endif
c
c     aromatic ring
c
      if(ffield(1:5).eq.'amber'.or.ffield(1:5).eq.'charm') then
      if(latm(2,i).ge.4.and.na.eq.3) then
      nimp=nimp+1
      limp(1,nimp)=i
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
c      write(*,'(6i5)') nimp,4,(limp(j,nimp),j=1,4)
      endif
      endif
c
      endif
    1 continue
c
      pre_imp=.true.
c
      return
      end
