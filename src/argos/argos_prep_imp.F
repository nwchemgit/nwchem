      logical function argos_prep_imp(latm,matm,natm,lbnd,mbnd,nbnd,
     + limp,mimp,nimp)
c
c $Id$
c
      implicit none
#include "argos_prep_common.fh"
c
      integer matm,natm,mbnd,nbnd,mimp,nimp
      integer latm(6,matm),lbnd(2,mbnd),limp(4,mimp)
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
      if(latm(2,i).ge.4.and.na.eq.3.and.
     + ((ffield(1:5).eq.'amber'.and.latm(2,i).le.6).or.
     + (ffield(1:5).eq.'charm'))) then
      nimp=nimp+1
      limp(1,nimp)=i
      limp(2,nimp)=ia(1)
      limp(3,nimp)=ia(2)
      limp(4,nimp)=ia(3)
      endif
c
      endif
    1 continue
c
      argos_prep_imp=.true.
c
      return
      end
