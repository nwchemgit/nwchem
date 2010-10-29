      logical function argos_prep_dihimp(latm,matm,natm,lang,mang,nang,
     + ldih,mdih,ndih,limp,mimp,nimp)
c
c $Id$
c
      implicit none
c
#include "argos_prep_common.fh"
c
      integer matm,natm,mang,nang,mdih,ndih,mimp,nimp
      integer latm(6,matm),lang(3,mang),ldih(4,mdih),limp(4,mimp)
c
      integer i,j,k,no
c
      if(ffield(1:6).eq.'gromos'.or.ffield(1:6).eq.'charmm') then
c
      no=nimp
c
c     dihedrals in aromatic ring are treated as improper by gromos
c                                            both proper and improper by charmm
c
      do 1 i=1,ndih
      if(latm(2,ldih(1,i)).lt.4) goto 1
      if(latm(2,ldih(2,i)).lt.4) goto 1
      if(latm(2,ldih(3,i)).lt.4) goto 1
      if(latm(2,ldih(4,i)).lt.4) goto 1
      do 2 j=1,nang
      if(latm(2,lang(1,j)).lt.4) goto 2
      if(latm(2,lang(2,j)).lt.4) goto 2
      if(latm(2,lang(3,j)).lt.4) goto 2
      if((ldih(1,i).eq.lang(1,j).and.ldih(4,i).eq.lang(3,j)).or.
     + (ldih(1,i).eq.lang(3,j).and.ldih(4,i).eq.lang(2,j))) then
      do 3 k=no+1,nimp
      if(limp(1,k).eq.i) goto 2
    3 continue
      nimp=nimp+1
      limp(1,nimp)=i
      endif
    2 continue
    1 continue
c
      do 4 i=1,ndih-1
      if(latm(2,ldih(1,i)).lt.4) goto 4
      if(latm(2,ldih(2,i)).lt.4) goto 4
      if(latm(2,ldih(3,i)).lt.4) goto 4
      if(latm(2,ldih(4,i)).lt.4) goto 4
      do 5 j=i+1,ndih
      if(latm(2,ldih(1,j)).lt.4) goto 5
      if(latm(2,ldih(2,j)).lt.4) goto 5
      if(latm(2,ldih(3,j)).lt.4) goto 5
      if(latm(2,ldih(4,j)).lt.4) goto 5
      if((ldih(1,i).eq.ldih(1,j).and.ldih(4,i).eq.ldih(4,j)).or.
     + (ldih(1,i).eq.ldih(4,j).and.ldih(4,i).eq.ldih(4,j))) then
      do 6 k=no+1,nimp
      if(limp(1,k).eq.i) goto 7
    6 continue
      nimp=nimp+1
      limp(1,nimp)=i
    7 continue
      do 8 k=no+1,nimp
      if(limp(1,k).eq.j) goto 5
    8 continue
      nimp=nimp+1
      limp(1,nimp)=j
      endif
    5 continue
    4 continue
c
      do 9 i=no+1,nimp
      k=limp(1,i)
      do 10 j=1,4
      limp(j,i)=ldih(j,k)
   10 continue
      if(ffield(1:6).eq.'gromos') ldih(1,k)=0
    9 continue
c
      no=ndih
      ndih=0
      do 11 i=1,no
      if(ldih(1,i).gt.0) then
      ndih=ndih+1
      do 12 j=1,4
      ldih(j,ndih)=ldih(j,i)
   12 continue
      endif
   11 continue
c
      endif
c
      argos_prep_dihimp=.true.
c
      return
      end
