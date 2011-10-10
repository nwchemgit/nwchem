      subroutine argos_cafe_lstmul(ixmul,imul,
     + msb,ibnd,msh,iang,msd,idih,mso,iimp)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer ixmul(mmuli),imul(mmult,4)
      integer msb,msh,msd,mso
      integer ibnd(msb,3),iang(msh,4),idih(msd,5),iimp(mso,5)
c
      integer i,j,k,m
c
      k=0
      m=0
c
      if(nmult(1).gt.0) then
      do 1 i=1,numb(2)
      if(ibnd(i,3).lt.0) then
      m=m+1
      ixmul(m)=i
      do 2 j=1,k
      if(ibnd(i,1).eq.imul(k,1)) goto 3
    2 continue
      k=k+1
      imul(k,1)=ibnd(i,1)
    3 continue
      do 4 j=1,k
      if(ibnd(i,2).eq.imul(k,1)) goto 1
    4 continue
      k=k+1
      imul(k,1)=ibnd(i,2)
c      if(me.eq.0) write(*,1000) ibnd(i,1),ibnd(i,2)
c 1000 format(' Multinode bond ',2i5)
      endif
    1 continue
      endif
c
      if(nmult(2).gt.0) then
      do 5 i=1,numh(2)
      if(iang(i,4).lt.0) then
      m=m+1
      ixmul(m)=i
      do 6 j=1,k
      if(iang(i,1).eq.imul(k,1)) goto 7
    6 continue
      k=k+1
      imul(k,1)=iang(i,1)
    7 continue
      do 8 j=1,k
      if(iang(i,2).eq.imul(k,1)) goto 9
    8 continue
      k=k+1
      imul(k,1)=iang(i,2)
    9 continue
      do 10 j=1,k
      if(iang(i,3).eq.imul(k,1)) goto 5
   10 continue
      k=k+1
      imul(k,1)=iang(i,3)
c      if(me.eq.0) write(*,1001) iang(i,1),iang(i,2),iang(i,3)
c 1001 format(' Multinode angle ',3i5)
      endif
    5 continue
      endif
c
      if(nmult(3).gt.0) then
      do 11 i=1,numd(2)
      if(idih(i,5).lt.0) then
      m=m+1
      ixmul(m)=i
      do 12 j=1,k
      if(idih(i,1).eq.imul(k,1)) goto 13
   12 continue
      k=k+1
      imul(k,1)=idih(i,1)
   13 continue
      do 14 j=1,k
      if(idih(i,2).eq.imul(k,1)) goto 15
   14 continue
      k=k+1
      imul(k,1)=idih(i,2)
   15 continue
      do 16 j=1,k
      if(idih(i,3).eq.imul(k,1)) goto 17
   16 continue
      k=k+1
      imul(k,1)=idih(i,3)
   17 continue
      do 18 j=1,k
      if(idih(i,4).eq.imul(k,1)) goto 11
   18 continue
      k=k+1
      imul(k,1)=idih(i,4)
c      if(me.eq.0) write(*,1002) idih(i,1),idih(i,2),idih(i,3),idih(i,4)
c 1002 format(' Multinode torsion ',4i5)
      endif
   11 continue
      endif
c
      if(nmult(4).gt.0) then
      do 19 i=1,numi(2)
      if(iimp(i,5).lt.0) then
      m=m+1
      ixmul(m)=i
      do 20 j=1,k
      if(iimp(i,1).eq.imul(k,1)) goto 21
   20 continue
      k=k+1
      imul(k,1)=iimp(i,1)
   21 continue
      do 22 j=1,k
      if(iimp(i,2).eq.imul(k,1)) goto 23
   22 continue
      k=k+1
      imul(k,1)=iimp(i,2)
   23 continue
      do 24 j=1,k
      if(iimp(i,3).eq.imul(k,1)) goto 25
   24 continue
      k=k+1
      imul(k,1)=iimp(i,3)
   25 continue
      do 26 j=1,k
      if(iimp(i,4).eq.imul(k,1)) goto 19
   26 continue
      k=k+1
      imul(k,1)=iimp(i,4)
c      if(me.eq.0) write(*,1003) iimp(i,1),iimp(i,2),iimp(i,3),iimp(i,4)
c 1003 format(' Multinode improper ',4i5)
      endif
   19 continue
      endif
c
      nmul=k
c
      return
      end
c $Id$
