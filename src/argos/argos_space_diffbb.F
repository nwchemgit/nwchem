      logical function argos_space_diffbb(boxsiz,ibxrng)
c
      implicit none
c
#include "argos_space_common.fh"
c
      integer ibxrng(maxbox,3,2)
      real*8 boxsiz(maxbox,3)
      integer i,ib,j,k,num
      real*8 rd
c
c     temporarily skip for non-rectangular boxes
c
      if(nbxtyp.eq.1) then
      argos_space_diffbb=.true.
      return
      endif
c
c     initialize
c
      if(nbbdif.lt.0) then
      do 1 k=1,2
      do 2 j=1,3
      do 3 i=1,maxbox
      ibxrng(i,j,k)=-1
    3 continue
    2 continue
    1 continue
      endif
c
      nbbdif=0
c
      do 4 i=1,nbx
      num=0
      rd=zero
      do 5 j=1,nbx
      ib=i+1
      if(ib.gt.nbx) ib=ib-nbx
      rd=rd+boxsiz(ib,1)
      if(rd.gt.rlong) goto 6
      num=num+1
    5 continue
    6 continue
      if(num.ne.ibxrng(i,1,1)) then
      ibxrng(i,1,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 7 j=1,nbx
      ib=i-1
      if(ib.lt.1) ib=ib+nbx
      rd=rd+boxsiz(ib,1)
      if(rd.gt.rlong) goto 8
      num=num+1
    7 continue
    8 continue
      if(num.ne.ibxrng(i,1,2)) then
      ibxrng(i,1,2)=num
      nbbdif=nbbdif+1
      endif
    4 continue
c
      do 9 i=1,nby
      num=0
      rd=zero
      do 10 j=1,nby
      ib=i+1
      if(ib.gt.nby) ib=ib-nby
      rd=rd+boxsiz(ib,2)
      if(rd.gt.rlong) goto 11
      num=num+1
   10 continue
   11 continue
      if(num.ne.ibxrng(i,2,1)) then
      ibxrng(i,2,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 12 j=1,nby
      ib=i-1
      if(ib.lt.1) ib=ib+nby
      rd=rd+boxsiz(ib,2)
      if(rd.gt.rlong) goto 13
      num=num+1
   12 continue
   13 continue
      if(num.ne.ibxrng(i,2,2)) then
      ibxrng(i,2,2)=num
      nbbdif=nbbdif+1
      endif
    9 continue
c
      do 14 i=1,nbz
      num=0
      rd=zero
      do 15 j=1,nbz
      ib=i+1
      if(ib.gt.nbz) ib=ib-nbz
      rd=rd+boxsiz(ib,3)
      if(rd.gt.rlong) goto 16
      num=num+1
   15 continue
   16 continue
      if(num.ne.ibxrng(i,3,1)) then
      ibxrng(i,3,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 17 j=1,nbz
      ib=i-1
      if(ib.lt.1) ib=ib+nbz
      rd=rd+boxsiz(ib,3)
      if(rd.gt.rlong) goto 18
      num=num+1
   17 continue
   18 continue
      if(num.ne.ibxrng(i,3,2)) then
      ibxrng(i,3,2)=num
      nbbdif=nbbdif+1
      endif
   14 continue
c
      argos_space_diffbb=nbbdif.ne.0
c
      return
      end
c $Id$
