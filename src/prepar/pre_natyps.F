      logical function pre_natyps(catm,matm,natm,natt)
c
c $Id: pre_natyps.F,v 1.2 1998-04-29 18:35:19 d3j191 Exp $
c
      implicit none
c
      integer matm,natm,natt
      character*6 catm(4,matm)
c
      integer i,j
c
      natt=1
c
      do 1 i=2,natm
      do 2 j=1,i-1
      if(catm(2,i).eq.catm(2,j).and.catm(3,i).eq.catm(3,j).and.
     + catm(4,i).eq.catm(4,j)) goto 1
    2 continue
      natt=natt+1
    1 continue
c
      pre_natyps=.true.
      return
      end
