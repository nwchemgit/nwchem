      logical function pre_natyps(catm,matm,natm,natt)
c
c $Id: pre_natyps.F,v 1.3 1998-05-05 22:09:22 d3j191 Exp $
c
      implicit none
c
      integer matm,natm,natt
      character*6 catm(4,matm)
c
      integer i,j,ii,jj,iii
c
      natt=0
c
      do 1 i=1,natm
      do 2 ii=1,3
      do 3 j=1,i-1
      do 4 jj=1,3
      if(catm(ii+1,i).eq.catm(jj+1,j)) goto 2
    4 continue
    3 continue
      do 5 iii=1,ii-1
      if(catm(ii+1,i).eq.catm(iii+1,i)) goto 2
    5 continue
      natt=natt+1
      write(*,'(i4,3a6)') natt,catm(ii+1,i)
    2 continue
    1 continue
c
      pre_natyps=.true.
      return
      end
