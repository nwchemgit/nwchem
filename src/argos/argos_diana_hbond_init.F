      subroutine argos_diana_hbond_init(snam,ibnd,ihbnd)
c
c $Id: argos_diana_hbond.F
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      character*16 snam(msa)
      integer ibnd(msb,2),ihbnd(msa)
c
      integer i,j,k
c
      if(np.gt.1)
     + call md_abort('h-bond analysis on single processor only',0)
c
      do 1 i=1,nsa
      ihbnd(i)=0
    1 continue
c
      do 2 j=1,2
      do 3 i=1,msb
      if(snam(ibnd(i,j))(6:7).eq.' O'.or.
     + snam(ibnd(i,j))(6:7).eq.' N'.or.
     + snam(ibnd(i,j))(6:7).eq.' S') ihbnd(ibnd(i,j))=-1
    3 continue
    2 continue
c
      do 4 j=1,2
      k=2
      if(j.eq.2) k=1
      do 5 i=1,msb
      if(ihbnd(ibnd(i,k)).eq.-1) then
      if(snam(ibnd(i,j))(6:7).eq.' H'.or.
     + snam(ibnd(i,j))(6:7).eq.'2H'.or.
     + snam(ibnd(i,j))(6:7).eq.'3H'.or.
     + snam(ibnd(i,j))(6:7).eq.'4H'.or.
     + snam(ibnd(i,j))(6:7).eq.'5H'.or.
     + snam(ibnd(i,j))(6:7).eq.'6H'.or.
     + snam(ibnd(i,j))(6:7).eq.'7H'.or.
     + snam(ibnd(i,j))(6:7).eq.'8H'.or.
     + snam(ibnd(i,j))(6:7).eq.'9H') ihbnd(ibnd(i,j))=ibnd(i,k)
      endif
    5 continue
    4 continue
c
      do 6 i=1,nwa
      ihbndw(i)=0
    6 continue
c
      do 7 j=1,2
      do 8 i=1,mwb
      if(wnam(ibndw(i,j))(6:7).eq.' O'.or.
     + wnam(ibndw(i,j))(6:7).eq.' N'.or.
     + wnam(ibndw(i,j))(6:7).eq.' S') ihbndw(ibnd(i,j))=-1
    8 continue
    7 continue
c
      do 9 j=1,2
      k=2
      if(j.eq.2) k=1
      do 10 i=1,mwb
      if(ihbndw(ibndw(i,k)).eq.-1) then
      if(wnam(ibndw(i,j))(6:7).eq.' H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'2H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'3H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'4H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'5H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'6H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'7H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'8H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'9H') ihbndw(ibndw(i,j))=ibndw(i,k)
      endif
   10 continue
    9 continue
c
      return
      end
