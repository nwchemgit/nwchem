      subroutine argos_diana_xaver(xtrj,xa)
c
c      x,xa,nx,mx,nave)
c
c $Id$
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "argos_diana_common.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      real*8 xtrj(ntrj,3,nfrdim),xa(nsa,3)
c
      integer i,j,k
c
      if(nsa.ne.ntrj) call md_abort('xaver dimension problem',0)
c
      do 1 j=1,3
      do 2 i=1,nsa
      xa(i,j)=zero
    2 continue
      do 3 k=1,nfrme
      do 4 i=1,ntrj
      xa(i,j)=xa(i,j)+xtrj(i,j,k)
    4 continue
    3 continue
    1 continue
c
c     global sum
c
      call ga_dgop(mag_d08,xa,3*nsa,'+')
c
      do 5 j=1,3
      do 6 i=1,ntrj
      xa(i,j)=xa(i,j)/dble(nfrtot)
    6 continue
    5 continue
c
      if(me.eq.0) then
      write(*,'(a)') ' Coordinates averaged'
      endif
c
      return
      end
