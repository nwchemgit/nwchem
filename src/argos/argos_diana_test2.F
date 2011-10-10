      subroutine argos_diana_test2(xtest,nxdim,nfr)
#include "argos_diana_common.fh"
      integer nxdim
      real*8 xtest(nxdim,3,nfr)
c
      integer i,j,k
c
      do 2 k=1,nfr
      write(*,'(a,i5)') 'FRAME ',me*nfr+k
      do 1 i=1,10
      write(*,'(i5,3f12.6)') i,(xtest(i,j,k),j=1,3)
    1 continue
    2 continue
c
      return
      end
c $Id$
