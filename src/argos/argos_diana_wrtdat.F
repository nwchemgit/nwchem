      subroutine argos_diana_wrtdat(tval,ntdim,nfr)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "argos_diana_common.fh"
#include "mafdecls.fh"
c
      integer ntdim,nfr
      real*8 tval(ntdim,nfrdim)
c
      integer i,j
c
      do 1 j=1,nfr
      write(lfndat,1000)
 1000 format('frame')
      write(lfndat,1001) (tval(i,j),i=1,ntdim)
 1001 format(6e12.5)
    1 continue
c
      return
      end
c $Id$
