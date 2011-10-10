      subroutine argos_space_bscale(boxsiz)
c
      implicit none
c
#include "argos_space_common.fh"
c
      real*8 boxsiz(maxbox,3)
c
      integer ibx,iby,ibz
      real*8 boxscl
c
c     scale the subbox boundaries
c
      boxscl=box(1)/boxsiz(nbx,1)
      do 1 ibx=1,nbx
      boxsiz(ibx,1)=boxscl*boxsiz(ibx,1)
    1 continue
      boxscl=box(2)/boxsiz(nby,2)
      do 2 iby=1,nby
      boxsiz(iby,2)=boxscl*boxsiz(iby,2)
    2 continue
      boxscl=box(3)/boxsiz(nbz,3)
      do 3 ibz=1,nbz
      boxsiz(ibz,3)=boxscl*boxsiz(ibz,3)
    3 continue
c
      return
      end
c $Id$
