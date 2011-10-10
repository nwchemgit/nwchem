      subroutine argos_space_wrtbxsz(lfn,boxsiz)
c
      implicit none
c
#include "argos_space_common.fh"
c
      real*8 boxsiz(maxbox,3)
      integer lfn
c
      integer i
c
      write(lfn,1000) 1,nbx,(boxsiz(i,1),i=1,nbx)
      write(lfn,1000) 2,nby,(boxsiz(i,2),i=1,nby)
      write(lfn,1000) 3,nbz,(boxsiz(i,3),i=1,nbz)
 1000 format(i3,i5,5f10.5,/,(8x,5f10.5))
c
      return
      end
c $Id$
