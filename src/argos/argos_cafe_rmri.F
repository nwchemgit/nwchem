      logical function argos_cafe_rmri(lfnmri,dera)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer lfnmri
      real*8 dera(6,nsatot)
c
      integer i,j
c
      if(me.eq.0) then
      read(lfnmri,err=9,end=9) ((dera(i,j),i=1,6),j=1,nsatot)
      else
      do 1 i=1,nsatot
      do 2 j=1,6
      dera(j,i)=zero
    2 continue
    1 continue
      endif
c
      argos_cafe_rmri=.true.
      return
c
    9 continue
      argos_cafe_rmri=.false.
      return
      end
c $Id$
