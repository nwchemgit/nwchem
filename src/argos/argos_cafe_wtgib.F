      subroutine argos_cafe_wtgib(lfngib,dera)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
c
      integer lfngib
      real*8 dera(6,nsatot)
c
      integer i,j
c
      call ga_dgop(mcf_68,dera,5*nsatot,'+')
c
      if(me.eq.0) then
      write(lfngib,1000) ((dera(i,j),i=1,6),j=1,nsatot)
 1000 format(4e20.10)
      else
      do 1 i=1,nsatot
      do 2 j=1,6
      dera(j,i)=zero
    2 continue
    1 continue
      endif
c
      return
      end
c $Id$
