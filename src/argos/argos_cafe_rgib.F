      subroutine argos_cafe_rgib(lfngib,dera,dec,dlambd)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer lfngib
      real*8 dera(6,nsatot),dec(6,nsatot),dlambd
c
      integer i,j
c
      if(me.eq.0) then
      read(lfngib,1000) ((dera(i,j),i=1,6),j=1,nsatot)
 1000 format(4e20.10)
      do 1 i=1,nsatot
      do 2 j=1,6
      dec(j,i)=dec(j,i)+dlambd*dera(j,i)
    2 continue
    1 continue
      else
      do 3 i=1,nsatot
      do 4 j=1,6
      dera(j,i)=zero
    4 continue
    3 continue
      endif
c
      return
      end
c $Id$
