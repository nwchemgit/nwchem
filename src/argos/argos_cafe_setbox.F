      subroutine argos_cafe_setbox(vl)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 vl(3,3)
c
      integer i,j
c
      do 1 i=1,3
      box(i)=vl(i,i)
      boxh(i)=half*box(i)
      do 2 j=1,3
      vlat(i,j)=vl(i,j)
    2 continue
    1 continue
      volume=box(1)*box(2)*box(3)
c
      return
      end
c $Id$
