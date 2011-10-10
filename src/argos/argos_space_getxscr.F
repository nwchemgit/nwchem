      subroutine argos_space_getxscr(xscr,xscrs)
c
      implicit none
c
#include "argos_space_common.fh"
c
      real*8 xscr(msm,3),xscrs(msm,3)
c
      integer i,j
c
      do 1 j=1,3
      do 2 i=1,nsm
      xscr(i,j)=xscrs(i,j)
    2 continue
    1 continue
c
      return
      end
c $Id$
