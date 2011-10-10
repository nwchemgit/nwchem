      subroutine argos_timer_init()
c
      implicit none
c
#include "argos.fh"
c
      integer i
c
      do 1 i=1,200
      call timer_reset(i)
    1 continue
c
      return
      end
c $Id$
