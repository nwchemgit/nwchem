      subroutine cputime(time)
      Implicit NONE
C     $Id$
      real*8 time
#include "tcgmsg.fh"
c -  convex
c      real*4 tx(2),tmptime
c      call etime(tx)
c      tmptime=tx(1)+tx(2)
c      time=dble(tmptime)
c - intel
c      real*8 dclock
c      time=dclock()
      time=tcgtime()
      return
      end
