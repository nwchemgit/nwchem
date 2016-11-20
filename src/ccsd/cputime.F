      subroutine cputime(time)
      Implicit NONE
C     $Id$
      real*8 time
#include "tcgmsg.fh"
      time=tcgtime()
      return
      end
