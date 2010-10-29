c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : termination routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function hnd_terminate()
      implicit none
#include "apiP.fh"
c
c     at this point nothing to do but just in case for later
c
      hnd_terminate = .true.
      return
      end
