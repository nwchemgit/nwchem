c $Id: hnd_init.F,v 1.2 1999-01-22 19:01:17 windus Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : initialization routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_init()
      implicit none
#include "hnd_tol.fh"
c
c     set default tolerance
c
      itol = 15
c
c     set up the HONDO hermit and wermit arrays
c
      call hnd_whermt
c
c     set up memory requirements
c
      call hnd_1emem
c
      return
      end
