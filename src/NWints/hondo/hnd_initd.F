c $Id: hnd_initd.F,v 1.1 1999-04-10 00:15:16 windus Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : initialization routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_initd()
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
      call hnd_d1emem
c
      return
      end
