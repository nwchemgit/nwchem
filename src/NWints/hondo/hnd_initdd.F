c $Id: hnd_initdd.F,v 1.3 2002-01-16 18:25:28 windus Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : initialization routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_initdd(bases,nbas)
      implicit none
#include "hnd_tol.fh"
c
      integer nbas         ! [input] number of basis sets to be used
      integer bases(nbas)  ! [input] basis set handles
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
      call hnd_dd1emem(bases,nbas)
c
      return
      end
