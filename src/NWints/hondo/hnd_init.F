c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : initialization routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_init(bases,nbas)
      implicit none
#include "hnd_tol.fh"
#include "bas.fh"
#include "errquit.fh"
#include "hnd_whermt.fh"
c
      integer nbas         ! [input] number of basis sets to be used
      integer bases(nbas)  ! [input] basis set handles
c
      integer ibas,maxang,angm
c
c     set default tolerance
c
      itol = 15
c
c     set up the HONDO hermit and wermit arrays
c
      hermit_done = .false.
      call hnd_whermt
c
c     set up memory requirements: get max L first
c
      maxang = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('hnd_init: angm error',angm, INT_ERR)
        maxang = max(maxang,angm)
00300 continue
c
      call hnd_1emem(maxang)
c
      return
      end
