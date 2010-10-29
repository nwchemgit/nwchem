      subroutine int_init_dip(maxbuf,maxscr,basis)
c
c $Id$
c
      implicit none
c
#include "hnd_pointers.fh"
#include "bas.fh"
#include "errquit.fh"
#include "hnd_tol.fh"
c
      integer maxbuf    ! [output] buffer space for integrals
      integer maxscr    ! [output] scratch space needed
      integer basis     ! [input] basis set info
c
      integer lmax, i
c
c     set default tolerance
c
      itol = 15
c
c     Define Hermite polynomial values
c
      call hnd_whermt
c
c     set up memory requirements: get max L first
c
      lmax = -1
      if (.not.bas_high_angular(basis,lmax))
     &   call errquit('int_init_dso: lmax error',lmax, INT_ERR)
      lmax = lmax + 1
c
c     Calculate maximum buffer needed
c
      maxbuf = lmax*(lmax+1)*lmax*(lmax+1)/4
c
c     Calculate maximum scratch space needed
c
      maxscr = 9*lmax*(lmax+1)
      maxscr = max(maxbuf,maxscr)
c
c     Set pointers for scratch space
c
      dip(1) = 1
      do i = 2, 9
         dip(i) = 1 + (i-1)*lmax*(lmax+1)
      enddo
c
c     If spherical basis, allocate scratch space for transformation
c     A little smaller than maxbuf: we will use maxbuf size
c
      maxscr = max(maxscr,maxbuf)
c
      return
      end
