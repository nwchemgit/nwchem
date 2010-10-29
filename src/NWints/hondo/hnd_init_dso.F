      subroutine int_init_dso(maxbuf,maxscr,basis,npt)
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
      integer npt       ! [input] number of grid points to be calculated
c
      integer lmax, iblk, i, mroots
c
c     set default tolerance
c
      itol = 15
c
c     Define Hermite polynomial values
c
      call hnd_whermt
      call hnd_wleg
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
      maxbuf = lmax*(lmax+1)*lmax*(lmax+1)/4*npt*9
c
c     Calculate maximum scratch space needed
c
      mroots = (2*lmax+2-2)/2+1
      iblk = lmax*lmax*4*mroots
      maxscr = 3*iblk
      maxscr = max(maxbuf,maxscr)
c
c     Set pointers for scratch space
c
      dsot(1) = 1
      do i = 2, 3
         dsot(i) = 1 + (i-1)*iblk
      enddo
c
c     If spherical basis, allocate scratch space for transformation
c     A little smaller than maxbuf: we will use maxbuf size
c
      maxscr = max(maxscr,maxbuf)
c
      return
      end
