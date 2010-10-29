      subroutine int_init_1eelec(maxbuf,maxscr,basis,nder,npt)
c
c $Id$
c
      implicit none
c
#include "hnd_pointers.fh"
#include "bas.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer maxbuf    ! [output] buffer space for integrals
      integer maxscr    ! [output] scratch space needed
      integer basis     ! [input] basis set info
      integer nder      ! [input] integral derivative to be calculated
      integer npt       ! [input] number of grid points to be calculated
c
      integer lmax, iblk, i, mroots
c
c     Define Hermite polynomial values
c
      call hnd_whermt
c
c     set up memory requirements: get max L first
c
      lmax = -1
      if (.not.bas_high_angular(basis,lmax))
     &   call errquit('int_init_1eelec: lmax error',lmax, INT_ERR)
      lmax = lmax + 1
c
c     Calculate maximum buffer needed
c
      maxbuf = lmax*(lmax+1)*lmax*(lmax+1)/4 * npt * max(nder*3,1)
c
c     Calculate maximum scratch space needed
c
      mroots = (2*lmax+max(nder,0)-2)/2+1
      maxscr = 3*lmax*lmax*mroots
      if (nder.ge.1) maxscr = (nder+1)*maxscr
c
c     Set pointers for scratch space
c
      iblk = lmax*lmax*mroots
      do i = 1, 9
         elpt(i) = 1
      enddo
      do i = 2, 2+max(nder*3+1,1)
         elpt(i) = 1 + (i-1)*iblk
      enddo
c
c     If spherical basis, allocate scratch space for transformation
c     A little smaller than maxbuf: we will use maxbuf size
c
      maxscr = max(maxscr,maxbuf)
c
      return
      end
