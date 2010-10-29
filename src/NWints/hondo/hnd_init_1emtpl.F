      subroutine int_init_1emtpl(maxbuf,maxscr,basis,mpole)
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
      integer mpole     ! [input] order of multipole
c
      integer lmax, iblk, i, mroots, lsz_buf, memsize
c
c     Define Hermite polynomial values
c
      call hnd_whermt
c
c     set up memory requirements: get max L first
c
      lmax = -1
      if (.not.bas_high_angular(basis,lmax))
     &   call errquit('int_mem_1emtpl: lmax error',lmax, INT_ERR)
      lmax = lmax + 1
c
c     Calculate maximum buffer needed
c
      maxbuf = lmax*(lmax+1)*lmax*(lmax+1)
      maxbuf = maxbuf * ((mpole+1)*(mpole+2)*(mpole+3))/6
c
c     Calculate maximum scratch space needed
c
      call emem_3ov(basis,basis,basis,lsz_buf,memsize)
      maxbuf = max(maxbuf,lsz_buf)
      maxscr = memsize
c
      return
      end
