c $Id: hnd_initdd.F,v 1.2 1999-05-06 17:25:50 windus Exp $
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
      write(6,*) 'before the call to hnd_dd1emem'
      call flush(6)
      call hnd_dd1emem(bases,nbas)
      write(6,*) 'after the call to hnd_dd1emem'
      call flush(6)
c
      return
      end
