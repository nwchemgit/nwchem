      subroutine argos_cafe_trjhdr(iunit)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer iunit
c
      call argos_cafe_wthdr(iunit,byte_mb(i_snam),
     + int_mb(i_ibnd(1)),mbt(1),int_mb(i_ibnd(2)),mbt(2))
c
      return
      end
c $Id$
