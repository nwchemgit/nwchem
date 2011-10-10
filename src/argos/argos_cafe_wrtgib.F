      subroutine argos_cafe_wrtgib(lfngib)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfngib
c
      call argos_cafe_wtgib(lfngib,dbl_mb(i_sti))
c
      return
      end
c $Id$
