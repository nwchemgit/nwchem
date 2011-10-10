      subroutine argos_cafe_dera_init()
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      call argos_cafe_dera_i(dbl_mb(i_sti))
c
      return
      end
c $Id$
