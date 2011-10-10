      subroutine argos_cafe_wrtmro(lfnmro)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfnmro
c
      call argos_cafe_wtmro(lfnmro,dbl_mb(i_sti))
c
      return
      end
c $Id$
