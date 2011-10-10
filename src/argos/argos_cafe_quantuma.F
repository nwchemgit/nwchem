      integer function argos_cafe_quantuma(i)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer argos_cafe_atotyp
      external argos_cafe_atotyp
c
      integer i
c
      argos_cafe_quantuma=argos_cafe_atotyp(int_mb(i_typ),i)
c
      return
      end
c $Id$
