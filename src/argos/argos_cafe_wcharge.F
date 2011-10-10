      real*8 function argos_cafe_wcharge(iq)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 argos_cafe_chg
      external argos_cafe_chg
c
      integer iq
c
      argos_cafe_wcharge=argos_cafe_chg(int_mb(i_iwq+iq-1),
     + dbl_mb(i_chg))
      return
      end
c $Id$
