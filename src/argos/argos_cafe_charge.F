      real*8 function argos_cafe_charge(iq)
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
      argos_cafe_charge=argos_cafe_chg(iq,dbl_mb(i_chg))
      return
      end
c $Id$
