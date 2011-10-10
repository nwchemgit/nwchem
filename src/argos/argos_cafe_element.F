      integer function argos_cafe_element(i)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer argos_cafe_atonum
      external argos_cafe_atonum
c
      integer i
c
      argos_cafe_element=argos_cafe_atonum(i,int_mb(i_num))
c
      return
      end
c $Id$
