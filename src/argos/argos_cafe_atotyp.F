      integer function argos_cafe_atotyp(ityp,ia)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer ia
      integer ityp(mat,3)
c
      argos_cafe_atotyp=ityp(ia,iset)
c
      return
      end
c $Id$
