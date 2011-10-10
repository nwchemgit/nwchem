      real*8 function argos_cafe_atmass(i,isat)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 argos_cafe_am
      external argos_cafe_am
c
      integer i,isat(msa)
c
      argos_cafe_atmass=argos_cafe_am(isat(i),dbl_mb(i_mas))
c
      return
      end
c $Id$
