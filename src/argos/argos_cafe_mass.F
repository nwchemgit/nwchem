      subroutine argos_cafe_mass(ww,ws,isat,nums)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nums,isat(msa)
      real*8 ww(mwa),ws(msa)
c
      call argos_cafe_m(ww,ws,isat,nums,int_mb(i_iwa),dbl_mb(i_mas))
c
      return
      end
c $Id$
