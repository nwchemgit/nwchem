      subroutine argos_cafe_coms(isat,isml,xs,xsm,gsm,isfr,isto)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer isfr,isto
      integer isat(msa),isml(msa)
      real*8 xs(msa,3),xsm(msm,3),gsm(msm,4,2)
c
      call argos_cafe_cms(isat,isml,xs,xsm,gsm,isfr,isto,dbl_mb(i_wsm),
     + dbl_mb(i_mas))
c
      return
      end
c $Id$
