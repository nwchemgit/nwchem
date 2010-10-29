      subroutine argos_cafe_weight(numw,nums,isat,isml,ishop,wbx)
c $Id$
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
c
      integer numw,nums,isat(msa),isml(msa),ishop(msa)
      real*8 wbx
c
      call argos_cafe_wbox(numw,nums,int_mb(i_iwa),isat,isml,ishop,
     + dbl_mb(i_wsm),dbl_mb(i_mas))
c
      wbx=wbox
c
      return
      end
