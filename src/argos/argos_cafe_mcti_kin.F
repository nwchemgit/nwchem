      subroutine argos_cafe_mcti_kin(isat,isga,vs,nums)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nums
      integer isat(msa),isga(msa)
      real*8 vs(msa,3)
c
      call argos_cafe_mcti_w(isat,isga,vs,nums,dbl_mb(i_mas),
     + dbl_mb(i_sti))
c
      return
      end
c $Id$
