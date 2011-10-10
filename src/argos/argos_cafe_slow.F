      subroutine argos_cafe_slow(xw,vw,numw,xs,vs,isat,nums)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer isat(msa)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
c
c      call argos_cafe_slowt(vw,numw,vs,nums,int_mb(i_iwa),isat,dbl_mb(i_mas))
      call argos_cafe_slowr(xw,vw,numw,xs,vs,nums,int_mb(i_iwa),isat,
     + dbl_mb(i_mas))
c
      return
      end
c $Id$
