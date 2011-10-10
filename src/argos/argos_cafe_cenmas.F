      subroutine argos_cafe_cenmas(numw,xw,xwm,nums,
     + isat,isml,xs,xsm,gsm)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),xwm(mwm,3)
      real*8 xs(msa,3),xsm(msm,3),gsm(msm,4,2)
      integer isat(msa),isml(msa)
c
      call argos_cafe_comw(xw,xwm,1,numw)
      call argos_cafe_coms(isat,isml,xs,xsm,gsm,1,nums)
c
      return
      end
c $Id$
