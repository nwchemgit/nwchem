      subroutine argos_cafe_vscale(tgauss,numw,nums,vw,vwt,vs,vst,iwdt,
     + isdt,isat)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 tgauss
      integer numw,nums
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),vs(msa,3),vst(msa,3)
      integer iwdt(mwm),isdt(msa),isat(msa)
c
      call argos_cafe_vscal(tgauss,
     + numw,nums,vw,vwt,vs,vst,iwdt,isdt,int_mb(i_iwa),
     + isat,dbl_mb(i_mas))
c
      return
      end
c $Id$
