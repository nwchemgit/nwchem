      subroutine argos_cafe_fnorm(iwdt,fw,numw,isat,isdt,fs,
     + nums,fnorm,fmax)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 fw(mwm,3,mwa),fs(msa,3),fnorm,fmax
      integer isat(msa),numw,nums
      integer isdt(msa),iwdt(mwm)
c
      call argos_cafe_fnrm(dbl_mb(i_mas),int_mb(i_iwa),iwdt,fw,numw,
     + isat,isdt,fs,nums,fnorm,fmax)
c
      return
      end
c $Id$
