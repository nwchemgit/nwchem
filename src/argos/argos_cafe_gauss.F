      subroutine argos_cafe_gauss(tgauss,frgaus,numw,nums,vw,vs,
     + iwdt,isdt,isat)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 tgauss,frgaus
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer iwdt(mwm),isdt(msa),isat(msa)
c
      call argos_cafe_gauss1(tgauss,frgaus,
     + numw,nums,vw,vs,iwdt,isdt,int_mb(i_iwa),
     + isat,dbl_mb(i_mas))
c
      return
      end
c $Id$
