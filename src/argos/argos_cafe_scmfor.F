      subroutine argos_cafe_scmfor(icmopt,ismol,isatt,vs,fs,nums,fcm)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer icmopt,ismol(msa),isatt(msa),nums
      real*8 vs(msa,3),fs(msa,3),fcm(msm,5)
c
      call argos_cafe_scmf(icmopt,ismol,isatt,vs,fs,nums,fcm,
     + dbl_mb(i_mas))
c
      return
      end
c $Id$
