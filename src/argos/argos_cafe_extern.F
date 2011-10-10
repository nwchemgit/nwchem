      subroutine argos_cafe_extern(stimei,nums,fs,isq,numw,fw)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nums,numw
      integer isq(msa)
      real*8 stimei,fs(msa,3),fw(mwm,3,mwa)
c
      stime=stimei
c
      call argos_cafe_external_field(nums,fs,isq,numw,fw,int_mb(i_iwq),
     + dbl_mb(i_chg))
c
      return
      end
c $Id$
