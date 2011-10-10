      subroutine argos_cafe_fold(numw,xw,xwm,nums,isat,isml,xs,xsm)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3),xsm(msm,3)
      integer isat(msa),isml(msa)
c
      if(.not.lscr) then
      if(.not.ma_push_get(mt_dbl,3*mscr,'s3r1',l_s3r1,i_s3r1))
     + call md_abort('Failed to allocate scratch array s3r1',0)
      endif
c
      call argos_cafe_fld(numw,xw,xwm,nums,isat,isml,xs,xsm,
     + dbl_mb(i_s3r1))
c
      if(.not.lscr) then
      if(.not.ma_pop_stack(l_s3r1))
     + call md_abort('Failed to deallocate scratch array s3r1',0)
      endif
c
      return
      end
c $Id$
