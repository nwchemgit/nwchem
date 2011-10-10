      subroutine argos_cafe_mdstep(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwgmn,isgan,issgm,tmpx,ishp)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa),isgan(msa),issgm(msa)
      integer iwgmn(mwm),ishp(msa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
      real*8 tmpx
c
      tmpext=tmpx
c
      call argos_cafe_tscale()
c
      if(mdalgo.eq.1) then
      call argos_cafe_lpfrog(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,int_mb(i_iwa),
     + dbl_mb(i_mas),iwgmn,ishp)
      elseif(mdalgo.eq.2) then
      call argos_cafe_lpfrog_bc(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,int_mb(i_iwa),
     + dbl_mb(i_mas),iwgmn,isgan)
      else
      call md_abort('Unknown integration algorithm',mdalgo)
      endif
c
      return
      end
c $Id$
