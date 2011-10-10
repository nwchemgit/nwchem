      subroutine argos_space_gagetixv(node,iwl,iwlp,xw,xwcr,
     + vw,numwm,isl,islp,
     + xs,vs,numsa,ipl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,ipl,mbox)
      numwm=ipl(1,2)
      numsa=ipl(2,2)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_get(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call ga_get(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      call argos_space_unpackw(numwm,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_get(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      call ga_get(ga_w,il,ih,jl+3*mwa,jl+6*mwa-1,vw,mwm)
      call ga_get(ga_w,il,ih,jl+6*mwa,jl+6*mwa+2,xwcr,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_get(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call ga_get(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      call argos_space_unpack(numsa,isl,islp)
      endif
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_get(ga_s,il,ih,jl,jl+2,xs,msa)
      call ga_get(ga_s,il,ih,jl+3,jl+5,vs,msa)
      endif
c
      return
      end
c $Id$
