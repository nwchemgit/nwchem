      subroutine argos_space_gaputix(node,iwl,iwlp,xw,numwm,
     + isl,islp,xs,numsa,ipl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),iwlp(mwm,npackw)
      integer isl(msa,mis2),islp(msa,npack)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      ipl(1,2)=numwm
      ipl(2,2)=numsa
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call argos_space_packw(numwm,iwl,iwlp)
      call ga_put(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      endif
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call argos_space_pack(numsa,isl,islp)
      call ga_put(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      endif
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl,jl+2,xs,msa)
      endif
c
      return
      end
c $Id$
