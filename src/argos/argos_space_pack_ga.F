      subroutine argos_space_pack_ga(numsa,isl,islp,numwm,iwl,iwlp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numsa,isl(msa,mis2),islp(msa,npack)
      integer numwm,iwl(mwm,miw2),iwlp(mwm,npackw)
c
      integer il,ih,jl,jh
c
      if(numsa.gt.0.and.npack.gt.0) then
      call ga_distribution(ga_is,me,il,ih,jl,jh)
      call ga_get(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      call argos_space_pack(numsa,isl,islp)
      call ga_put(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      endif
c
      if(numwm.gt.0.and.npackw.gt.0) then
      call ga_distribution(ga_iw,me,il,ih,jl,jh)
      call ga_get(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      call argos_space_packw(numwm,iwl,iwlp)
      call ga_put(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      endif
c
      call ga_sync()
c
      return
      end
c $Id$
