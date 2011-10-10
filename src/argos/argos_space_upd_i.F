      subroutine argos_space_upd_i(numsa,isl,islp,numwm,iwl,iwlp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numsa,isl(msa,mis2)
      integer numwm,iwl(mwm,miw2)
      integer iwlp(mwm,npackw),islp(msa,npack)
c
      integer il,ih,jl,jh
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,me,il,ih,jl,jh)
      if(npack.eq.0) then
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      else
      call argos_space_pack(numsa,isl,islp)
      call ga_put(ga_is,il,il+numsa-1,jl,jl+npack-1,islp,msa)
      endif
      endif
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,me,il,ih,jl,jh)
      if(npackw.eq.0) then
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      else
      call argos_space_packw(numwm,iwl,iwlp)
      call ga_put(ga_iw,il,il+numwm-1,jl,jl+npackw-1,iwlp,mwm)
      endif
      endif
c
      call ga_sync()
c
      return
      end
c $Id$
