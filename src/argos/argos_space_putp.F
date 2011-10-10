      subroutine argos_space_putp(node,pw,pwp,numwm,ps,psp,numsa,lpp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl+12*mwa+3,jl+18*mwa+2,pw,mwm)
      if(lpp) call ga_put(ga_w,il,ih,jl+18*mwa+3,jl+30*mwa+2,pwp,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl+12,jl+17,ps,msa)
      if(lpp) call ga_put(ga_s,il,ih,jl+18,jl+30,psp,msa)
      endif
c
      return
      end
c $Id$
