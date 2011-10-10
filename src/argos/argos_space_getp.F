      subroutine argos_space_getp(node,pw,pwp,numwm,ps,psp,
     + numsa,lpp,ndx)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa,ndx
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpp
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      if(ndx.eq.1) then
      call ga_get(ga_w,il,ih,jl+12*mwa+3,jl+15*mwa+2,pw,mwm)
      if(lpp) call ga_get(ga_w,il,ih,jl+18*mwa+3,jl+24*mwa+2,pwp,mwm)
      elseif(ndx.eq.2) then
      call ga_get(ga_w,il,ih,jl+15*mwa+3,jl+18*mwa+2,pw(1,1,1,2),mwm)
      if(lpp) call ga_get(ga_w,il,ih,jl+18*mwa+3,jl+24*mwa+2,pwp,mwm)
      else
      call ga_get(ga_w,il,ih,jl+12*mwa+3,jl+18*mwa+2,pw,mwm)
      if(lpp) call ga_get(ga_w,il,ih,jl+18*mwa+3,jl+24*mwa+2,pwp,mwm)
      endif
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      if(ndx.eq.1) then
      call ga_get(ga_s,il,ih,jl+12,jl+14,ps,msa)
      if(lpp) call ga_get(ga_s,il,ih,jl+15,jl+20,psp,msa)
      elseif(ndx.eq.2) then
      call ga_get(ga_s,il,ih,jl+15,jl+17,ps(1,1,2),msa)
      if(lpp) call ga_get(ga_s,il,ih,jl+15,jl+20,psp,msa)
      else
      call ga_get(ga_s,il,ih,jl+12,jl+17,ps,msa)
      if(lpp) call ga_get(ga_s,il,ih,jl+15,jl+20,psp,msa)
      endif
      endif
c
      return
      end
c $Id$
