      subroutine argos_space_gaputf(node,fw,numwm,fs,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer il,ih,jl,jh
c
      if(numwm.gt.0) then
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      endif
      if(numsa.gt.0) then
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      endif
c
      return
      end
c $Id$
