      subroutine argos_space_copyg(fw,fs)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa),fs(msa,3)
c
      integer il,ih,jl,jh
c
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      endif
c
      return
      end
c $Id$
