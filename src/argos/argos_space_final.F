      subroutine argos_space_final(fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,il,ih,jl,jh
c
      if(np.gt.0) then
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_acc(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm,one)
      if(llong) call ga_acc(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm,one)
      call ga_get(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      if(ltwin) call ga_get(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_acc(ga_s,il,ih,jl+6,jl+8,fs,msa,one)
      if(llong) call ga_acc(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa,one)
      call ga_get(ga_s,il,ih,jl+6,jl+8,fs,msa)
      if(ltwin) call ga_get(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa)
      endif
      endif
c
      if(lpair) then
      if(nwm.gt.0) then
      call ga_distribution(ga_iwz,me,il,ih,jl,jh)
      call ga_acc(ga_iwz,il,ih,1,1,iwz,mwm,1)
      call ga_get(ga_iwz,il,ih,1,1,iwz,mwm)
      do 1 i=1,nwmloc
      iwz(i)=min(1,iwz(i))
    1 continue
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_isz,me,il,ih,jl,jh)
      call ga_acc(ga_isz,il,ih,1,1,isz,msa,1)
      call ga_get(ga_isz,il,ih,1,1,isz,msa)
      do 2 i=1,nsaloc
      isz(i)=min(1,isz(i))
    2 continue
      endif
      endif
c
      return
      end
c $Id$
