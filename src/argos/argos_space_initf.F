      subroutine argos_space_initf(fw,fs,llng,iwz,isz,lpair)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer iwz(mwm),isz(msa)
      logical llng,lpair
c
      integer i,j,k,l,m,il,ih,jl,jh
c
      llong=llng
c
      m=1
      if(llong) m=2
c
      do 1 l=1,m
      if(nwm.gt.0) then
      do 2 k=1,mwa
      do 3 j=1,3
      do 4 i=1,mwm
      fw(i,j,k,l)=zero
    4 continue
    3 continue
    2 continue
      endif
      if(nsa.gt.0) then
      do 5 j=1,3
      do 6 i=1,msa
      fs(i,j,l)=zero
    6 continue
    5 continue
      endif
    1 continue
c
      if(nwm.gt.0) then
      call ga_distribution(ga_w,me,il,ih,jl,jh)
      call ga_put(ga_w,il,ih,jl+6*mwa+3,jl+9*mwa+2,fw,mwm)
      if(llong) call ga_put(ga_w,il,ih,jl+9*mwa+3,jl+12*mwa+2,
     + fw(1,1,1,2),mwm)
      endif
      if(nsa.gt.0) then
      call ga_distribution(ga_s,me,il,ih,jl,jh)
      call ga_put(ga_s,il,ih,jl+6,jl+8,fs,msa)
      if(llong) call ga_put(ga_s,il,ih,jl+9,jl+11,fs(1,1,2),msa)
      endif
c
      if(lpair) then
      do 7 i=1,mwm
      iwz(i)=0
    7 continue
      do 8 i=1,msa
      isz(i)=0
    8 continue
      call ga_zero(ga_iwz)
      call ga_zero(ga_isz)
      endif
c
      return
      end
c $Id$
