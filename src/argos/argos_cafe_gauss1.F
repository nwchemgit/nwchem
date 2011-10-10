      subroutine argos_cafe_gauss1(tgauss,frgaus,
     + numw,nums,vw,vs,iwdt,isdt,iwat,isat,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "util.fh"
c
      real*8 tgauss,frgaus
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3),wgt(mat,mset)
      integer iwdt(mwm),isdt(msa),iwat(mwa),isat(msa)
c
      real*8 facdof,rran,tdev
      integer iwa,iwm,ix,isa,iran
c
      tdev=zero
c
      if(numw.gt.0) then
      facdof=rgas*tgauss
      do 1 iwa=1,nwa
      if(facdof.gt.zero) tdev=sqrt(facdof/wgt(iwat(iwa),iset))
      do 2 iwm=1,numw
      if(iwdt(iwm).ne.0) then
      vw(iwm,1,iwa)=zero
      vw(iwm,2,iwa)=zero
      vw(iwm,3,iwa)=zero
      else
      do 3 ix=1,3
      rran=zero
      do 4 iran=1,12
      rran=rran+util_random(0)
    4 continue
      if(facdof.gt.zero) then
      vw(iwm,ix,iwa)=(one-frgaus)*vw(iwm,ix,iwa)+frgaus*(rran-six)*tdev
      else
      vw(iwm,ix,iwa)=zero
      endif
    3 continue
      endif
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      facdof=rgas*tgauss
      do 5 isa=1,nums
      if(isdt(isa).ne.0) then
      vs(isa,1)=zero
      vs(isa,2)=zero
      vs(isa,3)=zero
      else
      do 6 ix=1,3
      rran=zero
      do 7 iran=1,12
      rran=rran+util_random(0)
    7 continue
      if(facdof.gt.zero) then
      vs(isa,ix)=(one-frgaus)*vs(isa,ix)+
     + frgaus*(rran-six)*sqrt(facdof/wgt(isat(isa),iset))
      else
      vs(isa,ix)=zero
      endif
    6 continue
      endif
    5 continue
      endif
c
      return
      end
c $Id$
