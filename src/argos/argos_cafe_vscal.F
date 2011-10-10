      subroutine argos_cafe_vscal(tgauss,
     + numw,nums,vw,vwt,vs,vst,iwdt,isdt,iwat,isat,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "util.fh"
#include "msgids.fh"
c
      real*8 tgauss
      integer numw,nums
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),vs(msa,3),vst(msa,3)
      real*8 wgt(mat,mset)
      integer iwdt(mwm),isdt(msa),iwat(mwa),isat(msa)
c
      real*8 ekwcg,ekscg,tempgw,tempgs,gaufcw,gaufcs,ekcg(2)
      integer iwm,iwa,ix,isa
c
      ekwcg=zero
      ekscg=zero
c
      if(numw.gt.0) then
      do 1 iwm=1,numw
      if(iwdt(iwm).ne.0) then
      do 2 iwa=1,nwa
      do 3 ix=1,3
      vw(iwm,ix,iwa)=zero
      vwt(iwm,ix,iwa)=zero
    3 continue
    2 continue
      else
      do 4 iwa=1,nwa
      ekwcg=ekwcg+wgt(iwat(iwa),iset)*(vw(iwm,1,iwa)*vw(iwm,1,iwa)
     + +vw(iwm,2,iwa)*vw(iwm,2,iwa)+vw(iwm,3,iwa)*vw(iwm,3,iwa))
    4 continue
      endif
    1 continue
      endif
c
      if(nums.gt.0) then
      do 5 isa=1,nums
      if(isdt(isa).ne.0) then
      do 6 ix=1,3
      vs(isa,ix)=zero
      vst(isa,ix)=zero
    6 continue
      else
      ekscg=ekscg+wgt(isat(isa),iset)*(vs(isa,1)*vs(isa,1)+
     + vs(isa,2)*vs(isa,2)+vs(isa,3)*vs(isa,3))
      endif
    5 continue
      endif
c
c     globally sum kinetic energies in ekcg
c
      ekcg(1)=half*ekwcg
      ekcg(2)=half*ekscg
      call ga_dgop(mcf_48,ekcg,2,'+')
      ekwcg=ekcg(1)
      ekscg=ekcg(2)
c
c      if(ivreas.ge.2) then
c      tempgw=fcwtmp*ekwcg
c      tempgs=fcstmp*ekscg
c      else
      tempgw=factmp*(ekwcg+ekscg)
      tempgs=tempgw
c      endif
c
      if(numw.gt.0.and.tempgw.gt.tiny) then
      gaufcw=sqrt(tgauss/tempgw)
      do 7 iwa=1,nwa
      do 8 ix=1,3
      do 9 iwm=1,numw
      vw(iwm,ix,iwa)=gaufcw*vw(iwm,ix,iwa)
      vwt(iwm,ix,iwa)=vw(iwm,ix,iwa)
    9 continue
    8 continue
    7 continue
      endif
c
      if(nums.gt.0.and.tempgs.gt.tiny) then
      gaufcs=sqrt(tgauss/tempgs)
      do 10 ix=1,3
      do 11 isa=1,nums
      vs(isa,ix)=gaufcs*vs(isa,ix)
      vst(isa,ix)=vs(isa,ix)
   11 continue
   10 continue
      endif
c
      return
      end
c $Id$
