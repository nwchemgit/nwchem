      subroutine argos_cafe_lpfrog_bc(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwat,wgt,iwgmn,isgan)
c
c     leap frog variant by Brown and Clarke, Mol.Phys., 51, 1243 (1984)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa),isgan(msa),iwgmn(mwm)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
      integer iwat(mwa)
      real*8 wgt(mat,mset)
c
      integer iwa,iwm,isa
      real*8 factor,chi,ekin(2)
c
      if(numw.gt.0) then
      do 1 iwa=1,mwa
      factor=tstep/wgt(iwat(iwa),iset)
      do 2 iwm=1,numw
      yw(iwm,1,iwa)=xw(iwm,1,iwa)
      yw(iwm,2,iwa)=xw(iwm,2,iwa)
      yw(iwm,3,iwa)=xw(iwm,3,iwa)
      if(iand(iwdt(iwm),mfixed).eq.lfixed) then
      vw(iwm,1,iwa)=zero
      vw(iwm,2,iwa)=zero
      vw(iwm,3,iwa)=zero
      vwt(iwm,1,iwa)=zero
      vwt(iwm,2,iwa)=zero
      vwt(iwm,3,iwa)=zero
      else
      vwt(iwm,1,iwa)=vw(iwm,1,iwa)+half*fw(iwm,1,iwa)*factor
      vwt(iwm,2,iwa)=vw(iwm,2,iwa)+half*fw(iwm,2,iwa)*factor
      vwt(iwm,3,iwa)=vw(iwm,3,iwa)+half*fw(iwm,3,iwa)*factor
      endif
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 3 isa=1,nums
      factor=tstep/wgt(isat(isa),iset)
      ys(isa,1)=xs(isa,1)
      ys(isa,2)=xs(isa,2)
      ys(isa,3)=xs(isa,3)
      if(iand(isdt(isa),mfixed).eq.lfixed) then
      vs(isa,1)=zero
      vs(isa,2)=zero
      vs(isa,3)=zero
      vst(isa,1)=zero
      vst(isa,2)=zero
      vst(isa,3)=zero
      else
      vst(isa,1)=vs(isa,1)+half*fs(isa,1)*factor
      vst(isa,2)=vs(isa,2)+half*fs(isa,2)*factor
      vst(isa,3)=vs(isa,3)+half*fs(isa,3)*factor
      endif
    3 continue
      endif
c
      call argos_cafe_tw(vwt,numw,iwat,dbl_mb(i_mas),ekin(1))
      call argos_cafe_ts(vst,nums,isat,dbl_mb(i_mas),ekin(2))
c
      call ga_dgop(mcf_65,ekin,2,'+')
c
      chi=sqrt(tmpext/(factmp*(ekin(1)+ekin(2))))
c
      if(numw.gt.0) then
      do 4 iwa=1,mwa
      factor=tstep/wgt(iwat(iwa),iset)
      do 5 iwm=1,numw
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      vw(iwm,1,iwa)=(two*chi-one)*vw(iwm,1,iwa)+chi*fw(iwm,1,iwa)*factor
      vw(iwm,2,iwa)=(two*chi-one)*vw(iwm,2,iwa)+chi*fw(iwm,2,iwa)*factor
      vw(iwm,3,iwa)=(two*chi-one)*vw(iwm,3,iwa)+chi*fw(iwm,3,iwa)*factor
      xw(iwm,1,iwa)=yw(iwm,1,iwa)+vw(iwm,1,iwa)*tstep
      xw(iwm,2,iwa)=yw(iwm,2,iwa)+vw(iwm,2,iwa)*tstep
      xw(iwm,3,iwa)=yw(iwm,3,iwa)+vw(iwm,3,iwa)*tstep
      endif
    5 continue
    4 continue
      endif
c
      if(nums.gt.0) then
      do 6 isa=1,nums
      factor=tstep/wgt(isat(isa),iset)
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      vs(isa,1)=(two*chi-one)*vs(isa,1)+chi*fs(isa,1)*factor
      vs(isa,2)=(two*chi-one)*vs(isa,2)+chi*fs(isa,2)*factor
      vs(isa,3)=(two*chi-one)*vs(isa,3)+chi*fs(isa,3)*factor
      xs(isa,1)=xs(isa,1)+vs(isa,1)*tstep
      xs(isa,2)=xs(isa,2)+vs(isa,2)*tstep
      xs(isa,3)=xs(isa,3)+vs(isa,3)*tstep
      endif
    6 continue
      endif
c
      return
      end
c $Id$
