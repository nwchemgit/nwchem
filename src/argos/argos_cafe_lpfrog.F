      subroutine argos_cafe_lpfrog(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwat,wgt,iwgmn,ishp)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa),ishp(msa)
      integer iwgmn(mwm)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
      integer iwat(mwa)
      real*8 wgt(mat,mset)
c
      integer iwa,iwm,isa
      real*8 factor
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
      vwt(iwm,1,iwa)=vw(iwm,1,iwa)
      vwt(iwm,2,iwa)=vw(iwm,2,iwa)
      vwt(iwm,3,iwa)=vw(iwm,3,iwa)
      vw(iwm,1,iwa)=scalew*(vw(iwm,1,iwa)+fw(iwm,1,iwa)*factor)
      vw(iwm,2,iwa)=scalew*(vw(iwm,2,iwa)+fw(iwm,2,iwa)*factor)
      vw(iwm,3,iwa)=scalew*(vw(iwm,3,iwa)+fw(iwm,3,iwa)*factor)
      xw(iwm,1,iwa)=yw(iwm,1,iwa)+vw(iwm,1,iwa)*tstep
      xw(iwm,2,iwa)=yw(iwm,2,iwa)+vw(iwm,2,iwa)*tstep
      xw(iwm,3,iwa)=yw(iwm,3,iwa)+vw(iwm,3,iwa)*tstep
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
      vst(isa,1)=vs(isa,1)
      vst(isa,2)=vs(isa,2)
      vst(isa,3)=vs(isa,3)
      vs(isa,1)=scales*(vs(isa,1)+fs(isa,1)*factor)
      vs(isa,2)=scales*(vs(isa,2)+fs(isa,2)*factor)
      vs(isa,3)=scales*(vs(isa,3)+fs(isa,3)*factor)
      xs(isa,1)=xs(isa,1)+vs(isa,1)*tstep
      xs(isa,2)=xs(isa,2)+vs(isa,2)*tstep
      xs(isa,3)=xs(isa,3)+vs(isa,3)*tstep
      endif
c      write(lfnhop,222) isa,ishp(isa),ys(isa,1),ys(isa,2),ys(isa,3),
c     + xs(isa,1),xs(isa,2),xs(isa,3),fs(isa,1),fs(isa,2),fs(isa,3)
c  222 format(2i5,6f10.5,3e12.6)
    3 continue
c
      endif
c
      return
      end
c $Id$
