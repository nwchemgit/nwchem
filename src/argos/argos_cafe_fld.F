      subroutine argos_cafe_fld(numw,xw,xwm,nums,isat,isml,xs,xsm,rwx)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3),xsm(msm,3)
      real*8 rwx(mscr,3)
      integer isat(msa),isml(msa)
c
      integer iwm,iwa,isa,ism,ix,numw,nums
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
c     vacuo simulations do not need folding
c
      if(npbtyp.eq.0) return
c
      if(numw.gt.0) then
      call argos_cafe_pbc(0,xwm,mwm,rwx,mscr,0,1,numw)
      do 1 ix=1,3
      do 2 iwm=1,numw
      xwm(iwm,ix)=xwm(iwm,ix)-rwx(iwm,ix)
    2 continue
      do 3 iwa=1,nwa
      do 4 iwm=1,numw
      xw(iwm,ix,iwa)=xw(iwm,ix,iwa)-rwx(iwm,ix)
    4 continue
    3 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      call argos_cafe_pbc(0,xsm,nsm,rwx,mscr,0,1,nsm)
      do 5 ism=1,nsm
      do 6 ix=1,3
      xsm(ism,ix)=xsm(ism,ix)-rwx(ism,ix)
      do 7 isa=1,nums
      if(isml(isa).eq.ism) then
      xs(isa,ix)=xs(isa,ix)-rwx(ism,ix)
      endif
    7 continue
    6 continue
    5 continue
      call argos_cafe_pbc(0,xs,msa,rwx,mscr,0,1,nums)
      do 8 isa=1,nums
      if(isml(isa).eq.0) then
      xs(isa,1)=xs(isa,1)-rwx(isa,1)
      xs(isa,2)=xs(isa,2)-rwx(isa,2)
      xs(isa,3)=xs(isa,3)-rwx(isa,3)
      endif
    8 continue
      endif
c
      return
      end
c $Id$
