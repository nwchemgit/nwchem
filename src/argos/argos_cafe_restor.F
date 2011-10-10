      subroutine argos_cafe_restor(xw,yw,fw,vw,numw,xs,ys,fs,vs,nums)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),fw(mwm,3,mwa),vw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),fs(msa,3),vs(msa,3)
c
      integer i,j
c
      if(numw.gt.0) then
      do 1 j=1,nwa
      do 2 i=1,numw
      xw(i,1,j)=yw(i,1,j)
      xw(i,2,j)=yw(i,2,j)
      xw(i,3,j)=yw(i,3,j)
      fw(i,1,j)=vw(i,1,j)
      fw(i,2,j)=vw(i,2,j)
      fw(i,3,j)=vw(i,3,j)
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 4 i=1,nums
      xs(i,1)=ys(i,1)
      xs(i,2)=ys(i,2)
      xs(i,3)=ys(i,3)
      fs(i,1)=vs(i,1)
      fs(i,2)=vs(i,2)
      fs(i,3)=vs(i,3)
    4 continue
      endif
c
      return
      end
c $Id$
