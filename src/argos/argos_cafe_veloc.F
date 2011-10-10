      subroutine argos_cafe_veloc(numw,xw,yw,vw,nums,xs,ys,vs)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),vw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3)
      integer numw,nums
c
      integer i,j,k
c
      if(numw.gt.0) then
      do 1 k=1,mwa
      do 2 j=1,3
      do 3 i=1,numw
      vw(i,j,k)=(xw(i,j,k)-yw(i,j,k))*tstepi
    3 continue
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 4 j=1,3
      do 5 i=1,nums
      vs(i,j)=(xs(i,j)-ys(i,j))*tstepi
    5 continue
    4 continue
      endif
c
      return
      end
       subroutine argos_cafe_get_msb(msb)
         implicit none
#include "argos_cafe_common.fh"
         integer msb
 
         msb=mbt(2)
 
       end
       subroutine argos_cafe_get_psb(psb)
         implicit none
#include "argos_cafe_common.fh"
         integer psb
 
         psb=i_ibnd(2)
 
       end
 
c $Id$
