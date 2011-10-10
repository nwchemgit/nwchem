      subroutine argos_cafe_lampmf(rlam,slam,rpmf,lfnpmf)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 rlam(3),slam(3),rpmf(3,6,numpmf)
      integer lfnpmf
c
      integer i,j,k
c
      do 1 i=1,numpmf
      do 2 j=1,3
      do 3 k=1,3
      rpmf(j,k,i)=rlam(k)*rpmf(j,6,i)+slam(k)*rpmf(j,5,i)
    3 continue
    2 continue
      if(me.eq.0.and.mropt.ne.3) then
      write(-lfnpmf,1000) i,rpmf(1,1,i),rlam(1)
 1000 format(i3,2f12.6)
      endif
    1 continue
c
      return
      end
c $Id$
