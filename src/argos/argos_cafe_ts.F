      subroutine argos_cafe_ts(vs,nums,isat,wgt,ekins)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isat(msa)
      real*8 vs(msa,3),wgt(mat,mset)
c
      integer i,j
      real*8 ekins
c
      ekins=zero
      if(nums.gt.0) then
      do 1 j=1,3
      do 2 i=1,nums
      ekins=ekins+wgt(isat(i),iset)*vs(i,j)*vs(i,j)
    2 continue
    1 continue
      ekins=half*ekins
      endif
c
      return
      end
c $Id$
