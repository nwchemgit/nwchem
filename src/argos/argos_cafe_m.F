      subroutine argos_cafe_m(ww,ws,isat,nums,iwat,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nums,isat(msa)
      real*8 ww(mwa),ws(msa)
      integer iwat(msa)
      real*8 wgt(mat,mset)
c
      integer i
c
      do 1 i=1,nwa
      ww(i)=wgt(iwat(i),iset)
    1 continue
      do 2 i=1,nums
      ws(i)=wgt(isat(i),iset)
    2 continue
c
      return
      end
c $Id$
