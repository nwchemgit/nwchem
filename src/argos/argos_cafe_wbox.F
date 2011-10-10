      subroutine argos_cafe_wbox(numw,nums,iwat,isat,isml,ishop,wsm,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer numw,nums,iwat(mwa),isat(msa),isml(msa),ishop(msa)
      real*8 wsm(msm),wgt(mat,mset)
c
      integer i
c
      wbox=zero
c
      if(numw.gt.0) then
      do 1 i=1,nwa
      wbox=wbox+wgt(iwat(i),iset)
    1 continue
      wbox=dble(numw)*wbox
      endif
c
      do 2 i=1,msm
      wsm(i)=zero
    2 continue
c
      if(nums.gt.0) then
      do 3 i=1,nums
      if(iand(ishop(i),1).ne.1) then
      wbox=wbox+wgt(isat(i),iset)
      wsm(isml(i))=wsm(isml(i))+wgt(isat(i),iset)
      endif
    3 continue
      endif
c
      if(np.gt.1) then
      call ga_dgop(mcf_34,wbox,1,'+')
      if(msm.gt.0) call ga_dgop(mcf_35,wsm,msm,'+')
      endif
c
      return
      end
c $Id$
