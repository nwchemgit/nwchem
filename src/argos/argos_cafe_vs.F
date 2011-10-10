      subroutine argos_cafe_vs
     + (xs,ys,vs,vst,esk,nums,isat,ismol,isfrc,ishop,wsm,wgt,vsm)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isat(msa),ismol(msa),isfrc(msa),ishop(msa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3)
      real*8 wgt(mat,mset),vsm(msm,3)
      real*8 wsm(msm)
      real*8 esk(msf)
c
      integer i,j,k
      real*8 w
c
      do 11 i=1,msf
      esk(i)=zero
   11 continue
      if(nums.gt.0) then
      do 1 j=1,3
      if(mdalgo.eq.1) then
      do 2 i=1,nums
      vst(i,j)=half*(vs(i,j)+vst(i,j))
      if(iand(ishop(i),1).ne.1) esk(isfrc(i))=
     + esk(isfrc(i))+half*wgt(isat(i),iset)*vst(i,j)*vst(i,j)
    2 continue
      else
      do 3 i=1,nums
      if(iand(ishop(i),1).ne.1) esk(isfrc(i))=
     + esk(isfrc(i))+half*wgt(isat(i),iset)*vs(i,j)*vs(i,j)
    3 continue
      endif
    1 continue
      endif
c
      eskc(1,1)=zero
      eskc(2,1)=zero
      eskc(3,1)=zero
      eskc(1,2)=zero
      eskc(2,2)=zero
      eskc(3,2)=zero
      eskc(1,3)=zero
      eskc(2,3)=zero
      eskc(3,3)=zero
      do 5 i=1,msm
      vsm(i,1)=zero
      vsm(i,2)=zero
      vsm(i,3)=zero
    5 continue
      do 4 i=1,nums
      if(iand(ishop(i),1).ne.1) then
      k=ismol(i)
      w=wgt(isat(i),iset)
      if(mdalgo.eq.1) then
      do 6 j=1,3
      vsm(k,j)=vsm(k,j)+w*vst(i,j)
    6 continue
      else
      do 7 j=1,3
      vsm(k,j)=vsm(k,j)+w*vs(i,j)
    7 continue
      endif
      endif
    4 continue
      if(np.gt.1) call ga_dgop(mcf_43,vsm,3*msm,'+')
      if(me.eq.0) then
      do 8 i=1,nsm
      w=half/wsm(i)
      eskc(1,1)=eskc(1,1)+w*vsm(i,1)*vsm(i,1)
      eskc(2,1)=eskc(2,1)+w*vsm(i,2)*vsm(i,1)
      eskc(3,1)=eskc(3,1)+w*vsm(i,3)*vsm(i,1)
      eskc(1,2)=eskc(1,2)+w*vsm(i,1)*vsm(i,2)
      eskc(2,2)=eskc(2,2)+w*vsm(i,2)*vsm(i,2)
      eskc(3,2)=eskc(3,2)+w*vsm(i,3)*vsm(i,2)
      eskc(1,3)=eskc(1,3)+w*vsm(i,1)*vsm(i,3)
      eskc(2,3)=eskc(2,3)+w*vsm(i,2)*vsm(i,3)
      eskc(3,3)=eskc(3,3)+w*vsm(i,3)*vsm(i,3)
    8 continue
      endif
c
      return
      end
c $Id$
