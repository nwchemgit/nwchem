      subroutine argos_cafe_cms(isat,isml,xs,xsm,gsm,isfr,isto,wsm,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer isfr,isto,isat(msa),isml(msa)
      real*8 xs(msa,3),xsm(msm,3),gsm(msm,4,2),wsm(msm),wgt(mat,mset)
c
      integer i,j
      real*8 dist2
c
      do 1 i=1,msm
      xsm(i,1)=zero
      xsm(i,2)=zero
      xsm(i,3)=zero
      gsm(i,1,1)=zero
      gsm(i,2,1)=zero
      gsm(i,3,1)=zero
      gsm(i,4,1)=zero
      gsm(i,1,2)=zero
      gsm(i,2,2)=zero
      gsm(i,3,2)=zero
      gsm(i,4,2)=zero
    1 continue
c
      if(isfr.gt.0.and.isfr.le.isto) then
c
      do 2 j=1,3
      do 3 i=isfr,isto
      xsm(isml(i),j)=xsm(isml(i),j)+wgt(isat(i),iset)*xs(i,j)
    3 continue
      do 4 i=1,msm
      xsm(i,j)=xsm(i,j)/wsm(i)
    4 continue
    2 continue
c
      endif
c
      if(np.gt.0) call ga_dgop(mcf_36,xsm,3*msm,'+')
c
      if(iradgy.gt.0) then
c
      if(isfr.gt.0.and.isfr.le.isto) then
c
      do 5 i=isfr,isto
      dist2=zero
      do 6 j=1,3
      dist2=dist2+(xs(i,j)-xsm(isml(i),j))**2
    6 continue
      do 7 j=1,3
      gsm(isml(i),j,1)=gsm(isml(i),j,1)+wgt(isat(i),j)*dist2
      gsm(isml(i),j,2)=gsm(isml(i),j,2)+wgt(isat(i),j)
    7 continue
      gsm(isml(i),4,1)=gsm(isml(i),4,1)+dist2
      gsm(isml(i),4,2)=gsm(isml(i),4,2)+1.0d0
    5 continue
c
      endif
c
      if(np.gt.0) call ga_dgop(mcf_67,gsm,8*msm,'+')
c
      do 9 i=1,nsm
      do 8 j=1,3
      gsm(i,j,2)=sqrt(gsm(i,j,1)/gsm(i,j,2))
    8 continue
      gsm(i,4,1)=sqrt(gsm(i,4,1)/gsm(i,4,2))
    9 continue
c
      endif
c
      return
      end
c $Id$
