      subroutine argos_cafe_scmf(icmopt,ismol,isatt,vs,fs,nums,fcm,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer icmopt,ismol(msa),isatt(msa),nums
      real*8 vs(msa,3),fs(msa,3),fcm(msm,3)
      real*8 wgt(mat,mset)
c
      integer i
c
      do 1 i=1,msm
      fcm(i,1)=zero
      fcm(i,2)=zero
      fcm(i,3)=zero
    1 continue
c
c     remove the center of mass velocity in the z-direction
c
      do 2 i=1,nums
      fcm(ismol(i),1)=fcm(ismol(i),1)+wgt(isatt(i),iset)*vs(i,3)
      fcm(ismol(i),2)=fcm(ismol(i),2)+wgt(isatt(i),iset)
      fcm(ismol(i),3)=fcm(ismol(i),3)+one
    2 continue
      if(np.gt.1) call ga_dgop(mcf_82,fcm,3*msm,'+')
      do 3 i=1,nums
      if(fcm(ismol(i),3).gt.one) then
      vs(i,3)=vs(i,3)-fcm(ismol(i),1)/fcm(ismol(i),2)
      endif
    3 continue
c
      do 4 i=1,msm
      fcm(i,1)=zero
      fcm(i,2)=zero
      fcm(i,3)=zero
    4 continue
c
c     remove the center of mass acceleration in the z-direction
c
      do 5 i=1,nums
      fcm(ismol(i),1)=fcm(ismol(i),1)+fs(i,3)
      fcm(ismol(i),2)=fcm(ismol(i),2)+wgt(isatt(i),iset)
      fcm(ismol(i),3)=fcm(ismol(i),3)+one
    5 continue
      if(np.gt.1) call ga_dgop(mcf_83,fcm,3*msm,'+')
      if(icmopt.eq.1) then
      do 6 i=1,nums
      if(fcm(ismol(i),3).gt.one) then
      fs(i,3)=fs(i,3)-wgt(isatt(i),iset)*fcm(ismol(i),1)/fcm(ismol(i),2)
      endif
    6 continue
      else
      do 7 i=1,nums
      if(fcm(ismol(i),3).gt.one) then
      fs(i,3)=1.0d-01*fs(i,3)
      endif
    7 continue
      endif
c
      return
      end
c $Id$
