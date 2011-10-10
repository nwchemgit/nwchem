      subroutine argos_cafe_collapse(ncoll,fcoll,nums,numw,ismol,issgm,
     + xs,xsm,mst,tsm,fs,xwm,xw,fw)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ncoll,numw,nums,mst
      integer ismol(msa),issgm(msa)
      real*8 fcoll
      real*8 xsm(msm,3),xs(msa,3),fs(msa,3),tsm(mst,4)
      real*8 xwm(mwm,3),xw(mwm,3,mwa),fw(mwm,3,mwa)
c
      integer i,j,k,ki,kf
      real*8 d
c
      ki=1
      kf=3
      if(ncoll.eq.1.or.ncoll.eq.4) ki=3
      if(ncoll.eq.2.or.ncoll.eq.5) kf=2
c
      if(ncoll.le.3) then
c
      do 5 i=1,msm
      if(ncoll.eq.1) then
      d=one/sqrt(xsm(i,1)**2+xsm(i,2)**2+xsm(i,3)**2)
      else
      d=one/sqrt(xsm(i,1)**2+xsm(i,2)**2)
      endif
      tsm(i,1)=xsm(i,1)*d
      tsm(i,2)=xsm(i,2)*d
      tsm(i,3)=xsm(i,3)*d
    5 continue
c
      do 1 k=ki,kf
      do 2 i=1,nums
      fs(i,k)=fs(i,k)-fcoll*tsm(ismol(i),k)
    2 continue
    1 continue
      do 4 i=1,numw
      if(ncoll.eq.1) then
      d=one/sqrt(xwm(i,1)**2+xwm(i,2)**2+xwm(i,3)**2)
      else
      d=one/sqrt(xwm(i,1)**2+xwm(i,2)**2)
      endif
      do 3 j=1,nwa
      do 6 k=ki,kf
      fw(i,k,j)=fw(i,k,j)-fcoll*d*xwm(i,k)
    6 continue
    3 continue
    4 continue
c
      else
c
      do 7 i=1,mst
      tsm(i,1)=0.0d0
      tsm(i,2)=0.0d0
      tsm(i,3)=0.0d0
      tsm(i,4)=0.0d0
    7 continue
      do 8 i=1,nums
      tsm(issgm(i),1)=tsm(issgm(i),1)+xs(i,1)
      tsm(issgm(i),2)=tsm(issgm(i),2)+xs(i,2)
      tsm(issgm(i),3)=tsm(issgm(i),3)+xs(i,3)
      tsm(issgm(i),4)=tsm(issgm(i),4)+1.0d0
    8 continue
      do 9 k=ki,kf
      do 10 i=1,nums
      fs(i,k)=fs(i,k)-fcoll*tsm(issgm(i),k)/tsm(issgm(i),4)
   10 continue
    9 continue
c
      endif
c
      return
      end
c $Id$
