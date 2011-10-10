      subroutine argos_cafe_slowr(xw,vw,numw,xs,vs,nums,iwat,isat,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),wgt(mat,mset)
      integer iwat(mwa),isat(msa)
c
      integer i,j,k
      real*8 abc(3),vbc(3),tbc(3,3),tbci(3,3),vbca(3)
      real*8 txx,txy,txz,tyy,tyz,tzz,wfac,dxt,dyt,dzt,deter,deteri
      real*8 drx,dry,drz
c
      do 1 i=1,3
      abc(i)=zero
      vbc(i)=zero
      vbca(i)=zero
    1 continue
c
      if(numw.gt.0) then
      do 2 k=1,nwa
      do 4 i=1,numw
      abc(1)=abc(1)+wgt(iwat(k),iset)
     + *(xw(i,2,k)*vw(i,3,k)-xw(i,3,k)*vw(i,2,k))
      abc(2)=abc(2)+wgt(iwat(k),iset)
     + *(xw(i,3,k)*vw(i,1,k)-xw(i,1,k)*vw(i,3,k))
      abc(3)=abc(3)+wgt(iwat(k),iset)
     + *(xw(i,1,k)*vw(i,2,k)-xw(i,2,k)*vw(i,1,k))
    4 continue
      do 5 j=1,3
      do 6 i=1,numw
      vbc(j)=vbc(j)+wgt(iwat(k),iset)*vw(i,j,k)
    6 continue
    5 continue
    2 continue
      endif
c
      if(nums.gt.0) then
      do 7 i=1,nums
      abc(1)=abc(1)+wgt(isat(i),iset)*(xs(i,2)*vs(i,3)-xs(i,3)*vs(i,2))
      abc(2)=abc(2)+wgt(isat(i),iset)*(xs(i,3)*vs(i,1)-xs(i,1)*vs(i,3))
      abc(3)=abc(3)+wgt(isat(i),iset)*(xs(i,1)*vs(i,2)-xs(i,2)*vs(i,1))
    7 continue
      do 8 j=1,3
      do 9 i=1,nums
      vbc(j)=vbc(j)+wgt(isat(i),iset)*vs(i,j)
    9 continue
    8 continue
      endif
c
      call ga_dgop(mcf_50,abc,3,'+')
      call ga_dgop(mcf_51,vbc,3,'+')
c
      vbc(1)=vbc(1)/wbox
      vbc(2)=vbc(2)/wbox
      vbc(3)=vbc(3)/wbox
c
      abc(1)=abc(1)-wbox*(xbm(2)*vbc(3)-xbm(3)*vbc(2))
      abc(2)=abc(2)-wbox*(xbm(3)*vbc(1)-xbm(1)*vbc(3))
      abc(3)=abc(3)-wbox*(xbm(1)*vbc(2)-xbm(2)*vbc(1))
c
c     evaluate the moment of inertia
c
      txx=zero
      txy=zero
      txz=zero
      tyy=zero
      tyz=zero
      tzz=zero
c
      if(numw.gt.0) then
      do 10 k=1,nwa
      wfac=wgt(iwat(k),iset)
      do 11 i=1,numw
      dxt=xw(i,1,k)-xbm(1)
      dyt=xw(i,2,k)-xbm(2)
      dzt=xw(i,3,k)-xbm(3)
      txx=txx+wfac*dxt*dxt
      txy=txy+wfac*dxt*dyt
      txz=txz+wfac*dxt*dzt
      tyy=tyy+wfac*dyt*dyt
      tyz=tyz+wfac*dyt*dzt
      tzz=tzz+wfac*dzt*dzt
   11 continue
   10 continue
      endif
c
      if(nums.gt.0) then
      do 12 i=1,nums
      wfac=wgt(isat(i),iset)
      dxt=xs(i,1)-xbm(1)
      dyt=xs(i,2)-xbm(2)
      dzt=xs(i,3)-xbm(3)
      txx=txx+wfac*dxt*dxt
      txy=txy+wfac*dxt*dyt
      txz=txz+wfac*dxt*dzt
      tyy=tyy+wfac*dyt*dyt
      tyz=tyz+wfac*dyt*dzt
      tzz=tzz+wfac*dzt*dzt
   12 continue
      endif
c
      tbc(1,1)=tyy+tzz
      tbc(2,2)=txx+tzz
      tbc(3,3)=txx+tyy
      tbc(1,2)=-txy
      tbc(2,1)=-txy
      tbc(1,3)=-txz
      tbc(3,1)=-txz
      tbc(2,3)=-tyz
      tbc(3,2)=-tyz
c
      call ga_dgop(mcf_52,tbc,9,'+')
c
c     inversion of the intertia tensor
c
      deter=tbc(1,1)*(tbc(2,2)*tbc(3,3)-tbc(3,2)*tbc(2,3))
     +     -tbc(2,1)*(tbc(1,2)*tbc(3,3)-tbc(3,2)*tbc(1,3))
     +     +tbc(3,1)*(tbc(1,2)*tbc(2,3)-tbc(2,2)*tbc(1,3))
      if(abs(deter).gt.small) then
      deteri=one/deter
      tbci(1,1)= deteri*(tbc(2,2)*tbc(3,3)-tbc(3,2)*tbc(2,3))
      tbci(1,2)=(-deteri)*(tbc(1,2)*tbc(3,3)-tbc(3,2)*tbc(1,3))
      tbci(1,3)= deteri*(tbc(1,2)*tbc(2,3)-tbc(2,2)*tbc(1,3))
      tbci(2,1)=(-deteri)*(tbc(2,1)*tbc(3,3)-tbc(3,1)*tbc(2,3))
      tbci(2,2)= deteri*(tbc(1,1)*tbc(3,3)-tbc(3,1)*tbc(1,3))
      tbci(2,3)=(-deteri)*(tbc(1,1)*tbc(2,3)-tbc(2,1)*tbc(1,3))
      tbci(3,1)= deteri*(tbc(2,1)*tbc(3,2)-tbc(3,1)*tbc(2,2))
      tbci(3,2)=(-deteri)*(tbc(1,1)*tbc(3,2)-tbc(3,1)*tbc(1,2))
      tbci(3,3)= deteri*(tbc(1,1)*tbc(2,2)-tbc(2,1)*tbc(1,2))
c
c     evaluate angular velocity center of mass
c
      do 13 i=1,3
      vbca(i)=zero
      do 14 j=1,3
      vbca(i)=vbca(i)+tbci(i,j)*abc(j)
   14 continue
   13 continue
      endif
c
c     remove the center of mass of box translation
c
      do 15 j=1,3
      if(numw.gt.0) then
      do 16 k=1,nwa
      do 17 i=1,numw
      vw(i,j,k)=vw(i,j,k)-vbc(j)
   17 continue
   16 continue
      endif
      if(nums.gt.0) then
      do 18 i=1,nums
      vs(i,j)=vs(i,j)-vbc(j)
   18 continue
      endif
      vbc(j)=zero
   15 continue
c      return
c
c     remove the center of mass of box rotation
c
      if(numw.gt.0) then
      do 19 k=1,nwa
      do 20 i=1,numw
      drx=xw(i,1,k)-xbm(1)
      dry=xw(i,2,k)-xbm(2)
      drz=xw(i,3,k)-xbm(3)
      vw(i,1,k)=vw(i,1,k)-vbca(2)*drz+vbca(3)*dry
      vw(i,2,k)=vw(i,2,k)-vbca(3)*drx+vbca(1)*drz
      vw(i,3,k)=vw(i,3,k)-vbca(1)*dry+vbca(2)*drx
   20 continue
   19 continue
      endif
      if(nums.gt.0) then
      do 21 i=1,nums
      drx=xs(i,1)-xbm(1)
      dry=xs(i,2)-xbm(2)
      drz=xs(i,3)-xbm(3)
      vs(i,1)=vs(i,1)-vbca(2)*drz+vbca(3)*dry
      vs(i,2)=vs(i,2)-vbca(3)*drx+vbca(1)*drz
      vs(i,3)=vs(i,3)-vbca(1)*dry+vbca(2)*drx
   21 continue
      endif
c
      return
      end
c $Id$
