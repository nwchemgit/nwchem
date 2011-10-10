      subroutine argos_cafe_vw(xw,yw,vw,vwt,numw,iwat,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer numw
      integer iwat(mwa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),vw(mwm,3,mwa),vwt(mwm,3,mwa)
      real*8 wgt(mat,mset)
c
      integer i,j,k,l
      real*8 vvsum,ww,visum,vjsum
c
      ewk=zero
      if(numw.gt.0) then
      do 1 k=1,mwa
      vvsum=zero
      if(mdalgo.eq.1) then
      do 2 j=1,3
      do 3 i=1,numw
      vwt(i,j,k)=half*(vw(i,j,k)+vwt(i,j,k))
      vvsum=vvsum+vwt(i,j,k)*vwt(i,j,k)
    3 continue
    2 continue
      else
      do 4 j=1,3
      do 5 i=1,numw
      vvsum=vvsum+vw(i,j,k)*vw(i,j,k)
    5 continue
    4 continue
      endif
      ewk=ewk+wgt(iwat(k),iset)*vvsum
    1 continue
      ewk=half*ewk
      endif
c
      ewkc(1,1)=zero
      ewkc(2,1)=zero
      ewkc(3,1)=zero
      ewkc(1,2)=zero
      ewkc(2,2)=zero
      ewkc(3,2)=zero
      ewkc(1,3)=zero
      ewkc(2,3)=zero
      ewkc(3,3)=zero
      if(numw.gt.0) then
      ww=zero
      do 6 i=1,mwa
      ww=ww+wgt(iwat(i),iset)
    6 continue
      if(mdalgo.eq.1) then
      do 17 l=1,3
      do 7 j=1,3
      do 8 i=1,numw
      visum=zero
      vjsum=zero
      do 9 k=1,mwa
      visum=visum+wgt(iwat(k),iset)*vwt(i,j,k)
      vjsum=vjsum+wgt(iwat(k),iset)*vwt(i,l,k)
    9 continue
      ewkc(l,j)=ewkc(l,j)+visum*vjsum
    8 continue
    7 continue
   17 continue
      else
      do 110 l=1,3
      do 10 j=1,3
      do 11 i=1,numw
      visum=zero
      vjsum=zero
      do 12 k=1,mwa
      visum=visum+wgt(iwat(k),iset)*vw(i,j,k)
      vjsum=vjsum+wgt(iwat(k),iset)*vw(i,l,k)
   12 continue
      ewkc(l,j)=ewkc(l,j)+visum*vjsum
   11 continue
   10 continue
  110 continue
      endif
      ewkc(1,1)=half*ewkc(1,1)/ww
      ewkc(2,1)=half*ewkc(2,1)/ww
      ewkc(3,1)=half*ewkc(3,1)/ww
      ewkc(1,2)=half*ewkc(1,2)/ww
      ewkc(2,2)=half*ewkc(2,2)/ww
      ewkc(3,2)=half*ewkc(3,2)/ww
      ewkc(1,3)=half*ewkc(1,3)/ww
      ewkc(2,3)=half*ewkc(2,3)/ww
      ewkc(3,3)=half*ewkc(3,3)/ww
      endif
c
      return
      end
c $Id$
