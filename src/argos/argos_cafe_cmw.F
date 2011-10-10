      subroutine argos_cafe_cmw(xw,xwm,iwfr,iwto,iwa,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer iwfr,iwto,iwa(mwa)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),wgt(mat,mset)
c
      integer i,j,k
      real*8 factor
c
      if(iwfr.eq.0.or.iwfr.gt.iwto) return
c
      factor=zero
      do 1 i=1,mwa
      factor=factor+wgt(iwa(i),iset)
    1 continue
      factor=one/factor
c
      do 2 j=1,3
      do 3 i=iwfr,iwto
      xwm(i,j)=zero
      do 4 k=1,mwa
      xwm(i,j)=xwm(i,j)+wgt(iwa(k),iset)*xw(i,j,k)
    4 continue
      xwm(i,j)=factor*xwm(i,j)
    3 continue
    2 continue
c
      return
      end
c $Id$
