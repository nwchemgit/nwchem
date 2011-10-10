      subroutine argos_cafe_tw(vw,numw,iwat,wgt,ekinw)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer numw
      integer iwat(mwa)
      real*8 vw(mwm,3,mwa)
      real*8 wgt(mat,mset)
c
      integer i,j,k
      real*8 ekinwa,ekinw
c
      ekinw=zero
      if(numw.gt.0) then
      do 1 k=1,mwa
      ekinwa=zero
      do 2 j=1,3
      do 3 i=1,numw
      ekinwa=ekinwa+vw(i,j,k)*vw(i,j,k)
    3 continue
    2 continue
      ekinw=ekinw+wgt(iwat(k),iset)*ekinwa
    1 continue
      ekinw=half*ekinw
      endif
c
      return
      end
c $Id$
