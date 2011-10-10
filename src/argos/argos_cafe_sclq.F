      subroutine argos_cafe_sclq(iqfr,iqto,chg)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer iqfr,iqto
      real*8 chg(mqt,mqp,mset)
c
      integer i,j
c
      do 1 i=iqfr,iqto
      do 2 j=1,nparms
      chg(i,1,j)=scaleq*chg(i,1,j)
    2 continue
    1 continue
c
      return
      end
c $Id$
