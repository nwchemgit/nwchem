      subroutine argos_cafe_mcti(fdata)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 fdata(28)
c
      integer i
c
      do 1 i=1,24
      fdata(i)=deriv(i,1)
    1 continue
      fdata(25)=ep2(1)
      fdata(26)=ep3(1)
      fdata(27)=ep2(3)
      fdata(28)=ep3(3)
c
      return
      end
c $Id$
