      subroutine argos_cafe_dera_i(dera)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      real*8 dera(6,nsatot)
c
      integer i
c
      do 1 i=1,nsatot
      dera(1,i)=zero
      dera(2,i)=zero
      dera(3,i)=zero
      dera(4,i)=zero
      dera(5,i)=zero
      dera(6,i)=zero
    1 continue
c
      return
      end
c $Id$
