      subroutine argos_cafe_add_ebonds(nsb,msb,msp,ibnd,bnd,rbnd,esa)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nsb,msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,mset),rbnd(msb,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nsb
      if(ibnd(i,3).eq.0) then
      esa(ibnd(i,1),2)=esa(ibnd(i,1),2)+half*rbnd(i,2)
      esa(ibnd(i,2),2)=esa(ibnd(i,2),2)+half*rbnd(i,2)
      endif
    1 continue
c
      return
      end
c $Id$
