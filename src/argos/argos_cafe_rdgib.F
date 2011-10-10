      subroutine argos_cafe_rdgib(lfngib,dec,dlambd)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfngib
      real*8 dec(6,nsatot),dlambd
c
      if(npgdec.eq.0) call md_abort('rdgib error',0)
c
      call argos_cafe_rgib(lfngib,dbl_mb(i_sti),dec,dlambd)
c
      return
      end
c $Id$
