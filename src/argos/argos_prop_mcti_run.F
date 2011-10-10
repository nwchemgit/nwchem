      subroutine argos_prop_mcti_run(rlambd,dlambd,ndec)
c
      implicit none
c
#include "argos_prop.fh"
#include "mafdecls.fh"
c
      real*8 rlambd,dlambd
      integer ndec
c
      call argos_prop_mcti_r(rlambd,dlambd,ndec,dbl_mb(i_dfr),
     + dbl_mb(i_dfrm),
     + psumt(16),psumt(117),psumt(118),psumt(119),psumt(120),
     + psumt(51),psumt(41))
c
      return
      end
c $Id$
