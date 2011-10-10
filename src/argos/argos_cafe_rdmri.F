      logical function argos_cafe_rdmri(lfnmri)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      logical argos_cafe_rmri
      external argos_cafe_rmri
c
      integer lfnmri
c
      if(npgdec.eq.0) call md_abort('rdmri error',0)
c
      argos_cafe_rdmri=argos_cafe_rmri(lfnmri,dbl_mb(i_sti))
c
      return
      end
c $Id$
