      subroutine argos_space_putix(node,iwl,xw,numwm,isl,xs,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
c
      call argos_space_gaputix(node,iwl,int_mb(i_packw),xw,numwm,
     + isl,int_mb(i_pack),xs,numsa,int_mb(i_ipl))
c
      return
      end
c $Id$
