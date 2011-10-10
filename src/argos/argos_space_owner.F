      subroutine argos_space_owner(xw,iwl,numwm,xs,isl,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer numwm,numsa
c
c     determine solvent ownership
c
      call argos_space_ownerw(xw,iwl,numwm,dbl_mb(i_boxs),
     + int_mb(i_iown))
c
c     determine solute ownership
c
      call argos_space_owners(xs,isl,numsa,dbl_mb(i_boxs),
     + int_mb(i_iown))
c
      return
      end
c $Id$
