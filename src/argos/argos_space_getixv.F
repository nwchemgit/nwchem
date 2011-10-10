      subroutine argos_space_getixv(node,iwl,xw,xwcr,vw,numwm,
     + isl,xs,vs,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
c
      call argos_space_gagetixv(node,iwl,int_mb(i_packw),
     + xw,xwcr,vw,numwm,
     + isl,int_mb(i_pack),xs,vs,numsa,int_mb(i_ipl))
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      endif
c
      return
      end
c $Id$
