      logical function argos_space_rdmri(lfnmri,stime,pres,
     + temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical argos_space_rmri
      external argos_space_rmri
c
      integer lfnmri
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
c
      argos_space_rdmri=argos_space_rmri(lfnmri,stime,pres,
     + temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,xwcr,isl,int_mb(i_pack),xs,vs,
     + int_mb(i_ipl))
c
      return
      end
c $Id$
