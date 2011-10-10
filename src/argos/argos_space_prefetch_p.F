      subroutine argos_space_prefetch_p(n,iwl,xw,pw,pwp,isl,xs,ps,psp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer n
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
c
      n=nbbl
c
      call argos_space_lbbl_prefetch(int_mb(i_bb),int_mb(i_ipl),
     + int_mb(i_jpl),
     + iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),xs,ps,psp,
     + .true.)
c
      return
      end
c $Id$
