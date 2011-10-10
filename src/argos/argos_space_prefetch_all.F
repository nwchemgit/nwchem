      subroutine argos_space_prefetch_all(n,iwl,xw,isl,xs)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer n
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(2),ps(2)
      real*8 pwp(2),psp(2)
c
      n=nbbl
c
      call argos_space_lbbl_prefetch(int_mb(i_bb),int_mb(i_ipl),
     + int_mb(i_jpl),
     + iwl,int_mb(i_packw),xw,pw,pwp,isl,int_mb(i_pack),xs,ps,psp,
     + .false.)
c
      icbbl=nbbl
c
      return
      end
c $Id$
