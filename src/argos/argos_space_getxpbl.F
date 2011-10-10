      subroutine argos_space_getxpbl(ibbl,lhandl,
     + iwl,xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,xs,ps,psp,isfr,isto,jsfr,jsto,lself,lpbcs)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lself,lpbcs
c
      call argos_space_gtxpbl(ibbl,lhandl,
     + iwl,int_mb(i_packw),xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,int_mb(i_pack),xs,ps,psp,isfr,isto,jsfr,jsto,
     + int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),lself,.true.)
c
      lpbcs=lpbc9
c
      return
      end
c $Id$
