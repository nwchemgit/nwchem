      subroutine argos_space_getxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + lself,lpbcs)
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
      logical lself,lpbcs
      real*8 pw(2),pwp(2),ps(2),psp(2)
c
      call argos_space_gtxpbl(ibbl,lhandl,
     + iwl,int_mb(i_packw),xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,int_mb(i_pack),xs,ps,psp,isfr,isto,jsfr,jsto,
     + int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),lself,.false.)
c
      lpbcs=lpbc9
c
      return
      end
c $Id$
