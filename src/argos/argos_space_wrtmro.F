      subroutine argos_space_wrtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,xw,vw,xwcr,isl,xs,vs,prjct)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnmro
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 stime,pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      character*80 prjct
c
      project=prjct
c
      call argos_space_wtmro(lfnmro,stime,pres,temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,xwcr,isl,int_mb(i_pack),xs,vs,
     + int_mb(i_ipl))
c
      return
      end
c $Id$
