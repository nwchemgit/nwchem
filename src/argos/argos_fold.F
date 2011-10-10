      subroutine argos_fold(iwl,isl,xw,xwm,xs,xsm)
c
      implicit none
c
#include "argos.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3),xsm(msm,3)
c
      call argos_cafe_fold(nwmloc,xw,xwm,nsaloc,isl(1,lsatt),isl(1,
     + lsmol),
     + xs,xsm)
c
      return
      end
c $Id$
