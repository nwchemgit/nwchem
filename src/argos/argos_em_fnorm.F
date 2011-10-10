      subroutine argos_em_fnorm(iwl,isl,xw,xwm,xs,fw,fs)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
c
      integer i,j,k
c
      call timer_start(47)
      call argos_cafe_fnorm(iwl(1,lwdyn),fw,nwmloc,
     + isl(1,lsatt),isl(1,lsdyn),fs,nsaloc,fnorm,fmax)
      call timer_stop(47)
c
      return
      end
c $Id$
