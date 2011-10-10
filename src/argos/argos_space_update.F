      subroutine argos_space_update(node,vlatt,
     + iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),vlatt(3,3)
c
      integer i,j
c
      call argos_space_gaputixv(node,iwl,int_mb(i_packw),
     + xw,xwcr,vw,numwm,
     + isl,int_mb(i_pack),xs,vs,numsa,int_mb(i_ipl))
c
      do 1 j=1,3
      box(j)=vlat(j,j)
      boxh(j)=half*box(j)
      do 2 i=1,3
      vlat(i,j)=vlatt(i,j)
      vlati(i,j)=vlatt(i,j)
    2 continue
    1 continue
c
      call matinv(vlati,3,3)
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      endif
c
      return
      end
c $Id$
