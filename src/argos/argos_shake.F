      subroutine argos_shake(xw,yw,iwl,xs,ys,isl,dmax)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical argos_cafe_shakep
      external argos_cafe_shakep
c
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),xs(msa,3),ys(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 dmax
c
      integer i,j,lhandl,ibbl,iwfr,iwto,isfr,isto,nbbl
      logical lself
c
      dmax=zero
c
      if(nwmloc.gt.0) then
      call argos_cafe_shakew(xw,yw,iwl(1,lwgmn),iwl(1,lwdyn),nwmloc)
      do 1 j=1,nwa
      do 2 i=1,nwmloc
      dmax=max(dmax,(xw(i,1,j)-yw(i,1,j))**2+
     + (xw(i,2,j)-yw(i,2,j))**2+(xw(i,3,j)-yw(i,3,j))**2)
    2 continue
    1 continue
      endif
c
    3 continue
      if(nsaloc.gt.0) then
      call argos_space_nbbl(nbbl)
      do 4 ibbl=1,nbbl
      call argos_space_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
      if(lself) call argos_cafe_shakes(lhandl,xs,ys,isl(1,lsgan),isl(1,
     + lsatt),
     + isl(1,lssgm),isl(1,lsdyn),isl(1,lshop),isfr,isto)
    4 continue
      do 5 i=1,nsaloc
      dmax=max(dmax,(xs(i,1)-ys(i,1))**2+(xs(i,2)-ys(i,2))**2+
     + (xs(i,3)-ys(i,3))**2)
    5 continue
      endif
      if(npmf.eq.1.or.lpmfc) then
      if(.not.argos_cafe_shakep(xs,ys,isl(1,lsgan),isl(1,lsatt),isl(1,
     + lsdyn),
     + isl(1,lshop),nsaloc)) goto 3
      endif
c
      dmax=sqrt(dmax)
      call ga_dgop(mrg_d45,dmax,1,'max')
c
      return
      end
c $Id$
