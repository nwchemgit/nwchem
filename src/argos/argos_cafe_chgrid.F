      subroutine argos_cafe_chgrid(lset,ipset,xw,nwmloc,xs,
     + isq,isgm,nsaloc,
     + grid,mgz,lgrid,knode,kfrom,lsldo,node,mnx,mxx,mny,mxy,mnz,mxz)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      logical lset,lgrid
      integer ipset,nwmloc,nsaloc,node,mgz
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer isq(msa),isgm(msa)
      integer knode(*),kfrom(*)
      logical lsldo(*)
      real*8 grid(ngrx,ngry,mgz)
      integer mnx,mxx,mny,mxy,mnz,mxz
c
      if(lset.or..not.lscr) then
      if(.not.lscr) call argos_cafe_salloc
      call argos_cafe_pmeset(xw,dbl_mb(i_pmew),nwmloc,xs,dbl_mb(i_pmes),
     + nsaloc,dbl_mb(i_theta),dbl_mb(i_dtheta))
      lset=.false.
      endif
c
      call argos_cafe_grid(ipset,xw,dbl_mb(i_pmew),int_mb(i_iwq),nwmloc,
     + xs,dbl_mb(i_pmes),isq,isgm,nsaloc,
     + dbl_mb(i_chg),dbl_mb(i_theta),dbl_mb(i_dtheta),
     + grid,mgz,lgrid,knode,kfrom,lsldo,node,mnx,mxx,mny,mxy,mnz,mxz,
     + int_mb(i_lseq))
c
      return
      end
c $Id$
