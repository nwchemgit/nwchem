      subroutine argos_space_travel(bx,xw,vw,xwcr,gw,iwl,numwm,
     + xs,vs,gs,isl,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),bx(3)
      real*8 gw(mwm,3,mwa),gs(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer numwm,numsa
c
      integer lenx,i_ndx,l_ndx,i_itmp,l_itmp,i_rtmp,l_rtmp
c
      call ga_sync()
c
      box(1)=bx(1)
      box(2)=bx(2)
      box(3)=bx(3)
c
      call argos_space_bscale(dbl_mb(i_boxs))
c
      call argos_space_owner(xw,iwl,nwmloc,xs,isl,nsaloc)
c
      lenx=max(nwm,nsa)
      if(.not.ma_push_get(mt_int,lenx,'ndx',l_ndx,i_ndx))
     + call md_abort('Failed to allocate ndx',0)
      if(.not.ma_push_get(mt_int,lenx,'itmp',l_itmp,i_itmp))
     + call md_abort('Failed to allocate itmp',0)
      if(.not.ma_push_get(mt_dbl,lenx,'rtmp',l_rtmp,i_rtmp))
     + call md_abort('Failed to allocate rtmp',0)
c
      call argos_space_trvl(xw,vw,xwcr,gw,iwl,int_mb(i_packw),
     + xs,vs,gs,isl,int_mb(i_pack),
     + dbl_mb(i_boxs),int_mb(i_iown),int_mb(i_ipl),
     + int_mb(i_ndx),int_mb(i_itmp),dbl_mb(i_rtmp),lenx)
c
      call ga_sync()
c
      if(.not.ma_pop_stack(l_rtmp))
     + call md_abort('Failed to de-allocate rtmp',0)
      if(.not.ma_pop_stack(l_itmp))
     + call md_abort('Failed to de-allocate itmp',0)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Failed to de-allocate ndx',0)
c
      call argos_space_gaputixv(me,iwl,int_mb(i_packw),
     + xw,xwcr,vw,nwmloc,
     + isl,int_mb(i_pack),xs,vs,nsaloc,int_mb(i_ipl))
c
      call argos_space_lbbl_indices()
c
      numwm=nwmloc
      numsa=nsaloc
c
      return
      end
c $Id$
