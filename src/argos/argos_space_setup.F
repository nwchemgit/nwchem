      subroutine argos_space_setup(node,iwl,xw,xwcr,vw,gw,numwm,
     + isl,xs,xscr,vs,gs,numsa,lpack)
c
c $Id$
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),xscr(msm,3)
      real*8 gw(mwm,3,mwa),gs(msa,3)
      logical lpack
c
      integer i_qlst,l_qlst
c
      if(iguide.gt.0) then
      call argos_space_gagetixvf(node,iwl,int_mb(i_packw),
     + xw,xwcr,vw,gw,numwm,
     + isl,int_mb(i_pack),xs,vs,gs,numsa,int_mb(i_ipl))
      else
      call argos_space_gagetixv(node,iwl,int_mb(i_packw),
     + xw,xwcr,vw,numwm,
     + isl,int_mb(i_pack),xs,vs,numsa,int_mb(i_ipl))
      endif
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      call argos_space_qatoms(isl,numsa)
      if(.not.ma_push_get(mt_int,nsa,'qlst',l_qlst,i_qlst))
     + call md_abort('Failed to allocate memory for qlst',0)
      call argos_space_latoms(isl,int_mb(i_qlst),numsa)
      if(.not.ma_pop_stack(l_qlst))
     + call md_abort('Failed to deallocate memory for qlst',0)
      call argos_space_putix(me,iwl,xw,numwm,isl,xs,numsa)
      endif
c
      if(nsm.gt.0) call argos_space_getxscr(xscr,dbl_mb(i_xscr))
c
      return
      end
