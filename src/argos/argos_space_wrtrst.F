      subroutine argos_space_wrtrst(lfnrst,filrst,lveloc,
     + pres,temp,tempw,temps,iwl,xw,vw,fw,xwcr,isl,xs,vs,fs,xscr,prjct,
     + lseq)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst
      character*255 filrst
      logical lveloc
      integer iwl(mwm,miw2),isl(msa,mis2),lseq(mseq)
      real*8 pres,temp,tempw,temps
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3),xscr(msm,3)
      character*80 prjct
c
      integer lenscr
c
      project=prjct
c
      lenscr=ma_inquire_avail(mt_byte)/
     + ((9*mwa+3)*ma_sizeof(mt_dbl,1,mt_byte)+
     + (mis2+4)*ma_sizeof(mt_int,1,mt_byte))-1
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bx',l_bx,i_bx))
     + call md_abort('Failed to allocate bx',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bv',l_bv,i_bv))
     + call md_abort('Failed to allocate bv',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3*mwa,'bf',l_bf,i_bf))
     + call md_abort('Failed to allocate bf',0)
      if(.not.ma_push_get(mt_dbl,lenscr*3,'br',l_br,i_br))
     + call md_abort('Failed to allocate br',0)
      if(.not.ma_push_get(mt_int,lenscr*max(mis2,2),'bi',l_bi,i_bi))
     + call md_abort('Failed to allocate bi',0)
      if(.not.ma_push_get(mt_int,lenscr,'n',l_n,i_n))
     + call md_abort('Failed to allocate n',0)
c
      call argos_space_wtrst(lfnrst,filrst,lveloc,pres,temp,tempw,temps,
     + iwl,int_mb(i_packw),xw,vw,fw,xwcr,isl,int_mb(i_pack),xs,vs,fs,
     + xscr,int_mb(i_ipl),lenscr,int_mb(i_bi),dbl_mb(i_bx),dbl_mb(i_bv),
     + dbl_mb(i_bf),dbl_mb(i_br),int_mb(i_bi),dbl_mb(i_bx),
     + dbl_mb(i_bv),dbl_mb(i_bf),lseq)
c
      if(.not.ma_pop_stack(l_n))
     + call md_abort('Failed to deallocate n',0)
      if(.not.ma_pop_stack(l_bi))
     + call md_abort('Failed to deallocate bi',0)
      if(.not.ma_pop_stack(l_br))
     + call md_abort('Failed to deallocate br',0)
      if(.not.ma_pop_stack(l_bf))
     + call md_abort('Failed to deallocate bf',0)
      if(.not.ma_pop_stack(l_bv))
     + call md_abort('Failed to deallocate bv',0)
      if(.not.ma_pop_stack(l_bx))
     + call md_abort('Failed to deallocate bx',0)
c
      return
      end
c $Id$
