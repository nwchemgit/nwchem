      subroutine argos_finish()
c
c $Id$
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
c
c     deallocate memory for coordinates, velocities, etc.
c     ---------------------------------------------------
c
c      call ma_summarize_allocated_blocks()
      call argos_space_pack_finish()
      if(ipolt.gt.0) then
      if(lpert2.or.lpert3) then
      if(.not.ma_pop_stack(l_psp))
     + call md_abort('Failed to deallocate memory for psp',0)
      if(.not.ma_pop_stack(l_pwp))
     + call md_abort('Failed to deallocate memory for pwp',0)
      endif
      if(.not.ma_pop_stack(l_ps))
     + call md_abort('Failed to deallocate memory for ps',0)
      if(.not.ma_pop_stack(l_pw))
     + call md_abort('Failed to deallocate memory for pw',0)
      endif
c
      if(imembr.gt.0) then
      if(.not.ma_pop_stack(l_fm))
     + call md_abort('Failed to deallocate memory for fm',me)
      if(.not.ma_pop_stack(l_mm))
     + call md_abort('Failed to deallocate memory for mm',me)
      endif
      if(icmopt.gt.0) then
      if(.not.ma_pop_stack(l_fcm))
     + call md_abort('Failed to deallocate memory for fcm',0)
      endif
      if(iguide.gt.0) then
      if(.not.ma_pop_stack(l_gs))
     + call md_abort('Failed to deallocate memory for gs',0)
      if(.not.ma_pop_stack(l_gw))
     + call md_abort('Failed to deallocate memory for gw',0)
      endif
      if(.not.ma_pop_stack(l_esa))
     + call md_abort('Failed to deallocate memory for esa',0)
      if(.not.ma_pop_stack(l_wws))
     + call md_abort('Failed to deallocate memory for wws',0)
      if(.not.ma_pop_stack(l_esk))
     + call md_abort('Failed to deallocate memory for esk',0)
      if(.not.ma_pop_stack(l_fss))
     + call md_abort('Failed to deallocate memory for fss',0)
      if(.not.ma_pop_stack(l_ess))
     + call md_abort('Failed to deallocate memory for ess',0)
      if(.not.ma_pop_stack(l_esw))
     + call md_abort('Failed to deallocate memory for esw',0)
      if(.not.ma_pop_stack(l_zs))
     + call md_abort('Failed to deallocate memory for zs',0)
      if(.not.ma_pop_stack(l_dsr))
     + call md_abort('Failed to deallocate memory for dsr',0)
      if(.not.ma_pop_stack(l_xscr))
     + call md_abort('Failed to deallocate memory for xscr',0)
      if(.not.ma_pop_stack(l_gsm))
     + call md_abort('Failed to deallocate memory for gsm',0)
      if(.not.ma_pop_stack(l_xsmp))
     + call md_abort('Failed to deallocate memory for xsmp',0)
      if(.not.ma_pop_stack(l_tsm))
     + call md_abort('Failed to deallocate memory for tsm',0)
      if(.not.ma_pop_stack(l_xsm))
     + call md_abort('Failed to deallocate memory for xsm',0)
      if(.not.ma_pop_stack(l_xwcr))
     + call md_abort('Failed to deallocate memory for xwcr',0)
      if(.not.ma_pop_stack(l_fs))
     + call md_abort('Failed to deallocate memory for fs',0)
      if(.not.ma_pop_stack(l_fw))
     + call md_abort('Failed to deallocate memory for fw',0)
      if(.not.ma_pop_stack(l_vst))
     + call md_abort('Failed to deallocate memory for vst',0)
      if(.not.ma_pop_stack(l_vwt))
     + call md_abort('Failed to deallocate memory for vwt',0)
      if(.not.ma_pop_stack(l_vs))
     + call md_abort('Failed to deallocate memory for vs',0)
      if(.not.ma_pop_stack(l_vw))
     + call md_abort('Failed to deallocate memory for vw',0)
      if(.not.ma_pop_stack(l_ys))
     + call md_abort('Failed to deallocate memory for ys',0)
      if(.not.ma_pop_stack(l_yw))
     + call md_abort('Failed to deallocate memory for yw',0)
      if(.not.ma_pop_stack(l_xs))
     + call md_abort('Failed to deallocate memory for xs',0)
      if(.not.ma_pop_stack(l_rtos))
     + call md_abort('Failed to deallocate memory for rtos',0)
      if(.not.ma_pop_stack(l_xwm))
     + call md_abort('Failed to deallocate memory for xwm',0)
      if(.not.ma_pop_stack(l_xw))
     + call md_abort('Failed to deallocate memory for xw',0)
      if(.not.ma_pop_stack(l_isz))
     + call md_abort('Failed to deallocate memory for isz',0)
      if(.not.ma_pop_stack(l_iwz))
     + call md_abort('Failed to deallocate memory for iwz',0)
      if(.not.ma_pop_stack(l_is))
     + call md_abort('Failed to deallocate memory for is',0)
      if(.not.ma_pop_stack(l_iw))
     + call md_abort('Failed to deallocate memory for iw',0)
c
      call argos_prop_finish()
c
      call argos_cafe_finish()
c
      call argos_space_finish()
c
      call argos_timing
c
c     close files
c
      if(me.eq.0) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      close(unit=lfntrj,status='keep')
      endif
      if(nfhop.gt.0) then
      close(unit=lfnhop,status='keep')
      endif
      if(nfprop.gt.0) then
      close(unit=lfnprp,status='keep')
      endif
      if(ntype.eq.3) then
      close(unit=lfngib,status='keep')
      endif
      if(ntype.eq.0.and.nftri.gt.0) then
      close(unit=lfntri,status='keep')
      endif
      if(ntype.eq.2.and.iprpmf.ne.0) then
      close(unit=iabs(lfnpmf),status='keep')
      endif
      if(nftime.gt.0) then
      close(unit=lfntim,status='keep')
      endif
      if(nfsync.gt.0) then
      close(unit=lfnsyn,status='keep')
      endif
      if(lfnpmf.lt.0) then
      close(unit=-lfnpmf,status='keep')
      endif
      endif
c
      if(npg.gt.1) call argos_partition_end()
c
      return
      end
