      subroutine md_finish()
c
c $Id: md_finish.F,v 1.16 2000-04-11 15:20:40 d3j191 Exp $
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
c     deallocate memory for coordinates, velocities, etc.
c     ---------------------------------------------------
c
c      call ma_summarize_allocated_blocks()
      call sp_pack_finish()
      if(ipolt.gt.0) then
      if(lpert2.or.lpert3) then
      if(.not.ma_pop_stack(l_psp))
     + call errquit('Failed to deallocate memory for psp',0)
      if(.not.ma_pop_stack(l_pwp))
     + call errquit('Failed to deallocate memory for pwp',0)
      endif
      if(.not.ma_pop_stack(l_ps))
     + call errquit('Failed to deallocate memory for ps',0)
      if(.not.ma_pop_stack(l_pw))
     + call errquit('Failed to deallocate memory for pw',0)
      endif
      if(iguide.gt.0) then
      if(.not.ma_pop_stack(l_gs))
     + call errquit('Failed to deallocate memory for gs',0)
      if(.not.ma_pop_stack(l_gw))
     + call errquit('Failed to deallocate memory for gw',0)
      endif
      if(.not.ma_pop_stack(l_esk))
     + call errquit('Failed to deallocate memory for esk',0)
      if(.not.ma_pop_stack(l_ess))
     + call errquit('Failed to deallocate memory for ess',0)
      if(.not.ma_pop_stack(l_esw))
     + call errquit('Failed to deallocate memory for esw',0)
      if(.not.ma_pop_stack(l_zs))
     + call errquit('Failed to deallocate memory for zs',0)
      if(.not.ma_pop_stack(l_dsr))
     + call errquit('Failed to deallocate memory for dsr',0)
      if(.not.ma_pop_stack(l_xscr))
     + call errquit('Failed to deallocate memory for xscr',0)
      if(.not.ma_pop_stack(l_xsm))
     + call errquit('Failed to deallocate memory for xsm',0)
      if(.not.ma_pop_stack(l_xwcr))
     + call errquit('Failed to deallocate memory for xwcr',0)
      if(.not.ma_pop_stack(l_fs))
     + call errquit('Failed to deallocate memory for fs',0)
      if(.not.ma_pop_stack(l_fw))
     + call errquit('Failed to deallocate memory for fw',0)
      if(.not.ma_pop_stack(l_vst))
     + call errquit('Failed to deallocate memory for vst',0)
      if(.not.ma_pop_stack(l_vwt))
     + call errquit('Failed to deallocate memory for vwt',0)
      if(.not.ma_pop_stack(l_vs))
     + call errquit('Failed to deallocate memory for vs',0)
      if(.not.ma_pop_stack(l_vw))
     + call errquit('Failed to deallocate memory for vw',0)
      if(.not.ma_pop_stack(l_ys))
     + call errquit('Failed to deallocate memory for ys',0)
      if(.not.ma_pop_stack(l_yw))
     + call errquit('Failed to deallocate memory for yw',0)
      if(.not.ma_pop_stack(l_xs))
     + call errquit('Failed to deallocate memory for xs',0)
      if(.not.ma_pop_stack(l_rtos))
     + call errquit('Failed to deallocate memory for rtos',0)
      if(.not.ma_pop_stack(l_xwm))
     + call errquit('Failed to deallocate memory for xwm',0)
      if(.not.ma_pop_stack(l_xw))
     + call errquit('Failed to deallocate memory for xw',0)
      if(.not.ma_pop_stack(l_isz))
     + call errquit('Failed to deallocate memory for isz',0)
      if(.not.ma_pop_stack(l_iwz))
     + call errquit('Failed to deallocate memory for iwz',0)
      if(.not.ma_pop_stack(l_is))
     + call errquit('Failed to deallocate memory for is',0)
      if(.not.ma_pop_stack(l_iw))
     + call errquit('Failed to deallocate memory for iw',0)
c
      call prp_finish()
c
      call cf_finish()
c
      if(.not.lqmd) call sp_finish()
c
      call md_timing
c
c     close files
c
      if(me.eq.0) then
      if(ntype.eq.1) then
      if(nfem.gt.0) then
      close(unit=lfnemt,status='keep')
      endif
      endif
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      close(unit=lfncoo,status='keep')
      endif
      if(nfprop.gt.0) then
      close(unit=lfnprp,status='keep')
      endif
      if(ntype.eq.3) then
      close(unit=lfngib,status='keep')
      endif
      if(nftime.gt.0) then
      close(unit=lfntim,status='keep')
      endif
      if(nfsync.gt.0) then
      close(unit=lfnsyn,status='keep')
      endif
      endif
c
      return
      end
      subroutine md_timing
c
      implicit none
c
#include "md_common.fh"
c
      integer timer_calls
      real*8 timer_cpu_average,timer_wall_average
      external timer_calls
      external timer_cpu_average,timer_wall_average
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' TIMING ANALYSIS',//,t20,' Number of calls',
     + t44,' Average times per call',
     + //,t40,'        Wall',t55,'         CPU',/)
      write(lfnout,1001) timer_calls(1),timer_wall_average(1),
     + timer_cpu_average(1)
 1001 format(' MD step',t20,i10,t40,f12.3,t55,f12.3)
c
      if(timer_calls(5).gt.0) then
      write(lfnout,2001) timer_calls(5),timer_wall_average(5),
     + timer_cpu_average(5)
 2001 format(' FFT',t20,i10,t40,f12.3,t55,f12.3)
      endif
      endif
c
      return
      end
