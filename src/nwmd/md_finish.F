      subroutine md_finish()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
c     deallocate memory for coordinates, velocities, etc.
c     ---------------------------------------------------
c
      if(.not.ma_pop_stack(l_ess))
     + call errquit('Failed to deallocate memory for ess',0)
      if(.not.ma_pop_stack(l_esw))
     + call errquit('Failed to deallocate memory for esw',0)
      if(.not.ma_pop_stack(l_zs))
     + call errquit('Failed to deallocate memory for zs',0)
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
      call cf_finish()
c
      call sp_finish()
c
      return
      end
