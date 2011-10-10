      subroutine argos_space_free()
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      if(.not.ga_destroy(ga_isz))
     + call md_abort('Error ga_destroy isz',ga_isz)
      if(.not.ga_destroy(ga_iwz))
     + call md_abort('Error ga_destroy iwz',ga_iwz)
c
      if(.not.ga_destroy(ga_s))
     + call md_abort('Error ga_destroy s',ga_s)
      if(.not.ga_destroy(ga_is))
     + call md_abort('Error ga_destroy is',ga_is)
      if(.not.ga_destroy(ga_w))
     + call md_abort('Error ga_destroy w',ga_w)
      if(.not.ga_destroy(ga_iw))
     + call md_abort('Error ga_destroy iw',ga_iw)
      if(.not.ga_destroy(ga_ip))
     + call md_abort('Error ga_destroy iw',ga_ip)
c
      if(.not.ma_pop_stack(l_sndx))
     + call md_abort('Failed to deallocate sndx',0)
      if(.not.ma_pop_stack(l_lseq))
     + call md_abort('Failed to deallocate lseq',0)
      if(.not.ma_pop_stack(l_jpl))
     + call md_abort('Failed to deallocate jpl',0)
      if(.not.ma_pop_stack(l_ipl))
     + call md_abort('Failed to deallocate ipl',0)
      if(.not.ma_pop_stack(l_rng))
     + call md_abort('Failed to deallocate rng',0)
      if(.not.ma_pop_stack(l_boxr))
     + call md_abort('Failed to deallocate boxr',0)
      if(.not.ma_pop_stack(l_boxs))
     + call md_abort('Failed to deallocate boxs',0)
      if(.not.ma_pop_stack(l_iown))
     + call md_abort('Failed to deallocate iown',0)
      if(.not.ma_pop_stack(l_buren))
     + call md_abort('Failed to deallocate buren',0)
      if(.not.ma_pop_stack(l_bindex))
     + call md_abort('Failed to deallocate bindex',0)
c
      return
      end
c $Id$
