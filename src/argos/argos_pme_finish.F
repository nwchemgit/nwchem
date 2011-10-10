      subroutine argos_pme_finish
c
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      if(.not.ga_destroy(ga_flg))
     + call md_abort('Failed to deallocate flg',0)
c
      call ga_destroy_list(mgz,lslab)
c
      if(.not.ma_pop_stack(l_slab))
     + call md_abort('Failed to deallocate slab',0)
      if(.not.ma_pop_stack(l_map2))
     + call md_abort('Failed to deallocate map2',0)
      if(.not.ma_pop_stack(l_map1))
     + call md_abort('Failed to deallocate map1',0)
      if(.not.ma_pop_stack(l_slmax))
     + call md_abort('Failed to deallocate slmax',0)
      if(.not.ma_pop_stack(l_lsld))
     + call md_abort('Failed to deallocate lsldo',0)
      if(.not.ma_pop_stack(l_knod))
     + call md_abort('Failed to deallocate knode',0)
      if(.not.ma_pop_stack(l_nodp))
     + call md_abort('Failed to deallocate nodep',0)
      if(.not.ma_pop_stack(l_kto))
     + call md_abort('Failed to deallocate kto',0)
      if(.not.ma_pop_stack(l_kfr))
     + call md_abort('Failed to deallocate kfrom',0)
      if(.not.ma_pop_stack(l_lnod))
     + call md_abort('Failed to deallocate lnode',0)
      if(.not.ma_pop_stack(l_grdh))
     + call md_abort('Failed to deallocate gridh',0)
      if(.not.ma_pop_stack(l_grdt))
     + call md_abort('Failed to deallocate gridt',0)
      if(.not.ma_pop_stack(l_grd))
     + call md_abort('Failed to deallocate grid',0)
      if(.not.ma_pop_stack(l_grd2))
     + call md_abort('Failed to deallocate grid2',0)
      if(.not.ma_pop_stack(l_grd1))
     + call md_abort('Failed to deallocate grid1',0)
c
      if(.not.ma_pop_stack(l_barr))
     + call md_abort('Failed to deallocate barr',0)
      if(.not.ma_pop_stack(l_bmod))
     + call md_abort('Failed to deallocate bmod',0)
c
      return
      end
c $Id$
