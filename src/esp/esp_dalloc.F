      subroutine esp_dalloc
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
c
      if(.not.ma_pop_stack(l_v))
     + call md_abort('esp: ma_pop_stack v failed',me)
      if(.not.ma_pop_stack(l_g))
     + call md_abort('esp: ma_pop_stack g failed',me)
      if(.not.ma_pop_stack(l_c))
     + call md_abort('esp: ma_pop_stack c failed',me)
      if(.not.ma_pop_stack(l_b))
     + call md_abort('esp: ma_pop_stack b failed',me)
      if(.not.ma_pop_stack(l_ak))
     + call md_abort('esp: ma_pop_stack Ak failed',me)
      if(.not.ma_pop_stack(l_a))
     + call md_abort('esp: ma_pop_stack A failed',me)
      if(.not.ma_pop_stack(l_o))
     + call md_abort('esp: ma_pop_stack o failed',me)
      if(.not.ma_pop_stack(l_e))
     + call md_abort('esp: ma_pop_stack e failed',me)
      if(.not.ma_pop_stack(l_d))
     + call md_abort('esp: ma_pop_stack d failed',me)
      if(.not.ma_pop_stack(l_l))
     + call md_abort('esp: ma_pop_stack l failed',me)
      if(.not.ma_pop_stack(l_n))
     + call md_abort('esp: ma_pop_stack n failed',me)
      if(.not.ma_pop_stack(l_qf))
     + call md_abort('esp: ma_pop_stack qf failed',me)
      if(.not.ma_pop_stack(l_qk))
     + call md_abort('esp: ma_pop_stack qk failed',me)
      if(.not.ma_pop_stack(l_q))
     + call md_abort('esp: ma_pop_stack q failed',me)
      if(.not.ma_pop_stack(l_x))
     + call md_abort('esp: ma_pop_stack x failed',me)
      if(.not.ma_pop_stack(l_t))
     + call md_abort('esp: ma_pop_stack t failed',me)
c
      if(.not.bas_destroy(ibasis))
     + call md_abort('esp: bas_destroy failed',me)
      if(.not.geom_destroy(igeom))
     + call md_abort('esp: geom_destroy failed',me)
c
      return
      end

