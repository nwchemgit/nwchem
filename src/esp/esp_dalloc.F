      subroutine esp_dalloc
c
c $Id: esp_dalloc.F,v 1.4 1998-12-23 19:41:21 d3j191 Exp $
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
     + call errquit('esp: ma_pop_stack v failed',me)
      if(.not.ma_pop_stack(l_g))
     + call errquit('esp: ma_pop_stack g failed',me)
      if(.not.ma_pop_stack(l_c))
     + call errquit('esp: ma_pop_stack c failed',me)
      if(.not.ma_pop_stack(l_b))
     + call errquit('esp: ma_pop_stack b failed',me)
      if(.not.ma_pop_stack(l_ak))
     + call errquit('esp: ma_pop_stack Ak failed',me)
      if(.not.ma_pop_stack(l_a))
     + call errquit('esp: ma_pop_stack A failed',me)
      if(.not.ma_pop_stack(l_o))
     + call errquit('esp: ma_pop_stack o failed',me)
      if(.not.ma_pop_stack(l_e))
     + call errquit('esp: ma_pop_stack e failed',me)
      if(.not.ma_pop_stack(l_d))
     + call errquit('esp: ma_pop_stack d failed',me)
      if(.not.ma_pop_stack(l_n))
     + call errquit('esp: ma_pop_stack n failed',me)
      if(.not.ma_pop_stack(l_qf))
     + call errquit('esp: ma_pop_stack qf failed',me)
      if(.not.ma_pop_stack(l_qk))
     + call errquit('esp: ma_pop_stack qk failed',me)
      if(.not.ma_pop_stack(l_q))
     + call errquit('esp: ma_pop_stack q failed',me)
      if(.not.ma_pop_stack(l_x))
     + call errquit('esp: ma_pop_stack x failed',me)
      if(.not.ma_pop_stack(l_t))
     + call errquit('esp: ma_pop_stack t failed',me)
c
      if(.not.bas_destroy(ibasis))
     + call errquit('esp: bas_destroy failed',me)
      if(.not.geom_destroy(igeom))
     + call errquit('esp: geom_destroy failed',me)
c
      return
      end

