      subroutine esp_dalloc
c
c $Id: esp_dalloc.F,v 1.5 2003-10-17 22:54:52 carlfahl Exp $
c
      implicit none
#include "errquit.fh"
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
c
      if(.not.ma_pop_stack(l_v))
     + call errquit('esp: ma_pop_stack v failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_g))
     + call errquit('esp: ma_pop_stack g failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_c))
     + call errquit('esp: ma_pop_stack c failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_b))
     + call errquit('esp: ma_pop_stack b failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_ak))
     + call errquit('esp: ma_pop_stack Ak failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_a))
     + call errquit('esp: ma_pop_stack A failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_o))
     + call errquit('esp: ma_pop_stack o failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_e))
     + call errquit('esp: ma_pop_stack e failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_d))
     + call errquit('esp: ma_pop_stack d failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_n))
     + call errquit('esp: ma_pop_stack n failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_qf))
     + call errquit('esp: ma_pop_stack qf failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_qk))
     + call errquit('esp: ma_pop_stack qk failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_q))
     + call errquit('esp: ma_pop_stack q failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_x))
     + call errquit('esp: ma_pop_stack x failed',me, MA_ERR)
      if(.not.ma_pop_stack(l_t))
     + call errquit('esp: ma_pop_stack t failed',me, MA_ERR)
c
      if(.not.bas_destroy(ibasis))
     + call errquit('esp: bas_destroy failed',me, BASIS_ERR)
      if(.not.geom_destroy(igeom))
     + call errquit('esp: geom_destroy failed',me, GEOM_ERR)
c
      return
      end

