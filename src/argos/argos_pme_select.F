      subroutine argos_pme_select(value)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
c
      real*8 value
      integer i_trgt,l_trgt
c
      if(.not.ma_push_get(mt_dbl,np,'target',l_trgt,i_trgt))
     + call md_abort('Failed to allocate target',0)
c
      call argos_pme_selnod(value,dbl_mb(i_trgt),int_mb(i_nodp),
     + int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),int_mb(i_map1))
c
      if(.not.ma_pop_stack(l_trgt))
     + call md_abort('Failed to deallocate target',0)
c
      return
      end
c $Id$
