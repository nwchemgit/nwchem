      subroutine argos_space_finish()
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      if(.not.ma_pop_stack(l_xscr))
     + call md_abort('Failed to deallocate xscr',0)
      if(.not.ma_pop_stack(l_bb))
     + call md_abort('Failed to deallocate memory for bb',0)
c
      call argos_space_free()
c
      return
      end
c $Id$
