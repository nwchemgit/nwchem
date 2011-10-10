      subroutine argos_prop_finish()
c
      implicit none
c
#include "argos_prop.fh"
#include "mafdecls.fh"
c
      if(ntype.eq.3.and.me.eq.0) then
      if(.not.ma_pop_stack(l_dfrm))
     + call md_abort('Failed to deallocate dfrm',me)
      if(.not.ma_pop_stack(l_dfr))
     + call md_abort('Failed to deallocate dfr',me)
      endif
c
      return
      end
c $Id$
