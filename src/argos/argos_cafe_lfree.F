      subroutine argos_cafe_lfree()
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
c      if(.not.llst) call md_abort('Error 1 in argos_cafe_lfree',0)
c
      if(llst) then
      if(.not.ma_free_heap(l_list))
     + call md_abort('Failed to free list',0)
      endif
c
      llst=.false.
      llist=.false.
c
      return
      end
c $Id$
