      subroutine argos_cafe_lalloc()
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      if(llst) call md_abort('Error 1 in argos_cafe_lalloc',0)
c
      maxl=ma_inquire_avail(mt_int)/2
c
      if(.not.ma_alloc_get(mt_int,maxl,'list',l_list,i_list))
     + call md_abort('Failed to allocate list',0)
c
      llst=.true.
      llist=.false.
c
      return
      end
c $Id$
