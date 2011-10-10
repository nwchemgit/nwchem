      subroutine argos_space_wtrest(lfn)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer lfn
c
      integer i_ltmp,l_ltmp
c
      if(.not.ma_push_get(mt_int,(mbbl+1)*mbb2,'ltmp',l_ltmp,i_ltmp))
     + call md_abort('Failed to allocate memory for ltmp',0)
      call argos_space_wrest(lfn,int_mb(i_bb),int_mb(i_ltmp),mbbl+1,
     + dbl_mb(i_boxs))
      if(.not.ma_pop_stack(l_ltmp))
     + call md_abort('Failed to deallocate ltmp',0)
c
      return
      end
c $Id$
