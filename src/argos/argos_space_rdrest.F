      logical function argos_space_rdrest(lfn,fil,boxsiz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer lfn
      character*255 fil
      real*8 boxsiz(maxbox,3)
c
      logical argos_space_rrest
      external argos_space_rrest
c
      integer i_ltmp,l_ltmp
c
      if(.not.ma_push_get(mt_int,(mbbl+1)*mbb2,'ltmp',l_ltmp,i_ltmp))
     + call md_abort('Failed to allocate memory for ltmp',0)
      argos_space_rdrest=argos_space_rrest(lfn,fil,int_mb(i_bb),
     + int_mb(i_ltmp),
     + mbbl+1,boxsiz)
      if(.not.ma_pop_stack(l_ltmp))
     + call md_abort('Failed to deallocate ltmp',0)
c
      return
      end
c $Id$
