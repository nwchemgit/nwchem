      subroutine argos_space_getda(jnode,jsfr,jsto,
     + isl,isga,isq3,xs,ndxda)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer jnode,jsfr,jsto,ndxda
      integer isga(msa),isq3(msa),isl(msa,mis2)
      real*8 xs(msa,3)
c
      integer il,ih,jl,jh
      integer l_islt,i_islt
c
      if(.not.ma_push_get(mt_int,msa*mis2,'islt',l_islt,i_islt))
     + call md_abort('Failed to allocate memory for islt',0)
      call argos_space_gtda(jnode,jsfr,jsto,int_mb(i_islt),
     + int_mb(i_pack),isga,isq3,xs,ndxda)
      if(.not.ma_pop_stack(l_islt))
     + call md_abort('Failed to deallocate memory for islt',0)
c
      return
      end
c $Id$
