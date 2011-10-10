      subroutine argos_cafe_ndxtrd(id,idx,jdx,n)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,n
      integer idx(n),jdx(n)
c
      if(n.gt.mtt(id)) call md_abort('Error in ndxtrd',0)
c
      if(id.eq.1) then
      call argos_cafe_icopy(idx,jdx,n,int_mb(i_itrd(id)),mtt(id))
      else
c      call argos_cafe_jcopy(idx,jdx,n,int_mb(i_itrd(id)),mtt(id))
c      call argos_cafe_index(int_mb(i_itrd(id)),mtt(id))
      call argos_cafe_ncopy(idx,jdx,n,int_mb(i_itrd(id)),mtt(id))
      endif
c
      return
      end
c $Id$
