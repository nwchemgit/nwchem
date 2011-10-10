      subroutine argos_cafe_ndxxcl(id,idx,jdx,n)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,n
      integer idx(n),jdx(n)
c
      if(n.gt.mxt(id)) call md_abort('Error in ndxxcl',0)
c
      if(id.eq.1) then
      call argos_cafe_icopy(idx,jdx,n,int_mb(i_ixcl(id)),mxt(id))
      else
c      call argos_cafe_jcopy(idx,jdx,n,int_mb(i_ixcl(id)),mxt(id))
c      call argos_cafe_index(int_mb(i_ixcl(id)),mxt(id))
      call argos_cafe_ncopy(idx,jdx,n,int_mb(i_ixcl(id)),mxt(id))
      endif
c
      return
      end
c $Id$
