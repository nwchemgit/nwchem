      subroutine argos_cafe_parseq(id,idp)
c
c     argos_cafe_parseq
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,idp
c
      int_mb(i_mprot+id-1)=idp
c
      return
      end
c $Id$
