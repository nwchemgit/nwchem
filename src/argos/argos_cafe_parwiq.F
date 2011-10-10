      subroutine argos_cafe_parwiq(id,ia,iq)
c
c     argos_cafe_parwiq
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,ia,iq
c
      if(id.lt.1.or.id.gt.mwa) call md_abort('Error in parwiq',0)
      call argos_cafe_copwiq(id,ia,iq,int_mb(i_iwa),int_mb(i_iwq),
     + dbl_mb(i_mas),dbl_mb(i_chg),dbl_mb(i_vdw))
c
      return
      end
c $Id$
