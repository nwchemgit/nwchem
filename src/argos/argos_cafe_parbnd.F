      subroutine argos_cafe_parbnd(id,idp,ip,p)
c
c     argos_cafe_parbnd
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(4)
      real*8 p(4,24)
c
      if(id.le.0.or.id.gt.nbs) call md_abort('Error in parbnd',0)
      if(idp.lt.1.or.idp.gt.mbt(id)) call md_abort('Error in parbnd',0)
      if(id.eq.2.and.ip(3).lt.0) nmult(1)=nmult(1)+1
c
      call argos_cafe_pcopy(idp,ip,p,4,mbp(id),int_mb(i_ibnd(id)),
     + dbl_mb(i_bnd(id)),mbt(id),
     + ip2(6+(id-1)*12),ip3(6+(id-1)*12),ith(6+(id-1)*12))
c
      return
      end
c $Id$
