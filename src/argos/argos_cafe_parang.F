      subroutine argos_cafe_parang(id,idp,ip,p)
c
c     argos_cafe_parang
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(5)
      real*8 p(4,24)
c
      if(id.le.0.or.id.gt.nbs) call md_abort('Error 1 in parang',0)
      if(idp.lt.1.or.idp.gt.mht(id))
     + call md_abort('Error 2 in parang',0)
      if(id.eq.2.and.ip(4).lt.0) nmult(2)=nmult(2)+1
c
      call argos_cafe_pcopy(idp,ip,p,5,mhp(id),int_mb(i_iang(id)),
     + dbl_mb(i_ang(id)),mht(id),
     + ip2(8+(id-1)*12),ip3(8+(id-1)*12),ith(8+(id-1)*12))
c
      return
      end
c $Id$
