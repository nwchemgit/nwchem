      subroutine argos_cafe_parimp(id,idp,ip,p)
c
c     argos_cafe_parimp
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(6)
      real*8 p(4,24)
c
      if(id.le.0.or.id.gt.nbs) call md_abort('Error 1 in parimp',0)
      if(idp.lt.1.or.idp.gt.mit(id))
     + call md_abort('Error 2 in parimp',0)
      if(id.eq.2.and.ip(5).lt.0) nmult(4)=nmult(4)+1
c
      call argos_cafe_pcopy(idp,ip,p,6,mip(id),int_mb(i_iimp(id)),
     + dbl_mb(i_imp(id)),mit(id),
     + ip2(10+(id-1)*12),ip3(10+(id-1)*12),ith(10+(id-1)*12))
c
      return
      end
c $Id$
