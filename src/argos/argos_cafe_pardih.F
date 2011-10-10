      subroutine argos_cafe_pardih(id,idp,ip,p)
c
c     argos_cafe_pardih
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(6)
      real*8 p(4,24)
c
      if(id.le.0.or.id.gt.nbs) call md_abort('Error 1 in pardih',0)
      if(idp.lt.1.or.idp.gt.mdt(id))
     + call md_abort('Error 2 in pardih',0)
      if(id.eq.2.and.ip(5).lt.0) nmult(3)=nmult(3)+1
c
      call argos_cafe_pcopy(idp,ip,p,6,mdp(id),int_mb(i_idih(id)),
     + dbl_mb(i_dih(id)),mdt(id),
     + ip2(9+(id-1)*12),ip3(9+(id-1)*12),ith(9+(id-1)*12))
c
      return
      end
c $Id$
