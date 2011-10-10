      subroutine argos_cafe_links(ltype,nums,isgan,isdt,lst)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nums,ltype
      integer isgan(msa),isdt(msa),lst(nsa)
c
      call argos_cafe_cntlnk(ltype,mbt(2),int_mb(i_ibnd(2)),
     + nums,isgan,isdt,lst)
      return
      end
c $Id$
