      subroutine argos_cafe_add_edihes(nsd,msd,msp,idih,dih,rdih,esa)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nsd,msd,msp
      integer idih(msd,5)
      real*8 dih(msd,msp,mset),rdih(msd,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nsd
      esa(idih(i,1),2)=esa(idih(i,1),2)+fourth*rdih(i,2)
      esa(idih(i,2),2)=esa(idih(i,2),2)+fourth*rdih(i,2)
      esa(idih(i,3),2)=esa(idih(i,3),2)+fourth*rdih(i,2)
      esa(idih(i,4),2)=esa(idih(i,4),2)+fourth*rdih(i,2)
    1 continue
c
      return
      end
c $Id$
