      subroutine argos_cafe_add_eimprs(nso,mso,msp,ioop,oop,roop,esa)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nso,mso,msp
      integer ioop(mso,5)
      real*8 oop(mso,msp,mset),roop(mso,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nso
      esa(ioop(i,1),2)=esa(ioop(i,1),2)+fourth*roop(i,2)
      esa(ioop(i,2),2)=esa(ioop(i,2),2)+fourth*roop(i,2)
      esa(ioop(i,3),2)=esa(ioop(i,3),2)+fourth*roop(i,2)
      esa(ioop(i,4),2)=esa(ioop(i,4),2)+fourth*roop(i,2)
    1 continue
c
      return
      end
c $Id$
