      subroutine argos_cafe_add_eangles(nsh,msh,msp,iang,
     + ang,rang,rub,esa)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nsh,msh,msp
      integer iang(msh,4)
      real*8 ang(msh,msp,mset),rang(msh,2),rub(msh,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nsh
      esa(iang(i,1),2)=esa(iang(i,1),2)+third*rang(i,2)
      esa(iang(i,2),2)=esa(iang(i,2),2)+third*rang(i,2)
      esa(iang(i,3),2)=esa(iang(i,3),2)+third*rang(i,2)
    1 continue
c
      if(iffld.eq.2) then
      do 2 i=1,nsh
      esa(iang(i,1),2)=esa(iang(i,1),2)+half*rub(i,2)
      esa(iang(i,3),2)=esa(iang(i,3),2)+half*rub(i,2)
    2 continue
      endif
c
      return
      end
c $Id$
