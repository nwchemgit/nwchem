      subroutine argos_cafe_volume()
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer i
c
c     recalculate half boxsize and volume
c
c     reactangular box type
c
      if(nbxtyp.eq.0) volume=box(1)*box(2)*box(3)
c
c     non-rectangular box type
c
      if(nbxtyp.eq.1) then
      do 1 i=1,3
      box(i)=sqrt(vlat(1,i)*vlat(1,i)+vlat(2,i)*vlat(2,i)+
     + vlat(3,i)*vlat(3,i))
    1 continue
      volume=vlat(1,1)*(vlat(2,2)*vlat(3,3)-vlat(2,3)*vlat(3,2))+
     +       vlat(1,2)*(vlat(2,3)*vlat(3,1)-vlat(2,1)*vlat(3,3))+
     +       vlat(1,3)*(vlat(2,1)*vlat(3,2)-vlat(2,2)*vlat(3,1))
      endif
c
      if(nbxtyp.eq.2) then
      volume=half*box(1)*box(2)*box(3)
      endif
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      return
      end
c $Id$
