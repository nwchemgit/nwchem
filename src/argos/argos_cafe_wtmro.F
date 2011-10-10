      subroutine argos_cafe_wtmro(lfnmro,dera)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
c
      integer lfnmro
      real*8 dera(6,nsatot)
c
      integer i,j
c
      call ga_dgop(mcf_68,dera,5*nsatot,'+')
c
      if(me.eq.0) then
      write(lfnmro) ((dera(i,j),i=1,6),j=1,nsatot)
      else
      do 1 i=1,nsatot
      do 2 j=1,6
      dera(j,i)=zero
    2 continue
    1 continue
      endif
c
      return
      end
c $Id$
