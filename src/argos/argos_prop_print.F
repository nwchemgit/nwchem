      subroutine argos_prop_print()
c
      implicit none
c
#include "argos_prop.fh"
c
      character*10 pdate,ptime
      integer i,j
c
      if(me.ne.0) return
c
      call swatch(pdate,ptime)
      write(lfnout,1000) pdate,ptime
 1000 format(/,' SINGLE POINT PROPERTIES',t110,2a10,/)
      do 1  i=1,nprop
      j=ixp(i)
      if(lp(j)) write(lfnout,1001) pronam(j)(1:39),p(j),pronam(j)(40:50)
 1001 format(1x,a39,t41,1pe18.9,1x,a11)
    1 continue
c
      return
      end
c $Id$
