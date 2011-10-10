      subroutine argos_diana_listtors(sgmnam,last)
c
      implicit none
c
#include "argos_diana_common.fh"
c
      integer last
      character*16 sgmnam(msa)
      integer i,j
c
      if(me.eq.0) then
      do 5 i=last,ntors
      write(ctors(i),1000) i,itors(i,1,3),(sgmnam(itors(i,j,1))(11:16),
     + ':',sgmnam(itors(i,j,1))(1:10),j=1,4)
 1000 format(2i4,4(a6,a,a10))
    5 continue
      endif
c
      return
      end
c $Id$
