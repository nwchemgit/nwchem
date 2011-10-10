      subroutine argos_diana_listdist(sgmnam,last)
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
      do 5 i=last,ndist
      write(cdist(i),1000) i,idist(i,1,3),(sgmnam(idist(i,j,1))(11:16),
     + ':',sgmnam(idist(i,j,1))(1:10),j=1,2)
 1000 format(2i4,4(a6,a,a10))
    5 continue
      endif
c
      return
      end
c $Id$
