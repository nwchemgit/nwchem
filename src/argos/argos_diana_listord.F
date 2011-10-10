      subroutine argos_diana_listord(sgmnam,last,ndordr)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      integer last,ndordr
      character*16 sgmnam(msa)
      integer i,j
c
      if(me.eq.0) then
      do 5 i=last,nord
      write(cord(i),1000) (sgmnam(iord(i,j,1))(11:16),':',
     + sgmnam(iord(i,j,1))(1:10),j=1,ndordr+1)
 1000 format(4(a6,a1,a10))
    5 continue
      endif
c
      return
      end
 
c $Id$
