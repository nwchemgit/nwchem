      subroutine argos_cafe_print_dstr(lfnout,snam,idrs,rdrs)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout
      character*16 snam(nsatot)
      integer idrs(ndrs,2)
      real*8 rdrs(ndrs,6)
c
      integer i
c
      if(me.eq.0.and.ndrs.gt.0) then
      write(lfnout,1001)
 1001 format(' Distance restraints',/)
      do 1 i=1,ndrs
      write(lfnout,1002) snam(idrs(i,1)),snam(idrs(i,2)),
     + rdrs(i,1),rdrs(i,2),rdrs(i,3),rdrs(i,4),rdrs(i,5)
 1002 format(1x,2a16,3f12.6,2e12.5)
    1 continue
      endif
c
      return
      end
c $Id$
