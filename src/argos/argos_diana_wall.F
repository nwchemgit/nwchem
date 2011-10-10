      subroutine argos_diana_wall(iwsel,ival,ndxw)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer iwsel(mwm),ndxw(mwm)
      integer i,ival
c
      nwsel=0
      if(ival.ne.0) nwsel=nwm
      do 1 i=1,nwm
      iwsel(i)=ival
      ndxw(i)=i
    1 continue
c
      return
      end
c $Id$
