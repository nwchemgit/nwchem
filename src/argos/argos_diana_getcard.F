      logical function argos_diana_getcard()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(me.eq.0) then
      read(lfncmd,1000,end=9999) card
 1000 format(a)
      endif
c
      argos_diana_getcard=.true.
      return
c
 9999 continue
      argos_diana_getcard=.false.
c
      return
      end
c $Id$
