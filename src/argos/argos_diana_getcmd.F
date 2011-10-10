      subroutine argos_diana_getcmd()
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
      lesd=.false.
c
      if(me.eq.0) then
      read(lfncmd,1000,end=9999) card
 1000 format(a)
      endif
      call util_char_ga_brdcst(mag_d00,card,0)
      cmd=card(1:6)
c
      return
 9999 continue
      call md_abort('Error readin command file',0)
      return
      end
c $Id$
