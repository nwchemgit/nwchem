      subroutine tce_multipole_destroy(d_mo_mp1,l_mo_mp1_offset,
     &                                 multipole)
c
c $Id$
c 
c Written by Jeff Hammond, January 2008.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
c#include "geom.fh"
c#include "rtdb.fh"
c#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
#include "tce_ints.fh"
c#include "tce_prop.fh"
      integer maxaxis,multipole

      maxaxis=0
      if (multipole.eq.1) then
        maxaxis=3
      elseif (multipole.eq.1) then
        maxaxis=9
      elseif (multipole.eq.1) then
        maxaxis=19
      endif

      do axis = maxaxis,1,-1
        call deletefile(d_mo_mp1(axis))
        if (.not.ma_pop_stack(l_mo_mp1_offset(axis)))
     1    call errquit("tce_multipole_destroy: MA problem",
     2                 l_mo_mp1_offset(axis),MA_ERR)
      enddo

      return
      end
