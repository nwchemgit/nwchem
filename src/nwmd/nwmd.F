      logical function nwmd(irtdb0)
c
c $Id: nwmd.F,v 1.3 1999-05-26 05:37:58 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "md_common.fh"
c
      integer irtdb0,mdold
      logical nwargos
      external nwargos
c
      irtdb=irtdb0
c
      mdold=1
      if(.not.rtdb_get(irtdb,'md:mdold',mt_int,1,mdold)) mdold=1
      if(mdold.ne.0) then
      nwmd=nwargos(irtdb0)
      return
      endif
c
      call timer_init()
c
      call md_start()
      call md_main()
      call md_finish()
c
      nwmd=.true.
c
      return
      end



