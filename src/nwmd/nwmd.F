      logical function nwmd(irtdb0)
c
c $Id: nwmd.F,v 1.8 2001-05-09 17:51:23 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "md_common.fh"
c
      integer irtdb0,mdold
c
      irtdb=irtdb0
c
      call timer_init()
c
      call md_rdinp()
c
      call md_start()
      call md_main()
      call md_finish()
c
      nwmd=.true.
c
      return
      end



