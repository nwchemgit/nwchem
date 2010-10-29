      logical function nwmd(irtdb0)
c
c $Id$
c
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "md_common.fh"
c
      integer irtdb0
c
      irtdb=irtdb0
c
      call util_print_push
      call util_print_rtdb_load(irtdb,'md')
c
      call timer_init()
c
      call md_rdinp()
c
      call md_start()
      call md_main()
      call md_finish()
c
      call util_print_pop
c
      nwmd=.true.
c
      return
      end



