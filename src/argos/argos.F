      logical function argos(irtdb0)
c
c $Id: nwmd.F,v 1.10 2003-08-07 00:24:45 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "argos.fh"
c
      integer irtdb0
c
      irtdb=irtdb0
c
      call util_print_push
      call util_print_rtdb_load(irtdb,'argos')
c
      call timer_init()
c
      call argos_rdinp()
c
      call argos_start()
      call argos_main()
      call argos_finish()
c
      call util_print_pop
c
      argos=.true.
c
      return
      end
