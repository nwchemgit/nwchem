      subroutine gw_defaults(rtdb)
c
c     Set defaults for GW@DFT
c
c     Daniel Mejia-Rodriguez
c     November 2020
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      logical ok
      double precision, parameter :: ha2ev = 27.211386245988d0
c
c     Set default values and put them in RTDB
c
      ok = rtdb_put(rtdb,'gw:evgw0',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:evgw',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:cdgw',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:acgw',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:analytic',mt_log,1,.true.).and.
     &     rtdb_put(rtdb,'gw:analyticw',mt_log,1,.true.).and.
     &     rtdb_put(rtdb,'gw:eviter',mt_int,1,4).and.
     &     rtdb_put(rtdb,'gw:qpiter',mt_int,1,10).and.
     &     rtdb_put(rtdb,'gw:noqp',mt_int,1,1).and.
     &     rtdb_put(rtdb,'gw:nvqp',mt_int,1,0).and.
     &     rtdb_put(rtdb,'gw:ngl',mt_int,1,200).and.
     &     rtdb_put(rtdb,'gw:nomega',mt_int,1,100).and.
     &     rtdb_put(rtdb,'gw:threshold',mt_dbl,1,0.01/ha2ev)
c
      if (ok) return
c
      call errquit('gw_defaults','rtdb_put failed',1,RTDB_ERR)
c
      end subroutine gw_defaults
