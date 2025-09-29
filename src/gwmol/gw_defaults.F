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
#include "gw.fh"
c
      integer rtdb
      logical ok
c
c     Set default values and put them in RTDB
c
      ok = rtdb_put(rtdb,'gw:evgw0',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:evgw',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:cdgw',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:acgw',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:analytic',mt_log,1,.true.).and.
     &     rtdb_put(rtdb,'gw:analyticw',mt_log,1,.true.).and.
     &     rtdb_put(rtdb,'gw:newton',mt_log,1,.true.).and.
     &     rtdb_put(rtdb,'gw:graph',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:linear',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:eviter',mt_int,1,4).and.
     &     rtdb_put(rtdb,'gw:qpiter',mt_int,1,10).and.
     &     rtdb_put(rtdb,'gw:noqp_alpha',mt_int,1,1).and.
     &     rtdb_put(rtdb,'gw:nvqp_alpha',mt_int,1,0).and.
     &     rtdb_put(rtdb,'gw:noqp_beta',mt_int,1,1).and.
     &     rtdb_put(rtdb,'gw:nvqp_beta',mt_int,1,0).and.
     &     rtdb_put(rtdb,'gw:ngl',mt_int,1,200).and.
     &     rtdb_put(rtdb,'gw:nomega',mt_int,1,100).and.
     &     rtdb_put(rtdb,'gw:first',mt_int,1,1).and. 
     &     rtdb_put(rtdb,'gw:threshold',mt_dbl,1,0.005d0/ha2ev).and.
     &     rtdb_put(rtdb,'gw:ncap',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:eta',mt_dbl,1,0.0001d0).and.
     &     rtdb_put(rtdb,'gw:rpa',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:diag',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:core',mt_log,1,.false.).and.
     &     rtdb_put(rtdb,'gw:init',mt_log,1,.true.) 
c
      if (ok) return
c
      call errquit('gw_defaults','rtdb_put failed',1,RTDB_ERR)
c
      end subroutine gw_defaults
