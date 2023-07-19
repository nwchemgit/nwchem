      subroutine bse_defaults(rtdb)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

      integer rtdb
      logical ok

      ok = rtdb_put(rtdb,'bse:analytic',mt_log,1,.true.) .and.
     &     rtdb_put(rtdb,'bse:davidson',mt_log,1,.false.) .and. 
     &     rtdb_put(rtdb,'bse:lanczos',mt_log,1,.false.) .and. 
     &     rtdb_put(rtdb,'bse:tda',mt_log,1,.false.) .and. 
     &     rtdb_put(rtdb,'bse:ncap',mt_log,1,.false.) .and.
     &     rtdb_put(rtdb,'bse:ncapr',mt_log,1,.false.) .and.
     &     rtdb_put(rtdb,'bse:singlet',mt_log,1,.true.) .and.
     &     rtdb_put(rtdb,'bse:triplet',mt_log,1,.false.) .and.
     &     rtdb_put(rtdb,'bse:dowindow',mt_log,1,.false.) .and.
     &     rtdb_put(rtdb,'bse:elower',mt_dbl,1,0d0) .and.
     &     rtdb_put(rtdb,'bse:upper',mt_dbl,1,0d0) .and.
     &     rtdb_put(rtdb,'bse:nroots',mt_int,1,10) .and.
     &     rtdb_put(rtdb,'bse:init',mt_log,1,.true.)

      if (.not.ok)
     & call errquit('bse_defaults: could not put to rtdb',0,RTDB_ERR) 

      end subroutine bse_defaults
