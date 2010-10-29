      Logical Function PStat_Free( Handle )
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer Handle
C
#include "pstatP.fh"
C
      logical status
c
      status = handle.gt.0 .and. handle.le.pstat_max_items
      if (status) status = pstat_active(handle)
      if (.not. status) then
         write(0,*) pstat_name(handle)
         call errquit
     $     ('pstat_free: Timer not assigned', Handle, UNKNOWN_ERR)
      endif
c
      pstat_active(handle) = .false.
      pstat_start_flop(handle) = 0d0
      pstat_flop(handle) = 0d0
c
      PStat_Free = .TRUE.
c
      End
