      Logical Function PStat_Free( Handle )
C$Id: pstat_free.F,v 1.4 1995-02-02 23:20:55 d3g681 Exp $
      Implicit NONE
      Integer Handle
C
#include "pstatP.fh"
C
      logical status
c
      status = handle.gt.0 .and. handle.le.pstat_max_items
      if (status) status = pstat_active(handle)
      if (.not. status) Call ErrQuit
     $     ('pstat_free: Timer not assigned', Handle)
c
      pstat_active(handle) = .false.
c
      PStat_Free = .TRUE.
c
      End
