      Logical Function PStat_Free( Handle )
C$Id: pstat_free.F,v 1.5 2003-07-03 22:37:30 edo Exp $
      Implicit NONE
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
     $     ('pstat_free: Timer not assigned', Handle)
      endif
c
      pstat_active(handle) = .false.
c
      PStat_Free = .TRUE.
c
      End
