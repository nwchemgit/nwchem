      Subroutine PStat_Off( Handle )
C$Id: pstat_off.F,v 1.9 1995-10-10 06:14:47 d3g681 Exp $
      Implicit NONE
      Integer Handle
C
#include "pstatP.fh"
#include "util.fh"
C
      Double precision Zero
      Parameter (Zero = 0.0d0)
      logical status
C
      status = handle.gt.0 .and. handle.le.pstat_max_items
      if (status) status = pstat_active(handle)
CHF   if (.not. status) Call ErrQuit
CHF  $     ('pstat_off: Timer not assigned', Handle)
      if (.not. status) return
C
      If ( PStat_Unmatched_Starts( Handle ) .ne. 1  ) then
         write(6,*) ' pstat_off: timer is not on ', pstat_name(handle)
         call util_flush(6)
         Call ErrQuit('pstat_off: Timer is not on', Handle)
      Else
         PStat_Unmatched_Starts( Handle ) = 0
         PStat_Nr_Calls(Handle) = PStat_Nr_Calls(Handle) + 1
      EndIf
C
      If ( PStat_Functions( Handle, PStat_QStat) ) then
         PStat_Times( Handle, 1 ) = PStat_Times( Handle, 1 )
     $        + util_cpusec() - PStat_Start_Times( Handle, 1)
         PStat_Times( Handle, 2 ) = PStat_Times( Handle, 2 )
     $        + util_wallsec() - PStat_Start_Times( Handle, 2)
      EndIf
C
      Return
      End
