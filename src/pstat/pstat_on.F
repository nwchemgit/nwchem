      Subroutine PStat_On( Handle )
C$Id: pstat_on.F,v 1.9 1995-10-10 06:14:47 d3g681 Exp $
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
CHF  $     ('pstat_on: Timer not assigned', Handle)
      if (.not. status) return
c
      If ( PStat_Unmatched_Starts( Handle ) .ne. 0  ) then
         Call ErrQuit('pstat_on: Timer is already on', Handle)
      Else
         PStat_Unmatched_Starts( Handle ) = 1
      EndIf
C
      If ( PStat_Functions( Handle, PStat_QStat) ) then
         PStat_Start_Times( Handle, 1 ) = util_cpusec()
         PStat_Start_Times( Handle, 2 ) = util_wallsec()
      EndIf
C
      Return
      End
