      Subroutine PStat_On( Handle )
      Implicit NONE
      Integer Handle
C
#include "pstatP.fh"
#include "tcgmsg.fh"
C
      Double precision Zero
      Parameter (Zero = 0.0d0)
C
      If ( Handle .gt. PStat_Assigned )
     $   Call ErrQuit('pstat_on: Timer not assigned', Handle)
C
      If ( PStat_Unmatched_Starts( Handle ) .ne. 0  ) then
         Call ErrQuit('pstat_on: Timer is already on', Handle)
      Else
         PStat_Unmatched_Starts( Handle ) = 1
      EndIf
C
      If ( PStat_Functions( Handle, PStat_QStat) ) then
         PStat_Start_Times( Handle, 1 ) = Zero
         PStat_Start_Times( Handle, 2 ) = TCGTime()
      EndIf
C
      Return
      End
