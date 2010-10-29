      Logical Function PStat_Init(Max_Timers, NAcc, AccNames)
C$Id$
      Implicit NONE
      Integer Max_Timers, NAcc
      Character*(80) AccNames
C
#include "pstatP.fh"
C
      Double precision Zero
      Parameter (Zero = 0.0d0)
C
      Integer I, J
C
      If ( Max_Timers .gt. PStat_Max_Items ) then
         PStat_Init = .FALSE.
         Return
      EndIf
C
C     Clear out the names
C
      Do I = 1, PStat_Max_Items
         pstat_active(i) = .false. ! Mark as inactive
         PStat_Name( I ) = ' '
         Do J = 1, PStat_Nr_Stats
            PStat_Functions( I , J ) = .FALSE.
         EndDo
      EndDo
C
      Call IFill( PStat_Max_Items, 0, PStat_Nr_Calls, 1)
      Call IFill( PStat_Max_Items, 0, PStat_Unmatched_Starts, 1)
C
      Call DFill( 2*PStat_Max_Items, Zero, PStat_Times, 1)
      Call DFill( 2*PStat_Max_Items, Zero, PStat_Start_Times, 1)
C
      PStat_Init = .TRUE.
      Return
      End
