      Subroutine PStat_Off( Handle )
      Implicit NONE
      Integer Handle
C
#include "pstatP.fh"
#include "tcgmsg.fh"
#if defined(SUN) || defined(SGI) || defined(SGITFP)
#define ETIME
      real*4 etime, tarray(2)
      external etime
#endif
C
      Double precision Zero
      Parameter (Zero = 0.0d0)
      logical status
C
      status = handle.gt.0 .and. handle.le.pstat_max_items
      if (status) status = pstat_active(handle)
      if (.not. status) Call ErrQuit
     $     ('pstat_off: Timer not assigned', Handle)
C
      If ( PStat_Unmatched_Starts( Handle ) .ne. 1  ) then
         Call ErrQuit('pstat_off: Timer is not on', Handle)
      Else
         PStat_Unmatched_Starts( Handle ) = 0
         PStat_Nr_Calls(Handle) = PStat_Nr_Calls(Handle) + 1
      EndIf
C
      If ( PStat_Functions( Handle, PStat_QStat) ) then
#if defined(ETIME)
         PStat_Times( Handle, 1 ) = PStat_Times( Handle, 1 )
     $      + dble(etime(tarray)) - PStat_Start_Times( Handle, 1)
#else
         PStat_Times( Handle, 1 ) = PStat_Times( Handle, 1 )
     $      + Zero - PStat_Start_Times( Handle, 1)
#endif
         PStat_Times( Handle, 2 ) = PStat_Times( Handle, 2 )
     $      + TCGTime() - PStat_Start_Times( Handle, 2)
      EndIf
C
      Return
      End
