      Subroutine PStat_On( Handle )
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
#if defined(IBM_USELESS)
c     This timer has too large a granularity to be useful
      double precision ibm_cputime
      external ibm_cputime
#endif
C
      Double precision Zero
      Parameter (Zero = 0.0d0)
      logical status
C
      status = handle.gt.0 .and. handle.le.pstat_max_items
      if (status) status = pstat_active(handle)
      if (.not. status) Call ErrQuit
     $     ('pstat_on: Timer not assigned', Handle)
c
      If ( PStat_Unmatched_Starts( Handle ) .ne. 0  ) then
         Call ErrQuit('pstat_on: Timer is already on', Handle)
      Else
         PStat_Unmatched_Starts( Handle ) = 1
      EndIf
C
      If ( PStat_Functions( Handle, PStat_QStat) ) then
#if defined(ETIME)
         PStat_Start_Times( Handle, 1 ) = dble(etime(tarray))
#elif defined(IBM_USELESS)
         PStat_Start_Times( Handle, 1 ) = ibm_cputime()
#else
         PStat_Start_Times( Handle, 1 ) = Zero
#endif
         PStat_Start_Times( Handle, 2 ) = TCGTime()
      EndIf
C
      Return
      End
