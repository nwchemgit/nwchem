C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     progress_report -- Give the user a progress report for time-
C     consuming operations
C
C REVISION
C     $Id: prog_report.F,v 1.2 1996-05-26 20:54:00 gg502 Exp $
C
C SYNOPSIS
      Subroutine Progress_Report(LUN, Msg)
      Implicit NONE
      Integer LUN
      Character*(*) Msg
C
C INCLUDE FILES
#include "global.fh"
C
C EXTERNAL ROUTINES
      Integer LNBlnk
      External LNBlnk
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Character*50 Date
C
      If (GA_NodeID() .eq. 0) then
         Date = ' '
         Call Util_Date( Date )
         Write (LUN, 9000) Msg, Date(:LNBlnk(Date))
         Call Util_Flush(LUN)
      EndIf
C
 9000 Format(1X, A, ' @', A)
C
      Return
      End
