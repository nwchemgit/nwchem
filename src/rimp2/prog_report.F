C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     progress_report -- Give the user a progress report for time-
C     consuming operations
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Progress_Report(LUN, Msg)
      Implicit NONE
      Integer LUN
      Character*(*) Msg
C
C INCLUDE FILES
#include "global.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C PARAMETERS
      Integer LinWid ! Desired line width
      Parameter (LinWid = 80)
C
C LOCAL VARIABLES
      Character*50 Date, Fmt
      Integer DLen, MLen, Spc
C
C     Note: Util_Date returns a newline at the end of the date string
C     (normal C string terminator), so we strip it off here.
C
      If (GA_NodeID() .eq. 0) then
         Date = ' '
         Call Util_Date( Date )
         DLen = Inp_StrLen(Date)-1
         MLen = Inp_StrLen(Msg)
C
C        Format things nicely
C        1) Message flush left, date flush right, if fits LinWid
C        2) Regular (long line) format as fall-back
C
         Spc = LinWid - (MLen + DLen + 2) 
         If ( Spc .ge. 1 ) then
            Write (Fmt, '(A, I3, A)') '(1X, A, ', Spc, 'X, ''@'', A)'
            Write (LUN, Fmt)  Msg(:MLen),Date(:DLen)
         Else
            Write (LUN, 9010) Msg(:MLen),Date(:DLen)
         EndIf
         Call Util_Flush(LUN)
      EndIf
C
 9010 Format(1X, A, ' @', A)
C
      Return
      End
