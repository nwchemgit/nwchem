      Logical function Close_3ci_file( LU, Stat)
C$Id: close_3ci.F,v 1.3 1998-01-31 20:14:13 bernhold Exp $
      Implicit NONE
      Integer LU
      Character*1 Stat
C
      Logical LSame
      External LSame
C
#include "stdio.fh"
#ifdef RIMP2_USE_DA
#include "dra.fh"
#endif
C
      Logical Keep
      Integer Code
C
      Keep = LSame( Stat, 'K')
      If ( .NOT. ( Keep .OR. LSame( Stat, 'D'))) then
         Write (LuErr, *) 'Close_3cI_File: Bad status, must be ',
     $      '''KEEP'' or ''DELETE''.'
         Close_3ci_File = .FALSE.
         Return
      EndIf
C
#ifdef RIMP2_USE_DA
      If ( Keep ) then
         If ( DRA_Close(LU) .ne. 0) Goto 8000
      Else
         If ( DRA_Delete(LU) .ne. 0) Goto 8000
      EndIf
#else
      If ( Keep ) then
         Close(LU, STATUS='KEEP', ERR=8000, IOSTAT=Code)
      Else
         Close(LU, STATUS='DELETE', ERR=8000, IOSTAT=Code)
      EndIf
#endif
C
      Close_3ci_File = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Close_3cI_File: Error closing file: ', Code
      Call ERRquit('Close_3cI_File',code)
C
      Close_3ci_File = .FALSE.
C
      Return
      End
