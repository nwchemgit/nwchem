      Logical function Close_3ci_file( LU, Stat,actualname)
C$Id: close_3ci.F,v 1.4 2003-02-22 05:58:52 edo Exp $
      Implicit NONE
      Integer LU
      Character*1 Stat
      character*(*) actualname
C
      Logical LSame
      External LSame
C
#include "stdio.fh"
#include "eaf.fh"
#include "global.fh"
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
         if (eaf_close(LU) .ne. 0) call errquit
     $        ('close_3ci: eaf_close failed',0)
cold         Close(LU, STATUS='KEEP', ERR=8000, IOSTAT=Code)
      Else
         if (eaf_close(LU) .ne. 0) call errquit
     $        ('close_3ci: eaf_close failed',0)
cedo need to transfer actualname
         if (eaf_delete(actualname) .ne. 0) call errquit
     $        ('close3ci: delete failed', 0)
cold         Close(LU, STATUS='DELETE', ERR=8000, IOSTAT=Code)
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
