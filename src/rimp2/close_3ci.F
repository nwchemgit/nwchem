      Logical function Close_3ci_file( LU, Stat,actualname)
C$Id$
      Implicit NONE
#include "errquit.fh"
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
#include "inp.fh"
#ifdef RIMP2_USE_DA
#include "dra.fh"
#endif
C
      Logical Keep
      Integer Code
      Integer ierr
      character*80 errmsg      
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
      if (LU.eq.-1) then
c
c     this means the file was not opened since recl.eq.0
c
      else
         ierr=eaf_close(LU)
         if (ierr .ne. 0) then
            write(luout,*) ga_nodeid(),' closing FD =',LU,
     ,           actualname(1:inp_strlen(actualname)) 
            call eaf_errmsg(ierr, errmsg)
            write(LuOut,*) ' IO error message ',
     $           errmsg(1:inp_strlen(errmsg))
            call util_flush(luout)
            call errquit('close_3ci: eaf_close failed',0, DISK_ERR)
         endif
         If ( Keep ) then
         else
            if (eaf_delete(actualname(1:inp_strlen(actualname)))
     .           .ne. 0) call errquit
     $           ('close3ci: delete failed', 0, DISK_ERR)
         EndIf
      endif
#endif
C
      Close_3ci_File = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Close_3cI_File: Error closing file: ', Code
      Call ERRquit('Close_3cI_File',code,0)
C
      Close_3ci_File = .FALSE.
C
      Return
      End
