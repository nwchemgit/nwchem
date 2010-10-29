      Logical Function Get_3cInts(LU, G_I, NA, ILo, IHi)
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer LU, G_I, NA, ILo, IHi
C
#ifdef RIMP2_USE_DA
C
#include "dra.fh"
#include "global.fh"
C
C LOCAL VARIABLES
      Integer Type, AILen, FitLen, Req
C
      Call GA_Inquire(G_I, Type, AILen, FitLen)
C
C     Note:  We don't necessarily want to read in the whole GA.
C
      If ( DRA_Read_Section(.FALSE.,
     $   G_I, 1,            NA*(IHi-ILo+1), 1, FitLen,
     $   LU,  1+NA*(ILo-1), NA*IHi,         1, FitLen, Req) .ne. 0)
     $   Call ErrQuit( 'Get_3cInts: Error from DRA_Read_Section', ILo,
     &       UNKNOWN_ERR)
C
      If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $   'Get_3cInts: Error from DRA_Wait', ILo, UNKNOWN_ERR)
C
      Get_3cInts = .TRUE.
C
#else
C
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "eaf.fh"
#include "inp.fh"
C
C LOCAL VARIABLES
      Logical Status
      Integer AILo, AIHi, SLo, SHi, NS, H_Local, Local_Idx, I, LocalI
      Integer ierr
      character*80 errmsg      
      double precision offset
C
      Call GA_Distribution( G_I, GA_NodeID(), AIlo, AIhi, Slo, Shi)
      NS = SHi-SLo+1
C
      If (NS*NA .gt. 0) then
C
C        Since we need a contiguous block of memory to write to the
C        file, the simplest thing is to create a local array to hold it.
C
         Status = MA_Push_Get(MT_Dbl, NA*NS, 'Get_3cInts scratch',
     $      H_Local, Local_Idx)
         If (.NOT. Status) then
            Write (LuErr, *) 'Get_3cInts: Unable to allocate scratch ',
     $         NA*NS, ' space.'
            Get_3cInts = .FALSE.
            Return
         EndIf
C
         Do I = ILo, IHi
            LocalI = I - ILo + 1
               offset=8.d0*NA*NS*(I-1)
               ierr = eaf_read(LU, offset, dbl_mb(local_idx),
     .              8*NA*NS)
               if (ierr .ne. 0) then
                  call eaf_errmsg(ierr, errmsg)
                  write(LuOut,*) ' IO offset ', offset
                  write(LuOut,*) ' IO error message ',
     $                 errmsg(1:inp_strlen(errmsg))
                  call errquit('put_3cints: read failed',0, DISK_ERR)
               endif
               offset=offset+8.d0*NA*NS

            Call GA_Put( G_I, 1+NA*(LocalI-1), NA*LocalI, SLo, SHi,
     $         Dbl_MB( Local_Idx), NA)
         EndDo
C
         Status = MA_Pop_Stack( H_Local)
         If (.NOT. Status) then
            Write (LuErr, *)
     $         'Get_3cInts: Problems freeing scratch space.'
            Get_3cInts = .FALSE.
            Return
         EndIf
      EndIf
C
      Get_3cInts = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Get_3cInts: Error writing ', ILo, I, IHi
      Call Errquit('Get_3cInts',0, UNKNOWN_ERR)
C
      Get_3cInts = .FALSE.
C
#endif
C
      Return
      End
