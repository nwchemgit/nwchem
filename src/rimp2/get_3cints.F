      Logical Function Get_3cInts(LU, G_I, NA, ILo, IHi)
      Implicit NONE
      Integer LU, G_I, NA, ILo, IHi
C
#include "global.fh"
#include "mafdecls.h"
#include "stdio.fh"
C
C LOCAL VARIABLES
      Logical Status
      Integer AILo, AIHi, SLo, SHi, NS, H_Local, Local_Idx, I, LocalI
      Integer Junk
C
      Call GA_Distribution( G_I, GA_NodeID(), AIlo, AIhi, Slo, Shi)
      NS = SHi-SLo+1
C
C     Since we need a contiguous block of memory to write to the
C     file, the simplest thing is to create a local array to hold it.
C
      Status = MA_Push_Get(MT_Dbl, NA*NS, 'Get_3cInts scratch',
     $   H_Local, Local_Idx)
      If (.NOT. Status) then
         Write (LuErr, *) 'Get_3cInts: Unable to allocate scratch ',
     $      NA*NS, ' space.'
         Get_3cInts = .FALSE.
         Return
      EndIf
C
      Do I = ILo, IHi
         LocalI = I - ILo + 1
         Read(LU, REC=I, Err=8000) (Dbl_MB( Local_Idx+Junk),
     $      Junk=0,NA*NS-1)
c$$$         Call DGeWr(NA, NS, Dbl_MB(Local_Idx), NA, 0, SLo-1, LuOut,
c$$$     $      'Local integrals from disk', 80, 'E15.6', Junk)
         Call GA_Put( G_I, 1+NA*(LocalI-1), NA*LocalI, SLo, SHi,
     $      Dbl_MB( Local_Idx), NA)
      EndDo
C
      Status = MA_Pop_Stack( H_Local)
      If (.NOT. Status) then
         Write (LuErr, *) 'Get_3cInts: Problems freeing scratch space.'
         Get_3cInts = .FALSE.
         Return
      EndIf
C
      Get_3cInts = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Get_3cInts: Error writing ', ILo, I, IHi
      Call Errquit('Get_3cInts',0)
C
      Get_3cInts = .FALSE.
C
      Return
      End
