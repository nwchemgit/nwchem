      Logical Function Put_3cInts(LU, G_I, NA, ILo, IHi)
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
      Integer Junk, Dim1, Dim2
c$$$      Double precision ChkSum, MxElem
c$$$C
c$$$      External GA_DASum_Patch, GA_DAMax_Patch
c$$$      Double precision GA_DASum_Patch, GA_DAMax_Patch
C
      Call GA_Distribution( G_I, GA_NodeID(), AIlo, AIhi, Slo, Shi)
      NS = SHi-SLo+1
C
c$$$      Call GA_Inquire(G_I, Junk, Dim1, Dim2)
c$$$      Do I = ILo, IHi
c$$$         LocalI = I - ILo + 1
c$$$         ChkSum = GA_DASum_Patch( G_I, 1+NA*(LocalI-1), NA*LocalI,
c$$$     $      1, Dim2)
c$$$         If ( GA_NodeID() .eq. 0)
c$$$     $      Write (LuOut, 9000) I, LocalI, ChkSum
c$$$         MxElem = GA_DAMax_Patch( G_I, 1+NA*(LocalI-1), NA*LocalI,
c$$$     $      1, Dim2)
c$$$         If ( GA_NodeID() .eq. 0)
c$$$     $      Write (LuOut, 9010) I, LocalI, MxElem
c$$$      EndDo
c$$$ 9000 Format(1X, 'Put_3cInts: Writing record ', I6, ' from GA block ',
c$$$     $   I6, ' with checksum ', F30.15)
c$$$ 9010 Format(1X, 'Put_3cInts: Writing record ', I6, ' from GA block ',
c$$$     $   I6, ' with maximum ', F30.15)
C
      If (NS*NA .gt. 0) then
C
C        Since we need a contiguous block of memory to write to the
C        file, the simplest thing is to create a local array to hold it.
C
         Status = MA_Push_Get(MT_Dbl, NA*NS, 'Put_3cInts scratch',
     $      H_Local, Local_Idx)
         If (.NOT. Status) then
            Write (LuErr, *) 'Put_3cInts: Unable to allocate scratch ',
     $         NA*NS, ' space.'
            Put_3cInts = .FALSE.
            Return
         EndIf
C
         Do I = ILo, IHi
            LocalI = I - ILo + 1
            Call GA_Get( G_I, 1+NA*(LocalI-1), NA*LocalI, SLo, SHi,
     $         Dbl_MB( Local_Idx), NA)
            Write(LU, REC=I, Err=8000) (Dbl_MB( Local_Idx+Junk),
     $         Junk=0,NA*NS-1)
         EndDo
C
         Status = MA_Pop_Stack( H_Local)
         If (.NOT. Status) then
            Write (LuErr, *)
     $         'Put_3cInts: Problems freeing scratch space.'
            Put_3cInts = .FALSE.
            Return
         EndIf
      EndIf
C
      Put_3cInts = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Put_3cInts: Error writing ', ILo, I, IHi
      Call Errquit('Put_3cInts',0)
C
      Put_3cInts = .FALSE.
C
      Return
      End
