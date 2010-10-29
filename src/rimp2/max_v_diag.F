      Subroutine Max_V_Diagonal(Basis, Big)
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer Basis
      Double precision Big
C
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
C
      Logical GA_Create_One_By_Shells
      External GA_Create_One_By_Shells
C
      Logical Status
      Integer G_T, N, I, MySLo, MySHi, SShLo, SShHi, BFLo, BFHi, S
      Integer Jnk1, Jnk2
      Integer LInt, H_Int, Int_Idx, LScr, H_Scr, Scr_Idx
C
      Big = FP_Zero
C
      Status = .TRUE.
      Status = Status .AND. GA_Create_One_By_Shells( MT_Dbl, 1, Basis, '
     $   Template', G_T)
      Call GA_Distribution(G_T, GA_NodeID(), Jnk1, Jnk2, MySLo, MySHi)
      Status = Status .AND. GA_Destroy(G_T)
C
C     GA_Distribution will return nonsensical ranges if there is no data
C     on this node (0:-1).  They are designed to keep a Fortran DO loop
c     from going, not to be fed into basis function lookup routines!
C
      If ( MySLo .gt. 0 .AND. MySHi .gt. 0) then
         Status = Status .AND. Bas_BF2Cn( Basis, MySLo, SShlo)
         Status = Status .AND. Bas_BF2Cn( Basis, MySHi, SShhi)
      Else
         SShLo =  0
         SShHi = -1
         Return
      EndIf
C
      Call Int_Mem_2e2c(LInt, LScr)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LInt, 'Int buffer',
     $   H_Int, Int_Idx)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LScr, 'Int scratch',
     $   H_Scr, Scr_Idx)
C
      If ( .NOT. Status )
     $   Call ErrQuit('Max_V_Diagonal: Big problems!',0, UNKNOWN_ERR)
C
      Do S = SShLo, SShHi
         Status = Status .AND. Bas_Cn2BFR(Basis, S, BFLo, BFHi)
         N = BFHi - BFLo + 1
C
         Call Int_2e2c (Basis, S, Basis, S,
     $      LScr, dbl_mb(Scr_Idx), LInt, dbl_mb(Int_Idx) )
C
C        Since these guys had better be positive anyway, the absolute
C        value aspect is irrelevant.  And we only need to scan the
C        diagonals
C
         Do I = 1, N
            Big = Max( Big, Dbl_MB( Int_Idx + (I-1) * (N+1)) )
         EndDo
      EndDo
C
      Status = Status .AND. MA_Pop_Stack( H_Scr )
      Status = Status .AND. MA_Pop_Stack( H_Int )
C
      Return
      End



