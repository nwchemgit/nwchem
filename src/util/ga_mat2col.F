      Subroutine GA_Mat2Col( G_A, AILo, AIHi, AJLo, AJHi,
     $   G_B, BILo, BIHi, BJLo, BJHi)
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer G_A, G_B, AILo, AIHi, AJLo, AJHi, BILo, BIHi, BJLo, BJHi
#include "mafdecls.fh"
#include "global.fh"
C LOCAL VARIABLES
      Integer LenCol, Me, NProc, BOff, Col
      Integer H_Col, I_Col
C
      Call GA_Sync
      Me = GA_NodeID()
      NProc = GA_NNodes()
C
C     Assume handles are okay, types is dbl, and dimensions are sensible
C
C     Use local memory long enough to hold a colum
C
      LenCol = AIHi - AILo + 1
      If ( .NOT. MA_Push_Get(MT_Dbl, LenCol, 'ga_mat2col tmp', H_Col,
     $   I_Col) ) Call ErrQuit(
     $   'GA_Mat2Col: Unable to allocate local memory', LenCol, MA_ERR)
C
C     Simple static load balancing, ignore data distribution.
C     Get a column and put it into the target.
C
      BOff = BILo
      Do Col = AJLo, AJHi
         If ( Mod(Col, NProc) .eq. Me) then
            Call GA_Get(G_A, AILo, AIHi, Col, Col,
     $         DBL_MB( I_Col), LenCol)
C
            Call GA_Put(G_B, BOff, BOff+LenCol-1, BJLo, BJHi,
     $         DBL_MB( I_Col), LenCol)
         EndIf
         BOff = BOff + LenCol
      EndDo
C
      If ( .NOT. MA_Pop_Stack( H_Col )) Call ErrQuit(
     $   'GA_Mat2Col: Unable to free local memory', 0, MA_ERR)
C
      Call GA_Sync
C
      Return
      End
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      Subroutine GA_Col2Mat( G_A, AILo, AIHi, AJLo, AJHi,
     $   G_B, BILo, BIHi, BJLo, BJHi)
      Implicit NONE
#include "errquit.fh"
      Integer G_A, G_B, AILo, AIHi, AJLo, AJHi, BILo, BIHi, BJLo, BJHi
#include "mafdecls.fh"
#include "global.fh"
C LOCAL VARIABLES
      Integer LenCol, Me, NProc, AOff, Col
      Integer H_Col, I_Col
C
      Call GA_Sync
      Me = GA_NodeID()
      NProc = GA_NNodes()
C
C     Assume handles are okay, types is dbl, and dimensions are sensible
C
C     Use local memory long enough to hold a colum
C
      LenCol = BIHi - BILo + 1
      If ( .NOT. MA_Push_Get(MT_Dbl, LenCol, 'ga_col2mat tmp', H_Col,
     $   I_Col) ) Call ErrQuit(
     $   'GA_Col2Mat: Unable to allocate local memory', LenCol, MA_ERR)
C
C     Simple static load balancing, ignore data distribution.
C     Get a column and put it into the target.
C
      AOff = AILo
      Do Col = BJLo, BJHi
         If ( Mod(Col, NProc) .eq. Me) then
            Call GA_Get(G_A, AOff, AOff+LenCol-1, AJLo, AJHi,
     $         DBL_MB( I_Col), LenCol)
C
            Call GA_Put(G_B, BILo, BIHi, Col, Col,
     $         DBL_MB( I_Col), LenCol)
         EndIf
         AOff = AOff + LenCol
      EndDo
C
      If ( .NOT. MA_Pop_Stack( H_Col )) Call ErrQuit(
     $   'GA_Col2Mat: Unable to free local memory', 0, MA_ERR)
C
      Call GA_Sync
C
      Return
      End
