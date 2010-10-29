C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_P2ab -- Form P^{(2)}_{ab} for batch in memory
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_P2AB(SameSpin, IJSym, IJTri, ILo, IHi, NVirI,
     $     JLo, JHi, NVirJ, Mixed, G_T, G_P2, Vir1)
      Implicit NONE
#include "errquit.fh"
C
      Logical SameSpin ![in]
      Logical IJSym
      Logical IJTri
      Integer ILo
      Integer IHi
      Integer NVirI
      Integer JLo
      Integer JHi
      Integer NVirJ
      Logical Mixed
      Integer G_T
      Integer G_P2
      Integer Vir1
C
C DESCRIPTION
C     Evaluate P^{(2)}_{ab} = 1/2 Sum(ijc) t(ij,ac) t(ij,bc)
C
C     Non-zero spin cases:
C     a b   i j c
C     a a   a a a   t(IJ,AC) t(IJ,BC)

C     a a   a b b   t(Ij,Ac) t(Ij,Bc)
C     a a   b a b   t(iJ,Ac) t(iJ,Bc) --> t(Ji,Ac) t(Ji,Bc)
C
C     b b   b a a   t(iJ,aC) t(iJ,bC) --> t(Ji,Ca) t(Ji,Cb)
C     b b   a b a   t(Ij,aC) t(Ij,bC) --> t(Ij,Ca) t(Ij,Cb)
C
C     b b   b b b   t(ij,ac) t(ij,bc)
C
C     So we learn
C     1) Pure spin amplitudes contribute to matching P spin
C     2) Mixed-spin amplitudes contribute to both P spins
C     3) Mixed-spin contributions carry factor of 2
C
C EXTERNAL ROUTINES
C INCLUDED FILES
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Me, Owner
      Integer H_Block, I_Block, H_P2, I_P2
      Integer Type, LT1, LT2, NI, NJ, I, J, JTop
      Double precision Fact
C
C
      Me = GA_NodeID()
C
C     We need a local matrix to hold the amplitude block and one in which
C     to accumulate the contribution.
C
      If ( .NOT. MA_Push_Get( MT_Dbl, NVirI*NVirJ,
     $   'Local temp for RIMP2_P2ab', H_Block, I_Block) ) Call ErrQuit(
     $   'RIMP2_P2ab: can''t allocate local temp', NVirI*NVirJ, MA_ERR)
      If ( Mixed ) then
         If ( .NOT. MA_Push_Get( MT_Dbl, NVirJ*NVirJ,
     $      'Local temp for RIMP2_P2ab', H_P2, I_P2) ) Call ErrQuit(
     $      'RIMP2_P2ab: can''t allocate local temp', NVirJ*NVirJ,
     &       MA_ERR)
         Call DFill(NVirJ*NVirJ, 0.0d0, Dbl_MB(I_P2), 1)
      Else
         If ( .NOT. MA_Push_Get( MT_Dbl, NVirI*NVirI,
     $      'Local temp for RIMP2_P2ab', H_P2, I_P2) ) Call ErrQuit(
     $      'RIMP2_P2ab: can''t allocate local temp', NVirI*NVirI,
     &       MA_ERR)
         Call DFill(NVirI*NVirI, 0.0d0, Dbl_MB(I_P2), 1)
      EndIf
C
C
C     To divide up the work, we will do every block of the matrix for
C     which we own the upperleft corner.
C
      Call GA_Inquire( G_T, Type, LT1, LT2)
      NI = LT1 / NVirI
      NJ = LT2 / NVirJ
C
      Call GA_Sync
      Do I = 1, NI
         If (IJSym .AND. ILo .eq. JLo .AND. IHi. eq. JHi
     $      .AND. .not. IJTri) then
            JTop = I
         Else
            JTop = NJ
         EndIf
         Do J = 1, JTop
            If ( .NOT. GA_Locate(G_T, (I-1)*NVirI+1, (J-1)*NVirJ+1,
     $         Owner) ) Call ErrQuit(
     $         'RIMP2_P2ab: indexing problem with t amplitudes', 0,
     &       GA_ERR)
C
            If ( Owner .eq. Me) then
C
C              Get a local copy of the full block (can't use ga_access
C              since it is not guaranteed to live entirely on this node.
C
               Call GA_Get(G_T, (I-1)*NVirI+1, I*NVirI, (J-1)*NVirJ+1,
     $            J*NVirJ, Dbl_MB(I_Block), NVirI)
C
C              P^{(2)}_{ab} = t^{ac}_{ij} t^{bc}_{ij}
C
               If ( Mixed ) then
                  Call DGEMM('T', 'N', NVirJ, NVirJ, NVirI, 1.0d0,
     $               Dbl_MB(I_Block), NVirI, Dbl_MB(I_Block), NVirI,
     $               1.0d0, Dbl_MB(I_P2), NVirJ)
               Else
                  Call DGEMM('N', 'T', NVirI, NVirI, NVirJ, 1.0d0,
     $               Dbl_MB(I_Block), NVirI, Dbl_MB(I_Block), NVirI,
     $               1.0d0, Dbl_MB(I_P2), NVirI)
               EndIf
            EndIf
         EndDo
      EndDo
C
C     Accumulate our contribution into global P(2) array
C     To understand the factor, see the discussion of spin cases, above.
C
      Fact = 1.0d0
      If ( SameSpin) Fact = 0.5d0
c$$$      If ( .NOT. IJSym ) Fact = 1.0d0
C
      If ( Mixed ) then
         Call GA_Acc(G_P2, Vir1, Vir1+NVirJ-1, Vir1, Vir1+NVirJ-1,
     $      Dbl_MB(I_P2), NVirJ, Fact)
      Else
         Call GA_Acc(G_P2, Vir1, Vir1+NVirI-1, Vir1, Vir1+NVirI-1,
     $      Dbl_MB(I_P2), NVirI, Fact)
      EndIf
C
C     Free local memory
C
      If ( .NOT. MA_Pop_Stack(H_P2) ) Call ErrQuit(
     $   'RIMP2_P2ab: can''t free local temp', H_P2, MA_ERR)
      If ( .NOT. MA_Pop_Stack(H_Block) ) Call ErrQuit(
     $   'RIMP2_P2ab: can''t free local temp', H_Block, MA_ERR)
C        
      If (Util_Print('partial p2', Print_Debug) ) then
         Write (LuOut, *) 'After RIMP2_P2ab...'
         Call GA_Print(G_P2)
      EndIf
C
      Return
      End

