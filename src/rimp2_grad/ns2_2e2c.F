C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_NS2_2e2c -- 2c2e derivative integrals for non-sep 2-c
C     contribution to gradient
C
C REVISION./rimp2_grad/ns2_2e2c.F
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_NS2_2e2c(FitBas, Buf, LBuf, Scr, LScr, g_V,
     $     Grad, NCent)
      Implicit NONE
#include "errquit.fh"
      Integer FitBas            ![in]
      Integer LBuf              ![in]
      Double Precision Buf(LBuf) ![in]
      Integer LScr              ![in]
      Double precision Scr(LScr) ![in]
      Integer g_V               ![in]
      Integer NCent             ![in]
      Double precision Grad(3, NCent) ![inout]
C
C DESCRIPTION
C PRINT CONTROLS
C EXTERNAL ROUTINES
C INCLUDE FILES
#include "stdio.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
C LOCAL VARIABLES
      Integer Me, First, Last
      Integer ILo, IHi, JLo, JHi, NI, NJ
      Integer ILoCn, IHiCn, JLoCn, JHiCn, ICn, JCn
      Integer H_Patch, I_Patch
      Integer Atoms(2), A
      Double precision Grad_Contrib(3, 2) ! xyz for each atom
C
C STATEMENT FUNCTIONS
C
      Me = GA_NodeID()
C
C     Identify our patch of g_V
C
      Call GA_Distribution(g_V, Me, ILo, IHi, JLo, JHi)
C
C     Figure out which contractions the patch correponds to
C
      If ( .NOT. Bas_BF2Cn(FitBas, ILo, ILoCn) ) Call ErrQuit(
     $     'RIMP2_NS2_2e2c: bas_bf2cn failed ', ILo, BASIS_ERR)
      If ( .NOT. Bas_BF2Cn(FitBas, IHi, IHiCn) ) Call ErrQuit(
     $     'RIMP2_NS2_2e2c: bas_bf2cn failed ', IHi, BASIS_ERR)
      If ( .NOT. Bas_BF2Cn(FitBas, JLo, JLoCn) ) Call ErrQuit(
     $     'RIMP2_NS2_2e2c: bas_bf2cn failed ', JLo, BASIS_ERR)
      If ( .NOT. Bas_BF2Cn(FitBas, JHi, JHiCn) ) Call ErrQuit(
     $     'RIMP2_NS2_2e2c: bas_bf2cn failed ', JHi, BASIS_ERR)
C
      Write (LuOut, *) 'RIMP2_NS2_2e2c: initial patch ',
     $     Me, ILo, IHi, JLo, JHi, ILoCn, IHiCn, JLoCn, JHiCn
C
C     Adjust patch to match shell boundaries.
C
C     Left and top edges: If they don't start at the beginning of a
C     shell already, shift to start of NEXT contraction.
C
      If ( .NOT. Bas_Cn2BFR(FitBas, ILoCn, First, Last) ) Call
     $     ErrQuit( 'RIMP2_NS2_2e2c: bas_cn2bfr failed', ILoCn,
     &       BASIS_ERR)
      If ( ILo .ne. First) then
         ILoCn = ILoCn + 1
         ILo = First + 1
      EndIf
C
      If ( .NOT. Bas_Cn2BFR(FitBas, JLoCn, First, Last) ) Call
     $     ErrQuit( 'RIMP2_NS2_2e2c: bas_cn2bfr failed', JLoCn,
     &       BASIS_ERR)
      If ( JLo .ne. First) then
         JLoCn = JLoCn + 1
         JLo = First + 1
      EndIf
C
C     Right and bottom edges: simply extend to end of current contraction
C
      If ( .NOT. Bas_Cn2BFR(FitBas, IHiCn, First, Last) ) Call
     $     ErrQuit( 'RIMP2_NS2_2e2c: bas_cn2bfr failed', IHiCn,
     &       BASIS_ERR)
      If ( IHi .ne. Last) IHi = Last
C
      If ( .NOT. Bas_Cn2BFR(FitBas, JHiCn, First, Last) ) Call
     $     ErrQuit( 'RIMP2_NS2_2e2c: bas_cn2bfr failed', JHiCn,
     &       BASIS_ERR)
      If ( JHi .ne. Last) JHi = Last
C
      Write (LuOut, *) 'RIMP2_NS2_2e2c: adjusted patch ',
     $     Me, ILo, IHi, JLo, JHi, ILoCn, IHiCn, JLoCn, JHiCn
C
C     Allocate for, and get our patch.
C     If Buf(LBuf) is sufficient for up to 6 derivatives of a shell
C     block, it is certainly large enough to hold a patch of g_V
C     to contract them with.
C
      If ( .NOT. MA_Push_Get(MT_Dbl, LBuf, 'local g_V patch',
     $     H_Patch, I_Patch) ) Call ErrQuit(
     $     'RIMP2_NS2_2e2c: failed to allocate for local patch', LBuf,
     &       MA_ERR)
C
C     Loop over our shells
C
      Do ICn = ILoCn, IHiCn
         Do JCn = JLoCn, JHiCn
C
C           Compute integrals
C
            Call IntD_2e2c(FitBas, JCn, FitBas, ICn, LScr, Scr,
     $           LBuf, Buf, Atoms)
C
C           Check to see if we have integrals to deal with
C           (same-center derivatives are zero by translational
C           invariance.  The integral routine indicates this by
C           setting Atoms to something negative)
C
            If (Atoms(1) .gt. 0 .AND. Atoms(2) .gt. 0) then
C
C              Figure out the basis function ranges.
C
               If ( .NOT. Bas_Cn2BFR(FitBas, ICn, ILo, IHi) ) Call
     $              ErrQuit( 'RIMP2_NS2_2e2c: bas_cn2bfr failed', ICn,
     &       BASIS_ERR)
               If ( .NOT. Bas_Cn2BFR(FitBas, JCn, JLo, JHi) ) Call
     $              ErrQuit( 'RIMP2_NS2_2e2c: bas_cn2bfr failed', JCn,
     &       BASIS_ERR)
C
               NI = IHi - ILo + 1
               NJ = JHi - JLo + 1
C
C              Get corresponding patch
C
               Call GA_Get(g_V, ILo, IHi, JLo, JHi, Dbl_MB(I_Patch), NI)
c$$$               Write (LuOut, *) 'RIMP2_NS2_2e2c: ICn, JCn ', ICn, JCn,
c$$$     $              ILo, IHi, JLo, JHi, NI, NJ
c$$$               Write (LuOut, *) 'RIMP2_NS2_2e2c: Atoms 0', Atoms
C
C              Contract patch with integrals
C
               Call DGEMV('T', NI*NJ, 6, 1.0d0, Buf, NI*NJ,
     $              Dbl_MB(I_Patch), 1, 0.0d0, Grad_Contrib, 1)
c$$$               Write (LuOut, *) 'RIMP2_NS2_2e2c: Grad_Contrib ',
c$$$     $              Grad_Contrib
c$$$               Write (LuOut, *) 'RIMP2_NS2_2e2c: Atoms 1', Atoms
C
C              Put the gradient contributions where they belong
C
               Do A = 1, 2
c$$$               Write (LuOut, *) 'RIMP2_NS2_2e2c: Atoms 2', Atoms
c$$$                  Write (LuOut, *) 'Copying grad_contrib(1, ', A,')'
c$$$                  Write (LuOut, *) 'to grad(1, ', Atoms(a), ')'
                  Call DAXPY(3, 1.0d0, Grad_Contrib(1, A), 1,
     $                 Grad(1, Atoms(A)), 1)
               EndDo ! A
c$$$               Write (LuOut, *) 'RIMP2_NS2_2e2c: Grad ', Grad
C
            EndIf               ! Integrals to process
         EndDo                  ! JCn
      EndDo                     ! ICn
C
      If ( .NOT. MA_Pop_Stack(H_Patch) ) Call ErrQuit(
     $     'RIMP2_NS2_2e2c: failed to free patch', 0, MA_ERR)
C
C     Synchronize and combine gradient contributions on each node
C
      Call GA_Sync
      Call GA_DGOp(9234, Grad, 3*NCent, '+')
C
      Return
      End
