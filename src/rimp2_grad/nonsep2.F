C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_NonSep2 -- Non-separable gradient contrib from 2c integrals
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_NonSep2(BraKetSame, TopSpin, NFrzO, NAct, NVir,
     $     FitBas, D_Int, D_Gam, Grad, NCent )
      Implicit NONE
#include "errquit.fh"
      Logical BraKetSame        ![in]
      Integer TopSpin           ![in]
      Integer NFrzO(TopSpin)    ![in]
      Integer NAct(TopSpin)     ![in]
      Integer NVir(TopSpin)     ![in]
      Integer FitBas            ![in]
      Integer D_Int(TopSpin)    ![in]
      Integer D_Gam(TopSpin)    ![in]
      Integer NCent             ![in]
      Double Precision Grad(3, NCent) ![inout]
C
C DESCRIPTION
C
C     2 (ia|D) V-1(D,D') V'(D',T) G(T,ia)
C
C NOTES
C     Equations have been derived only for the 'V' approximation
C     Assumes that the ket integrals are actually (ia|D) V-1(D,D')
C
C PRINT CONTROLS
C EXTERNAL ROUTINES
c$$$      Logical  GA_Create
c$$$      External GA_Create
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "bas.fh"
#include "util.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk
      Parameter (MinChunk = -64)
C
C LOCAL VARIABLES
      Integer Me
      Integer NFit, Spin
      Integer g_V, H_Scr, I_Scr, H_Buf, I_Buf, LBuf, LScr
      Integer NTotal, IBlock, GBlock
      Double Precision Initial, Final, Block
C
C STATEMENT FUNCTIONS
C
C
      Me = GA_NodeID()
C
      If ( BraKetSame ) then
         If ( Me .eq. 0) Write (LuOut, 9000)
         Call ErrQuit(
     $     'RIMP2_NonSep2: invocation inconsistent with implementation',
     $        0, INPUT_ERR)
      EndIf
 9000 Format(/1X, 'ERROR: RIMP2_NonSep2 implementation requires an ',
     $            'asymmetric integral',
     $       /1X, 'ERROR: approximation with the V^{-1} subsumed ',
     $            'into the ket fitting',
     $       /1X, 'ERROR: basis index (BraKetSame must be .FALSE.).')
C
      If ( .NOT. Bas_NumBF( FitBas, NFit) ) Call ErrQuit(
     $     'RIMP2_NonSep2: can''t get fitting basis size', FitBas,
     &       BASIS_ERR)
C
C     *******************************************************
C     * Create target GA for 3c integral-gamma3 contraction *
C     *******************************************************
C     NOTE: Should probably distribute this on shell boundaries
C
      If ( .NOT. GA_Create(MT_Dbl, NFit, NFit, 'V-like temp',
     $     MinChunk, MinChunk, g_V) ) Call ErrQuit(
     $     'RIMP2_NonSep2: can''t allocate V-like temporary', 0, GA_ERR)
      Call GA_Zero(g_V)
C
C     ************************************
C     * Contract 3c integrals and Gamma3 *
C     ************************************
C
      Do Spin = 1, TopSpin
         Call RIMP2_NS2_I3G3(NFrzO(Spin), NAct(Spin), NVir(Spin),
     $        NFit, D_Int(Spin), D_Gam(Spin), g_V)
      EndDo
c$$$      Call GA_Print( g_V )
C
C     **************************************
C     * Contract with derivative integrals *
C     **************************************
C
      Call Int_Mem_2e2c(LBuf, LScr)
      LBuf = 2*LBuf                                      !!! CHEATING!!!
      If (.NOT. MA_Push_Get(MT_DBL, LBuf, 'deriv int buf',
     $     H_Buf, I_Buf ) ) Call ErrQuit(
     $     'RIMP2_NonSep2: ma failed (buf)', LBuf, MA_ERR)
      If (.NOT. MA_Push_Get(MT_DBL, LScr, 'deriv int scr',
     $     H_Scr, I_Scr ) ) Call ErrQuit(
     $     'RIMP2_NonSep2: ma failed (scr)', LScr , MA_ERR)
C
      Call RIMP2_NS2_2e2c(FitBas, Dbl_MB(I_Buf), LBuf,
     $     Dbl_MB(I_Scr), LScr, g_V, Grad, NCent)
C
      If (.NOT. MA_Pop_Stack( H_Scr) ) Call ErrQuit(
     $     'RIMP2_NonSep2: pop of scr failed', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack( H_Buf) ) Call ErrQuit(
     $     'RIMP2_NonSep2: pop of buf failed', 0, MA_ERR)
C
C     Free remaining GAs
C
      If ( .NOT. GA_Destroy( g_V ) ) Call ErrQuit(
     $     'RIMP2_NonSep2: can''t destroy V-like GA', 0, GA_ERR)
C
C     Print partial results if requested
C
      Call GA_Sync
      If ( Util_Print('partial grad', Print_Debug) .AND. Me .eq. 0)
     $     Call DGeWr(3, NCent, Grad, 3, 0, 0, LuOut,
     $     'Partial Non-Separable Gradient (2c contributions)',
     $     80, 'E12.6')
C
      Return
      End
