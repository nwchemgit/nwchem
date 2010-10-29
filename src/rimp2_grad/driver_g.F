C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_G -- Master routine for RI-MP2 gradient evaluation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_Driver_G(RTDB, Geom, AOBas, NAO,
     $     MaxSpin, TopSpin, NMO, NFrzO, NAct, NVir, FitBas,
     $     BraKetSame, g_P2, g_W2, Eig, LDEig, g_AOXF,
     $     Tol2e_Fock, JobTitle, D_Int, D_Gam, NCent)
      Implicit NONE
#include "errquit.fh"
C
      Integer RTDB          ![in]
      Integer Geom          ![in]
      Integer AOBas         ![in]
      Integer NAO           ![in]
      Integer MaxSpin       ![in]
      Integer TopSpin       ![in]
      Integer NMO(TopSpin)  ![in]
      Integer NFrzO(TopSpin)![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer FitBas        ![in]
      Logical BraKetSame    ![in]
      Integer g_P2(TopSpin) ![in]
      Integer g_W2(TopSpin) ![in]
      Integer LDEig ![in]
      Double Precision Eig(LDEig, TopSpin) ![in]
      Integer g_AOXF(TopSpin) ![in]
      Double precision Tol2e_Fock ![in]
      Character*(*) JobTitle ![in]
      Integer D_Int(MaxSpin, 2, 2) ![in]
      Integer D_Gam(MaxSpin) ![in]
      Integer NCent ![in]
C
C
C DESCRIPTION
C
C NOTES
C EXTERNAL ROUTINES
      Logical  File_Write_GA, File_Read_GA, CPHF
      External File_Write_GA, File_Read_GA, CPHF
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
#include "dra.fh"
#include "bas.fh"
#include "rtdb.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MyMaxSpin
      Parameter (MyMaxSpin = 2)
C
C LOCAL VARIABLES
      Integer ISpin, NFit, Length, C, A, V
      Integer g_L(MyMaxspin), g_RHS, H_Grad, I_Grad
      Logical PrInfo, PrPrgRpt
      Character*256 String1, FN_CPHF_RHS, FN_CPHF_Sol
      Integer Bases(2)
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
C
C     ************************
C     * Initialization stuff *
C     ************************
C
      If ( TopSpin .gt. MaxSpin) Call ErrQuit(
     $     'RIMP2_Driver_G: fatal program error: TopSpin > Maxspin',
     $     MaxSpin, INPUT_ERR)
      If ( MaxSpin .ne. MyMaxSpin ) Call ErrQuit(
     $     'RIMP2_Driver_G: fatal program error: Maxspin != MyMaxSpin',
     $     MaxSpin, INPUT_ERR)
C
      If ( .NOT. Bas_NumBF( FitBas, NFit) ) Call ErrQuit(
     $     'RIMP2_Driver_G: can''t get fitting basis size', FitBas,
     &       BASIS_ERR)
C
      PrInfo =   Util_Print('information', Print_Low)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
      If ( (PrInfo .AND. GA_NodeID() .eq. 0) .or. PrPrgRpt ) then
         Write (LuOut, 9000) 
      EndIf
 9000 Format(/1X, 'Gradient evaluation')
C
C     Prepare gradient array.   This could be done as a GA instead
C     of being replicated, but then we would have nodes constantly
C     sending tiny messages to a tiny array -- seems undesirable!
C
      If ( .NOT. MA_Push_Get( MT_Dbl, 3*NCent, 'RI-MP2 gradient',
     $     H_Grad, I_Grad) ) Call ErrQuit(
     $     'RIMP2_Driver_G: Failed to allocate gradient', 3*NCent,
     &       MA_ERR)
C
      Call DFill(3*NCent, 0.0d0, Dbl_MB(I_Grad), 1)
C
C     ***************
C     * Complete W2 *
C     ***************
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Completing energy-weighted density')
      Call RIMP2_Mk_W2(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct, NVir,
     $     g_P2, g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
C
C     ***********************************
C     * Do non-separable gradient terms *
C     ***********************************
C     There are terms involving both the derivative two-center integrals
C     and derivative three-center integrals.
C
      If ( PrPrgRpt ) Call Progress_Report( LuOut,
     $     '    Non-separable gradient contrib (2-c deriv ints)')
C
C     Initialize for derivative integrals
C
      Bases(1) = AOBas
      Bases(2) = FitBas
      Call Int_Init(RTDB, 2, Bases)
      Call Int_Terminate()
      Call IntD_Init(RTDB, 2, Bases)
C
C     Needs ket integrals only, and gamma3
C
      Call RIMP2_NonSep2(BraKetSame, TopSpin, NFrzO, NAct, NVir, FitBas,
     $     D_Int(1, 2, 1), D_Gam, Dbl_MB(I_Grad), NCent)
C
      If ( PrPrgRpt ) Call Progress_Report( LuOut,
     $     '    Non-separable gradient contrib (3-c deriv ints)')
C
      Call IntD_Terminate()
C
C     Separable gradient contribution expect to see the non-separable
C     component on the RTDB
C
      If ( .NOT. RTDB_Put( RTDB, 'mp2:nonseparable gradient', MT_Dbl,
     $     3*NCent, Dbl_MB(I_Grad) ) ) Call ErrQuit(
     $     'RIMP2_Driver_G: failed to store gradient in RTDB', 0,
     &       RTDB_ERR)
C
C     *******************************
C     * Do separable gradient terms *
C     *******************************
C     Non-separable contribution must already be complete
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Separable gradient contributions')
C
      Call RIMP2_Prep_Sep(TopSpin, g_P2, g_W2, g_AOXF)
      Call Grad_Force( RTDB, AOBas, Geom)
      Call RIMP2_Clean_Sep(TopSpin)
C
C     Clean up
C
      If ( .NOT. MA_Pop_Stack(H_Grad) ) Call ErrQuit(
     $     'RIMP2_Driver_G: ma_pop_stack failed', H_Grad, MA_ERR)
C
      Return
      End
