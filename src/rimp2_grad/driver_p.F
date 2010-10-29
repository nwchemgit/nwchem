C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_P -- Master routine for RI-MP2 density evaluation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_Driver_P(RTDB, Geom, AOBas, NAO,
     $     MaxSpin, TopSpin, NMO, NFrzO, NAct, NVir, FitBas,
     $     g_P2, g_W2, Eig, LDEig, g_AOXF,
     $     Tol2e_Fock, JobTitle, D_Int, D_Gam)
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
      Integer g_P2(TopSpin) ![in]
      Integer g_W2(TopSpin) ![in]
      Integer LDEig ![in]
      Double Precision Eig(LDEig, TopSpin) ![in]
      Integer g_AOXF(TopSpin) ![in]
      Double precision Tol2e_Fock ![in]
      Character*(*) JobTitle ![in]
      Integer D_Int(MaxSpin, 2, 2)
      Integer D_Gam(MaxSpin)
C
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
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk             ! Minimum panel size that is reasonable
      Parameter ( MinChunk = -64)  ! for global distribution of data
C
      Integer MyMaxSpin
      Parameter (MyMaxSpin = 2)
C
C LOCAL VARIABLES
      Integer ISpin, NFit, Length, C, A, V
      Integer g_L(MyMaxspin), g_RHS
      Logical PrInfo, PrPrgRpt
      Character*256 String1, FN_CPHF_RHS, FN_CPHF_Sol
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
     $     'RIMP2_Driver_P: fatal program error: TopSpin > Maxspin',
     $     MaxSpin, INPUT_ERR)
      If ( MaxSpin .ne. MyMaxSpin ) Call ErrQuit(
     $     'RIMP2_Driver_P: fatal program error: Maxspin != MyMaxSpin',
     $     MaxSpin, INPUT_ERR)
C
      If ( .NOT. Bas_NumBF( FitBas, NFit) ) Call ErrQuit(
     $     'RIMP2_Driver_P: can''t get fitting basis size', FitBas,
     &       BASIS_ERR)
C
      PrInfo =   Util_Print('information', Print_Low)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
      If ( (PrInfo .AND. GA_NodeID() .eq. 0) .or. PrPrgRpt ) then
         Write (LuOut, 9000) 
      EndIf
 9000 Format(/1X, 'Density matrix evaluation')
C
C     ****************************
C     * Compute Lagrangian terms *
C     ****************************
C
      Do ISpin = 1, TopSpin
C
C        Create L(a,m)
C
         String1 = 'L spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_P: can''t allocate L', ISpin, GA_ERR)
C
C        Create W2(p,q) (initial contributions come from L)
C
         String1 = 'W2 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_W2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_P: can''t allocate W2', ISpin, GA_ERR)
      EndDo ! ISpin
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Computing Lagrangian (incl. Fock build)')
C
      Call RIMP2_Mk_L(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct, NVir,
     $     NFit, D_Int(1, 1, 1), D_Int(1, 1, 2), D_Gam, g_L, g_P2,
     $     g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
C
C     ************************
C     * Solve CPHF Equations *
C     ************************
C     Produce the relaxed part of the density
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Computing P2(a,m) (CPHF equations)')
C
      Call RIMP2_P2am(RTDB, TopSpin, NFrzO, NAct, NVir, g_L, g_P2)
C
C     Now we can clean up the Lagrangian
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_L(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_P: can''t free L', ISpin, GA_ERR)
      EndDo
C
      If ( Util_Print('partial p2', Print_Debug) .OR.
     $     Util_Print('p2', Print_Debug) ) then
         Call Banner(LuOut, 'Correlated Density (MO basis)', '-',
     $        .TRUE., .TRUE., .FALSE.)
         Do ISpin = 1, TopSpin
            Call GA_Print( g_P2(ISpin) )
         EndDo
      EndIf
C
C     *************************
C     * Save correlated 1-PDM *
C     *************************
C     Correlated one-particle density is now complete, so let's
C     save it, and do some other things with it before going on
C     to complete the gradient.
C
      If ( TopSpin .eq. 2) then
         String1 = 'UHF'
      ElseIf ( TopSpin .eq. 1) then
         String1 = 'RHF'
      Else
         Call ErrQuit( 'RIMP2_Driver_P: TopSpin out of range', TopSpin,
     &       INPUT_ERR)
      EndIf
C
      Call MP2_Make_Dip_Mom(RTDB, AOBas, NAO,
     $     Max(NMO(1), NMO(TopSpin)), NFrzO(1)+NAct(1),
     $     NFrzO(TopSpin)+NAct(TopSpin), g_P2(1), g_P2(TopSpin),
     $     g_AOXF(1), g_AOXF(TopSpin) )
C
      Call MP2_Make_NOs(RTDB, AOBas, NAO, Max(NMO(1), NMO(TopSpin)),
     $     NFrzO(1)+NAct(1), NFrzO(TopSpin)+NAct(TopSpin),
     $     g_P2(1), g_P2(TopSpin), g_AOXF(1), g_AOXF(TopSpin),
     $     String1, Max(NFrzO(1), NFrzO(TopSpin)),
     $     JobTitle(:Inp_StrLen(JobTitle) ) )
C
      Return
      End
