C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_G -- Master routine for RI-MP2 gradient evaluation
C
C REVISION
C     $Id: driver_g.F,v 1.8 1998-10-06 17:19:11 bernhold Exp $
C
C SYNOPSIS
      Logical Function RIMP2_Driver_G(RTDB, Geom, AOBas, NAO,
     $     MaxSpin, TopSpin, NMO, NFrzO, NAct, NVir, FitBas,
     $     FN_Int, FN_Gam, BraKetSame, g_P2, Eig, LDEig, g_AOXF,
     $     Tol2e_Fock, JobTitle)
      Implicit NONE
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
      Character*(*) FN_Int(MaxSpin, 2, 2) ![in]
      Character*(*) FN_Gam(TopSpin) ![in]
      Logical BraKetSame    ![in]
      Integer g_P2(TopSpin) ![in]
      Integer LDEig ![in]
      Double Precision Eig(LDEig, TopSpin) ![in]
      Integer g_AOXF(TopSpin) ![in]
      Double precision Tol2e_Fock ![in]
      Character*(*) JobTitle ![in]
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
      Integer g_L(MyMaxspin), g_W2(MyMaxSpin), g_RHS
      Logical PrInfo, PrPrgRpt
      Character*256 String1, FN_CPHF_RHS, FN_CPHF_Sol
      Integer D_Int(MyMaxSpin, 2, 2), D_Gam(MyMaxSpin)
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
     $     MaxSpin)
      If ( MaxSpin .ne. MyMaxSpin ) Call ErrQuit(
     $     'RIMP2_Driver_G: fatal program error: Maxspin != MyMaxSpin',
     $     MaxSpin)
C
      If ( .NOT. Bas_NumBF( FitBas, NFit) ) Call ErrQuit(
     $     'RIMP2_Driver_G: can''t get fitting basis size', FitBas)
C
      PrInfo =   Util_Print('information', Print_Low)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
      If ( (PrInfo .AND. GA_NodeID() .eq. 0) .or. PrPrgRpt ) then
         Write (LuOut, 9000) 
      EndIf
 9000 Format(/1X, 'Gradient evaluation')
c$$$C
c$$$C     **************
c$$$C     * Open files *
c$$$C     **************
c$$$C     Open both bra and ket files if they are distinct, otherwise
c$$$C     simply duplicate the bra handle.
c$$$C
c$$$      Do ISpin = 1, TopSpin
c$$$C
c$$$C        Integrals
c$$$C
c$$$         If ( DRA_Open(
c$$$     $        FN_Int(ISpin, 1, 1)(:Inp_StrLen(FN_Int(ISpin, 1, 1))),
c$$$     $        DRA_R, D_Int(ISpin, 1, 1) ) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G:  can''t open integral file for spin ',
c$$$     $        ISpin)
c$$$         If ( BraKetSame ) then
c$$$            D_Int(ISpin, 2, 1) = D_Int(ISpin, 1, 1)
c$$$         Else
c$$$            If ( DRA_Open(
c$$$     $           FN_Int(ISpin, 2, 1)(:Inp_StrLen(FN_Int(ISpin, 2, 1))),
c$$$     $           DRA_R, D_Int(ISpin, 2, 1) ) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G:  can''t open integral file for spin ',
c$$$     $        ISpin)
c$$$         EndIf
c$$$         If ( DRA_Open(
c$$$     $        FN_Int(ISpin, 1, 2)(:Inp_StrLen(FN_Int(ISpin, 1, 2))),
c$$$     $        DRA_R, D_Int(ISpin, 1, 2) ) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G:  can''t open integral file for spin ',
c$$$     $        ISpin)
c$$$         If ( DRA_Open(FN_Gam(ISpin)(:Inp_StrLen(FN_Gam(ISpin))),
c$$$     $        DRA_R, D_Gam(ISpin) ) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G:  can''t open gamma file for spin ',
c$$$     $        ISpin)
c$$$      EndDo                     ! ISpin
c$$$C
c$$$C     ****************************
c$$$C     * Compute Lagrangian terms *
c$$$C     ****************************
c$$$C
c$$$      Do ISpin = 1, TopSpin
c$$$C
c$$$C        Create L(a,m)
c$$$C
c$$$         String1 = 'L spin ' // SpinItoA(ISpin)
c$$$         If ( .NOT. GA_Create(MT_Dbl,
c$$$     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
c$$$     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
c$$$     $        String1(:Inp_StrLen(String1)), MinChunk,
c$$$     $        MinChunk, g_L(ISpin) ) ) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G: can''t allocate L', ISpin)
c$$$C
c$$$C        Create W2(p,q) (initial contributions come from L)
c$$$C
c$$$         String1 = 'W2 spin ' // SpinItoA(ISpin)
c$$$         If ( .NOT. GA_Create(MT_Dbl,
c$$$     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
c$$$     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
c$$$     $        String1(:Inp_StrLen(String1)), MinChunk,
c$$$     $        MinChunk, g_W2(ISpin) ) ) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G: can''t allocate W2', ISpin)
c$$$      EndDo ! ISpin
c$$$C
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $     '    Computing Lagrangian (incl. Fock build)')
c$$$C
c$$$      Call RIMP2_Mk_L(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct, NVir,
c$$$     $     NFit, D_Int(1, 1, 1), D_Int(1, 1, 2), D_Gam, g_L, g_P2,
c$$$     $     g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
c$$$C
c$$$C     ************************
c$$$C     * Solve CPHF Equations *
c$$$C     ************************
c$$$C     Produce the relaxed part of the density
c$$$C
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $     '    Computing P2(a,m) (CPHF equations)')
c$$$C
c$$$      Call RIMP2_P2am(RTDB, TopSpin, NFrzO, NAct, NVir, g_L, g_P2)
c$$$C
c$$$C     Now we can clean up the Lagrangian
c$$$C
c$$$      Do ISpin = 1, TopSpin
c$$$         If ( .NOT. GA_Destroy( G_L(ISpin) ) ) Call ErrQuit(
c$$$     $        'RIMP2_Driver_G: can''t free L', ISpin)
c$$$      EndDo
c$$$C
c$$$C     *************************
c$$$C     * Save correlated 1-PDM *
c$$$C     *************************
c$$$C     Correlated one-particle density is now complete, so let's
c$$$C     save it, and do some other things with it before going on
c$$$C     to complete the gradient.
c$$$C
c$$$      If ( TopSpin .eq. 2) then
c$$$         String1 = 'UHF'
c$$$      ElseIf ( TopSpin .eq. 1) then
c$$$         String1 = 'RHF'
c$$$      Else
c$$$         Call ErrQuit( 'RIMP2_Driver_G: TopSpin out of range', TopSpin)
c$$$      EndIf
c$$$C
c$$$      Call MP2_Make_Dip_Mom(RTDB, AOBas, NAO,
c$$$     $     Max(NMO(1), NMO(TopSpin)), NFrzO(1)+NAct(1),
c$$$     $     NFrzO(TopSpin)+NAct(TopSpin), g_P2(1), g_P2(TopSpin),
c$$$     $     g_AOXF(1), g_AOXF(TopSpin) )
c$$$C
c$$$      Call MP2_Make_NOs(RTDB, AOBas, NAO, Max(NMO(1), NMO(TopSpin)),
c$$$     $     NFrzO(1)+NAct(1), NFrzO(TopSpin)+NAct(TopSpin),
c$$$     $     g_P2(1), g_P2(TopSpin), g_AOXF(1), g_AOXF(TopSpin),
c$$$     $     String1, Max(NFrzO(1), NFrzO(TopSpin)),
c$$$     $     JobTitle(:Inp_StrLen(JobTitle) ) )
C
C     ***************
C     * Complete W2 *
C     ***************
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Completing energy-weighted density (incl. Fock build)')
      Call RIMP2_Mk_W2(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct, NVir,
     $     g_P2, g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
C
C     *******************
C     * Clean up memory *
C     *******************
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_W2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free W2', ISpin)
C
C        Created in RIMP2_Driver_E
C
         If ( .NOT. GA_Destroy( G_P2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free P2', ISpin)
      EndDo
c$$$C
c$$$C     ************************
c$$$C     * Close & delete files *
c$$$C     ************************
c$$$C
c$$$      Do ISpin = 1, TopSpin
c$$$         If ( DRA_Delete(D_Int(ISpin, 1, 1)) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_E: can''t delete integral file', ISpin)
c$$$         If ( .NOT. BraKetSame ) then
c$$$            If ( DRA_Delete(D_Int(ISpin, 2, 1)) .ne. 0) Call ErrQuit(
c$$$     $           'RIMP2_Driver_E: can''t delete integral file', ISpin)
c$$$         EndIf
c$$$         If ( DRA_Delete(D_Int(ISpin, 1, 2)) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_E: can''t delete integral file', ISpin)
c$$$         If ( DRA_Delete(D_Gam(ISpin)) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_Driver_E: can''t delete gamma file', ISpin)
c$$$      EndDo    ! ISpin
C
      Return
      End
