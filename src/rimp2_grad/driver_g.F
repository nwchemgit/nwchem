C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_G -- Master routine for RI-MP2 gradient evaluation
C
C REVISION
C     $Id: driver_g.F,v 1.5 1998-07-19 00:08:30 bernhold Exp $
C
C SYNOPSIS
      Logical Function RIMP2_Driver_G(RTDB, Geom, AOBas, MaxSpin,
     $     TopSpin, NFrzO, NAct, NVir, FitBas, FN_Int, FN_Gam,
     $     BraKetSame, g_P2, Eig, LDEig, g_AOXF, Tol2e_Fock)
      Implicit NONE
C
      Integer RTDB          ![in]
      Integer Geom          ![in]
      Integer AOBas         ![in]
      Integer MaxSpin       ![in]
      Integer TopSpin       ![in]
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
C
C     **************
C     * Open files *
C     **************
C     Open both bra and ket files if they are distinct, otherwise
C     simply duplicate the bra handle.
C
      Do ISpin = 1, TopSpin
C
C        Integrals
C
         If ( DRA_Open(
     $        FN_Int(ISpin, 1, 1)(:Inp_StrLen(FN_Int(ISpin, 1, 1))),
     $        DRA_R, D_Int(ISpin, 1, 1) ) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_G:  can''t open integral file for spin ',
     $        ISpin)
         If ( BraKetSame ) then
            D_Int(ISpin, 2, 1) = D_Int(ISpin, 1, 1)
         Else
            If ( DRA_Open(
     $           FN_Int(ISpin, 2, 1)(:Inp_StrLen(FN_Int(ISpin, 2, 1))),
     $           DRA_R, D_Int(ISpin, 2, 1) ) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_G:  can''t open integral file for spin ',
     $        ISpin)
         EndIf
         If ( DRA_Open(
     $        FN_Int(ISpin, 1, 2)(:Inp_StrLen(FN_Int(ISpin, 1, 2))),
     $        DRA_R, D_Int(ISpin, 1, 2) ) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_G:  can''t open integral file for spin ',
     $        ISpin)
         If ( DRA_Open(FN_Gam(ISpin)(:Inp_StrLen(FN_Gam(ISpin))),
     $        DRA_R, D_Gam(ISpin) ) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_G:  can''t open gamma file for spin ',
     $        ISpin)
      EndDo                     ! ISpin
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
     $        'RIMP2_Driver_G: can''t allocate L', ISpin)
C
C        Create W2(p,q) (initial contributions come from L)
C
         String1 = 'W2 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_W2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t allocate W2', ISpin)
      EndDo ! ISpin
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Computing Lagrangian (incl. Fock build)')
      Call RIMP2_Mk_L(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct, NVir,
     $     NFit, D_Int(1, 1, 1), D_Int(1, 1, 2), D_Gam, g_L, g_P2,
     $     g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
C
C     ************************
C     * Solve CPHF Equations *
C     ************************
C
      If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $     '    Solving CPHF equations')
C
C     Create array to hold RHS
C
      Length = 0
      Do ISpin = 1, TopSpin
         Length = Length + (NVir(ISpin) * (NFrzO(ISpin)+NAct(ISpin) ) )
      EndDo
C
      If ( .NOT. GA_Create(MT_Dbl, Length, 1, 'RI-MP2 CPHF RHS/Soln', 
     $     1, 1, g_RHS ) ) Call ErrQuit(
     $     'RIMP2_Driver_G: can''t allocate RHS', 0)
C
C     Copy appropriate bits of L into RHS
C
      Length = 0
      Do ISpin = 1, TopSpin
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
         Call GA_Copy_Patch('N', g_L(ISpin), C+A+1, C+A+V, 1, C+A,
     $        g_RHS, 1+Length, (V*(C+A))+Length, 1, 1)
C
C        Need length of alpha spin for next spin case
C
         Length = Length + (V * (C+A) )
      EndDo
C
      If ( Util_Print('partial p2', Print_Debug) .OR.
     $     Util_Print('cphf', Print_Debug) ) Call GA_Print(g_RHS)
C
C     CPHF solver expects RHS on a file.  Name is fixed by CPHF module.
C
      Call Util_File_name('cphf_rhs', .TRUE., .TRUE., FN_CPHF_RHS)
      If (.NOT. File_Write_GA(FN_CPHF_RHS, g_RHS) ) Call ErrQuit(
     $     'RIMP2_Driver_G: can''t write CPHF RHS to file', 0)
C
C     Invoke the CPHF solver
C
      If ( .NOT. CPHF( RTDB ) ) Call ErrQuit(
     $     'RIMP2_Driver_G: cphf solution failed', 0)
C
C     Once again, solution is returned in a file, name fixed by CPHF module
C
      Call Util_File_Name('cphf_sol', .TRUE., .TRUE., FN_CPHF_Sol)
      If ( .NOT. File_Read_GA(FN_CPHF_Sol, g_rhs) ) Call ErrQuit(
     $     'RIMP2_Driver_G: can''t read CPHF sol''n file', 0)
C
      If ( Util_Print('partial p2', Print_Debug) .OR.
     $     Util_Print('cphf', Print_Debug) ) Call GA_Print(g_RHS)
c$$$C
c$$$C     Copy appropriate bits of L into RHS
c$$$C
c$$$      Length = 0
c$$$      Do ISpin = 1, TopSpin
c$$$         C = NFrzO(ISpin)
c$$$         A = NAct(ISpin)
c$$$         V = NVir(ISpin)
c$$$C
c$$$         Call GA_Copy_Patch('N', g_L(ISpin), C+A+1, C+A+V, 1, C+A,
c$$$     $        g_RHS, 1+Length, (V*(C+A))+Length, 1, 1)
c$$$C
c$$$C        Need length of alpha spin for next spin case
c$$$C
c$$$         Length = Length + (V * (C+A) )
c$$$      EndDo
C
C     Free temporary memory
C
      If ( .NOT. GA_Destroy( g_RHS) ) Call ErrQuit(
     $     'RIMP2_Driver_G: can''t destroy CPHF RHS GA', 0)
C
C     *******************
C     * Clean up memory *
C     *******************
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_L(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L', ISpin)
         If ( .NOT. GA_Destroy( G_W2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free W2', ISpin)
C
C        Created in RIMP2_Driver_E
C
         If ( .NOT. GA_Destroy( G_P2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free P2', ISpin)
      EndDo
C
C     ************************
C     * Close & delete files *
C     ************************
C
      Do ISpin = 1, TopSpin
         If ( DRA_Delete(D_Int(ISpin, 1, 1)) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_E: can''t delete integral file', ISpin)
         If ( .NOT. BraKetSame ) then
            If ( DRA_Delete(D_Int(ISpin, 2, 1)) .ne. 0) Call ErrQuit(
     $           'RIMP2_Driver_E: can''t delete integral file', ISpin)
         EndIf
         If ( DRA_Delete(D_Int(ISpin, 1, 2)) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_E: can''t delete integral file', ISpin)
         If ( DRA_Delete(D_Gam(ISpin)) .ne. 0) Call ErrQuit(
     $        'RIMP2_Driver_E: can''t delete gamma file', ISpin)
      EndDo    ! ISpin
C
      Return
      End
