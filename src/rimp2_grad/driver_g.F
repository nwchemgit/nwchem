C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_G -- Master routine for RI-MP2 gradient evaluation
C
C REVISION
C     $Id: driver_g.F,v 1.3 1998-07-14 16:18:34 bernhold Exp $
C
C SYNOPSIS
      Logical Function RIMP2_Driver_G(MaxSpin, TopSpin, NFrzO, NAct,
     $     NVir, FitBas, FN_Int, FN_Gam, BraKetSame)
      Implicit NONE
C
      Integer MaxSpin       ![in]
      Integer TopSpin       ![in]
      Integer NFrzO         ![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer FitBas        ![in]
      Character*(*) FN_Int(MaxSpin, 2, 2) ![in]
      Character*(*) FN_Gam(TopSpin) ![in]
      Logical BraKetSame    ![in]
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
      Integer ISpin, NFit
      Integer g_L1(MyMaxSpin), g_L2(MyMaxSpin), g_L3(MyMaxSpin),
     $     g_L4(MyMaxSpin), g_L(MyMaxspin)
      Logical PrInfo, PrPrgRpt
      Character*256 String1
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
C        Create L(p,q)
C
         String1 = 'L spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl, NFrzO+NAct(ISpin)+NVir(ISpin),
     $        NFrzO+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t allocate L', ISpin)
      EndDo ! ISpin
C
      Call RIMP2_Mk_L(TopSpin, NFrzO, NAct, NVir, NFit, D_Int(1, 1, 1),
     $     D_Int(1, 1, 2), D_Gam, g_L)
C
C     *******************
C     * Clean up memory *
C     *******************
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_L(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L', ISpin)
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
