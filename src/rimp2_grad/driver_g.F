C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_G -- Master routine for RI-MP2 gradient evaluation
C
C REVISION
C     $Id: driver_g.F,v 1.2 1998-07-13 02:38:08 bernhold Exp $
C
C SYNOPSIS
      Logical Function RIMP2_Driver_G(TopSpin, NFrzO, NAct, NVir,
     $     FitBas)
      Implicit NONE
C
      Integer TopSpin ![in]
      Integer NFrzO ![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer FitBas ![in]
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk             ! Minimum panel size that is reasonable
      Parameter ( MinChunk = -64)  ! for global distribution of data
C
      Integer MaxSpin
      Parameter (MaxSpin = 2)
C
C LOCAL VARIABLES
      Integer ISpin
      Integer g_L1(MaxSpin)
      Logical PrInfo, PrPrgRpt
      Character*256 String1
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
C
      PrInfo =   Util_Print('information', Print_Low)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
      If ( (PrInfo .AND. GA_NodeID() .eq. 0) .or. PrPrgRpt ) then
         Write (LuOut, 9000) 
      EndIf
 9000 Format(1X, 'Gradient evaluation')
C
C     ******************************
C     * Create and compute L1(p,i) *
C     ******************************
C
      Do ISpin = 1, TopSpin
         String1 = 'L1 spin ' // SpinItoA(ISpin)
C
         If ( .NOT. GA_Create(MT_Dbl, NFrzO+NAct(ISpin)+NVir(ISpin),
     $        NAct(ISpin), String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L1(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t allocate L1', ISpin)
         Call GA_Zero( g_L1(ISpin) )
C
         Call RIMP2_L1()
      EndDo ! ISpin for L1
C
C     Clean up memory we no longer need
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_L1(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L1', ISpin)
      EndDo
C
      Return
      End
