C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Mk_W2 -- Complete formation of W2
C
C REVISION
C     $Id: mk_w2.F,v 1.1 1998-07-22 00:46:48 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_Mk_W2(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct,
     $     g_P2, g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
      Implicit NONE
C
      Integer RTDB          ![in]
      Integer Geom          ![in]
      Integer AOBas         ![in]
      Integer TopSpin       ![in]
      Integer NFrzO(TopSpin)![in]
      Integer NAct(TopSpin) ![in]
      Integer g_P2(TopSpin) ![in]
      Integer g_W2(TopSpin) ![in]
      Integer LDEig ![in]
      Double Precision Eig(LDEig, TopSpin) ![in]
      Integer g_AOXF(TopSpin) ![in]
      Double precision Tol2e_Fock ![in]
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "rtdb.fh"
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
      Integer g_W2ij(MaxSpin)
      Character*(256) String1
      Integer C, A, V, PLo, PHi, QLo, QHi, I, Y, Index, LD
      Integer Me, Beta
      Double Precision Scale
      Logical Status
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
C
C     Sanity checks and initialization
C
      If ( TopSpin .gt. MaxSpin) Call ErrQuit(
     $     'RIMP2_Mk_W2: fatal program error: TopSpin > MaxSpin',
     $     MaxSpin)
C
      Me = GA_NodeID()
C
C     *****************************************
C     * W2(ij) is computed using a Fock build *
C     *****************************************
C
      Do ISpin = 1, MaxSpin
         String1 = 'W2(ij) ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NFrzO(ISpin)+NAct(ISpin), NFrzO(ISpin)+NAct(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_W2ij(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t allocate W2ij', ISpin)
         Call GA_Zero(g_W2ij(ISpin) )
      EndDo
C
C     In RHF cases, we need to play a little trick here, since
C     the MP2 routine demands both spin cases.  Use the RHF
C     values for all inputs of both spin cases.  Get two separate
C     L matrices back, but we can later throw one away.
C
      Beta = Min(TopSpin, MaxSpin)
      Call MP2_Wij_Fock_UHF(RTDB, Geom, AOBas, Tol2e_Fock, g_P2(1),
     $     g_P2(Beta), NFrzO(1)+NAct(1), NFrzO(Beta)+NAct(Beta),
     $     g_AOXF(1), g_AOXF(Beta), g_W2ij(1), g_W2ij(2) )
C
C     Print results if requested
C
      If ( Util_Print('partial w2', Print_Debug) ) then
         Do ISpin = 1, TopSpin
            Call GA_Print( g_W2ij(ISpin) )
         EndDo
      EndIf
C
      Do ISpin = 1, MaxSpin
         If ( .NOT. GA_Destroy( g_W2ij(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_W2: can''t destroy L_tr', ISpin)
      EndDo
C
      Return
      End











