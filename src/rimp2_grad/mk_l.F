C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Mk_L -- Form Lagrangian terms & their contributions elsewhere
C
C REVISION
C     $Id: mk_l.F,v 1.1 1998-07-14 16:18:36 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_Mk_L(TopSpin, NFrzO, NAct, NVir, NFit,
     $     D_Int_ai, D_Int_ij, D_Gam, g_L)
      Implicit NONE
C
      Integer TopSpin       ![in]
      Integer NFrzO         ![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer NFit          ![in]
      Integer D_Int_ai(TopSpin) ![in]
      Integer D_Int_ij(TopSpin) ![in]
      Integer D_Gam(TopSpin) ![in]
      Integer g_L(TopSpin) ![in]
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c$$$#include "stdio.fh"
#include "inp.fh"
c$$$#include "dra.fh"
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
      Integer g_L1(MaxSpin), g_L2(MaxSpin), g_L3(MaxSpin),g_L4(MaxSpin)
      Character*(256) String1
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
C
C     Sanity checks and initialization
C
      If ( TopSpin .gt. MaxSpin) Call ErrQuit(
     $     'RIMP2_Mk_L: fatal program error: TopSpin > MaxSpin',
     $     MaxSpin)
C
C     Allocate matrices for partial results
C
      Do ISpin = 1, TopSpin
C
C        Create L1(p,i)
C
         String1 = 'L1 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl, NFrzO+NAct(ISpin)+NVir(ISpin),
     $        NAct(ISpin), String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L1(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L1', ISpin)
C
C        Create L2(a,p)
C
         String1 = 'L2 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl, NVir(ISpin),
     $        NFrzO+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L2', ISpin)
C
C        Create L3(p,n)
C
         String1 = 'L3 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl, NFrzO+NAct(ISpin)+NVir(ISpin),
     $        NFrzO+NAct(ISpin), String1(:Inp_StrLen(String1)),
     $        MinChunk, MinChunk, g_L3(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L3', ISpin)
C
C        Create L4(p,q)
C
         String1 = 'L4 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl, NFrzO+NAct(ISpin)+NVir(ISpin),
     $        NFrzO+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L4(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L4', ISpin)
      EndDo ! ISpin
C
C     ********************************
C     * Evaluate L1(p,i) and L2(a,p) *
C     ********************************
C
      Do ISpin = 1, TopSpin
         Call RIMP2_L1(NFrzO, NAct(ISpin), NVir(ISpin), NFit,
     $        D_Int_ai(ISpin), D_Gam(ISpin), g_L1(ISpin) )
C
         Call RIMP2_L2(NFrzO, NAct(ISpin), NVir(ISpin), NFit,
     $        D_Int_ai(ISpin), D_Int_ij(ISpin), D_Gam(ISpin),
     $        g_L2(ISpin) )
      EndDo ! ISpin
C
C     ******************************
C     * Create and compute L3(p,n) *
C     ******************************
C
      Call GA_Zero( g_L3(1) )
      Call GA_Zero( g_L3(2) )
C
C     ******************************
C     * Create and compute L4(p,q) *
C     ******************************
C
      Call GA_Zero( g_L4(1) )
      Call GA_Zero( g_L4(2) )
C
C     ********************************
C     * Put the pieces of L together *
C     ********************************
C
      Do ISpin = 1, TopSpin
         Call GA_Zero( g_L(ISpin) )
C
         If ( NFrzO .gt. 0) then
C
C           L(iy) = -L1(yi)
C
            Call GA_Copy('T',
     $           g_L1(ISpin), 1, NFrzO, NFrzO+1, NFrzO+NAct(ISpin),
     $           g_L(ISpin), NFrzO+1, NFrzO+NAct(ISpin), 1, NFrzO)
            Call GA_Scale_Patch(g_L(ISpin), NFrzO+1, NFrzO+NAct(ISpin),
     $           1, NFrzO, -1.0d0 )
C
C           L(ay) = L2(ay) + L3(ay) + L4(ay)
C           Total size ap       pn       pq
C
            Call GA_Add_Patch(
     $           1.0d0, g_L2(ISpin), 1, NVir(ISpin), 1, NFrzO,
     $           1.0d0, g_L3(ISpin), NFrzO+NAct(ISpin)+1,
     $           NFrzO+NAct(ISpin)+NVir(ISpin), 1, NFrzO,
     $           g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $           NFrzO+NAct(ISpin)+NVir(ISpin), 1, NFrzO)
            Call GA_Add_Patch( 1.0d0, g_L4(ISpin), NFrzO+NAct(ISpin)+1,
     $           NFrzO+NAct(ISpin)+NVir(ISpin), 1, NFrzO,
     $           1.0d0, g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $           NFrzO+NAct(ISpin)+NVir(ISpin), 1, NFrzO,
     $           g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $           NFrzO+NAct(ISpin)+NVir(ISpin), 1, NFrzO)
         EndIf                  ! NFrzO .gt. 0
C
C        L(ai) = L1(ai) + L2(ai) + L3(ai) + L4(ai)
C        Total size pi       ap       pn       pq
C
            Call GA_Add_Patch(
     $        1.0d0, g_L1(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), 1, NAct(ISpin),
     $        1.0d0, g_L2(ISpin), 1, NVir(ISpin), NFrzO+1,
     $        NFrzO+NAct(ISpin),
     $        g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin))
            Call GA_Add_Patch(
     $        1.0d0, g_L3(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin),
     $        1.0d0, g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin),
     $        g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin))
            Call GA_Add_Patch(
     $        1.0d0, g_L4(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin),
     $        1.0d0, g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin),
     $        g_L(ISpin), NFrzO+NAct(ISpin)+1,
     $        NFrzO+NAct(ISpin)+NVir(ISpin), NFrzO+1, NFrzO+NAct(ISpin))
      EndDo ! ISpin
C
C     Print results if requested
C
      If ( Util_Print('partial l', Print_Debug) ) then
         Do ISpin = 1, TopSpin
            Call GA_Print( g_L(ISpin) )
         EndDo
      EndIf
C
C     Clean up memory
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_L4(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L4', ISpin)
         If ( .NOT. GA_Destroy( G_L3(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L3', ISpin)
         If ( .NOT. GA_Destroy( G_L2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L2', ISpin)
         If ( .NOT. GA_Destroy( G_L1(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L1', ISpin)
      EndDo
C
      Return
      End
