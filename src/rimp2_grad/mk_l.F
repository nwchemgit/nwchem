C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Mk_L -- Form Lagrangian terms & their contributions elsewhere
C
C REVISION
C     $Id: mk_l.F,v 1.5 1998-07-22 00:46:48 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_Mk_L(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct,
     $     NVir, NFit, D_Int_ai, D_Int_ij, D_Gam, g_L, g_P2, g_W2,
     $     Eig, LDEig, g_AOXF, Tol2e_Fock)
      Implicit NONE
C
      Integer RTDB          ![in]
      Integer Geom          ![in]
      Integer AOBas         ![in]
      Integer TopSpin       ![in]
      Integer NFrzO(TopSpin)![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer NFit          ![in]
      Integer D_Int_ai(TopSpin) ![in]
      Integer D_Int_ij(TopSpin) ![in]
      Integer D_Gam(TopSpin) ![in]
      Integer g_L(TopSpin) ![in]
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
      Integer g_L1(MaxSpin), g_L2(MaxSpin), g_L34(MaxSpin),
     $     g_L_tr(MaxSpin)
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
     $     'RIMP2_Mk_L: fatal program error: TopSpin > MaxSpin',
     $     MaxSpin)
C
      Me = GA_NodeID()
C
C     Allocate matrices for partial results
C
      Do ISpin = 1, TopSpin
C
C        Create L1(p,i)
C
         String1 = 'L1 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        NAct(ISpin), String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L1(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L1', ISpin)
C
C        Create L2(a,p)
C
         String1 = 'L2 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl, NVir(ISpin),
     $        NFrzO(ISpin)+NAct(ISpin)+NVir(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L2', ISpin)
C
C        Create L34(a,m)
C
         String1 = 'L34 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NVir(ISpin), NFrzO(ISpin)+NAct(ISpin),
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L34(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t allocate L34', ISpin)
      EndDo ! ISpin
C
C     ********************************
C     * Evaluate L1(p,i) and L2(a,p) *
C     ********************************
C
      Do ISpin = 1, TopSpin
         Call RIMP2_L1(NFrzO(ISpin), NAct(ISpin), NVir(ISpin), NFit,
     $        D_Int_ai(ISpin), D_Gam(ISpin), g_L1(ISpin) )
C
         Call RIMP2_L2(NFrzO(ISpin), NAct(ISpin), NVir(ISpin), NFit,
     $        D_Int_ai(ISpin), D_Int_ij(ISpin), D_Gam(ISpin),
     $        g_L2(ISpin) )
C
C        WARNING:  This is only being done at the moment in order to 
C        match the exact MP2 code!
C
         Call GA_Scale(g_L1(ISpin), -1.0d0) ! CHEATING!
         Call GA_Scale(g_L2(ISpin), -1.0d0) ! CHEATING!
C
      EndDo ! ISpin
C
C     ************************
C     * L1 contributes to P2 *
C     ************************
C     This is separate from the CPHF equations
C
C     P2(iy) = 1/2 L1(yi) (e(i) - e(y))^{-1}
C
C
      If ( NFrzO(ISpin) .gt. 0) then
         Do ISpin = 1, TopSpin
C
            C = NFrzO(ISpin)
            A = NAct(ISpin)
            V = NVir(ISpin)
C
            Call GA_Copy_Patch('T', g_L1(ISpin), 1,   C,   1, A,
     $                              g_P2(ISpin), C+1, C+A, 1, C)
            Call GA_Scale_Patch(g_P2(ISpin), C+1, C+A, 1, C, 0.5d0 )
C
            If ( Util_Print('partial p2', Print_Debug) )
     $           Call GA_Print( g_P2(ISpin) )
C
C           Find out what portion of matrix we hold
C
            Call GA_Distribution(g_P2(ISpin), Me, PLo, PHi, QLo, QHi)
C
C           If we own part of the AC patch, we have work to do
C
            PLo = Max(PLo, C+1)
            PHi = Min(PHi, C+A)
            QLo = Max(QLo, 1)
            QHi = Min(QHi, C)
            If ( (PHi-PLo+1) * (QHi-QLo+1) .gt. 0 ) then
               Call GA_Access(g_P2(ISpin), PLo, PHi, QLo, QHi,
     $              Index, LD)
               Do Y = QLo, QHi
                  Do I = PLo, PHi
                     Scale = 1/( Eig(I, ISpin) - Eig(Y, ISpin) )
                     Dbl_MB( Index ) = Dbl_MB(Index) * Scale
                     Index = Index + 1
                  EndDo
                  Index = Index + LD - 1
               EndDo
               Call GA_Release_Update(g_P2(ISpin), PLo, PHi, QLo, QHi)
            EndIf
            Call GA_Sync
C
C           Print results if requested
C
            If ( Util_Print('partial p2', Print_Debug) )
     $           Call GA_Print( g_P2(ISpin) )
C     
         EndDo
      EndIf
C
C     ******************************************
C     * L3 + L4 is computed using a Fock build *
C     ******************************************
C     Unfortunately, the interface routine from the exact MP2 code
C     produces the transpose of what we want.
C
      Do ISpin = 1, MaxSpin
         String1 = 'transpose L34 spin ' // SpinItoA(ISpin)
         If ( .NOT. GA_Create(MT_Dbl,
     $        NFrzO(ISpin)+NAct(ISpin), NVir(ISpin), 
     $        String1(:Inp_StrLen(String1)), MinChunk,
     $        MinChunk, g_L_tr(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t allocate L_tr', ISpin)
         Call GA_Zero(g_L_tr(ISpin) )
      EndDo
C
C     In RHF cases, we need to play a little trick here, since
C     the MP2 routine demands both spin cases.  Use the RHF
C     values for all inputs of both spin cases.  Get two separate
C     L matrices back, but we can later throw one away.
C
      Beta = Min(TopSpin, MaxSpin)
      Call MP2_Lai_Fock_UHF(Geom, AOBas, g_P2(1), g_P2(Beta),
     $     g_AOXF(1), g_AOXF(Beta), NFrzO(1)+NAct(1),
     $     NFrzO(Beta)+NAct(Beta), NVir(1), NVir(Beta),
     $     g_L_tr(1), g_L_tr(2), RTDB, Tol2E_Fock)
C
      Do ISpin = 1, MaxSpin
         If ( ISpin .le. TopSpin)
     $        Call GA_Transpose(g_L_tr(ISpin), g_L34(ISpin) )
C
         If ( .NOT. GA_Destroy( g_L_tr(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_L: can''t destroy L_tr', ISpin)
      EndDo
C
C     Print results if requested
C
      If ( Util_Print('partial l', Print_Debug) ) then
         Do ISpin = 1, TopSpin
            Call GA_Print( g_L34(ISpin) )
         EndDo
      EndIf
C
C     WARNING:  This is only being done at the moment in order to 
C     match the exact MP2 code!
C
      Do ISpin = 1, TopSpin                 ! CHEATING!
         Call GA_Scale(g_L34(ISpin), 2.0d0) ! CHEATING!
      EndDo                                 ! CHEATING!
C
C
C
C     ********************************
C     * Put the pieces of L together *
C     ********************************
C
      Do ISpin = 1, TopSpin
         Call GA_Zero( g_L(ISpin) )
C
C        Shorthand, so we can actually read the code!
C
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
         If ( C .gt. 0) then
C
C           L(iy) = -L1(yi)
C
            Call GA_Copy_Patch('T', g_L1(ISpin), 1,   C,   1, A,
     $                              g_L(ISpin),  C+1, C+A, 1, C   )
            Call GA_Scale_Patch(g_L(ISpin), C+1, C+A, 1, C, -1.0d0 )
C
C           L(ay) = L2(ay) + L3(ay) + L4(ay)
C           Total size ap       am       am
C
            Call GA_Add_Patch( 1.0d0, g_L2(ISpin),  1, V, 1, C,
     $                         1.0d0, g_L34(ISpin), 1, V, 1, C,
c$$$     $                         0.0d0, g_L34(ISpin), 1, V, 1, C,
     $                                g_L(ISpin),  C+A+1, C+A+V, 1, C)
         EndIf                  ! C .gt. 0
C
C        L(ai) = L1(ai) + L2(ai) + L3(ai) + L4(ai)
C        Total size pi       ap       am       am
C
         Call GA_Add_Patch( 1.0d0, g_L1(ISpin), C+A+1, C+A+V, 1,   A,
     $                      1.0d0, g_L2(ISpin), 1,     V,     C+1, C+A,
     $                             g_L(ISpin),  C+A+1, C+A+V, C+1, C+A)
         Call GA_Add_Patch( 1.0d0, g_L34(ISpin), 1,     V,     C+1, C+A,
     $                      1.0d0, g_L(ISpin),   C+A+1, C+A+V, C+1, C+A,
     $                             g_L(ISpin),   C+A+1, C+A+V, C+1, C+A)
C
C        Print results if requested
C
         If ( Util_Print('partial l', Print_Debug) )
     $        Call GA_Print( g_L(ISpin) )
C
      EndDo ! ISpin
C
C     ***********************
C     * L contributes to W2 *
C     ***********************
C
      Do ISpin = 1, TopSpin
         Call GA_Zero( g_W2(ISpin) )
C
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
C        NOTE: All L contributions to W2 come with a 1/2 factor attached.
C        Also, all L contributions to the (C+A)x(C+A) region carry a minus
C        sign. These will be handled at the end, after the matrix is 
C        assembled.
c$$$C
c$$$C        W2(mi) <-- L1(mi) 
c$$$C
c$$$                           ! + on 0.5d0 is CHEATING!
c$$$         Call GA_Add_Patch(0.5d0, g_L1(ISpin), 1, C+A, 1,   A,
c$$$     $                     1.0d0, g_W2(ISpin), 1, C+A, C+1, C+A,
c$$$     $                            g_W2(ISpin), 1, C+A, C+1, C+A)
C
C        W2(ap) <-- L2(ap)
C
                           ! - on 0.5d0 is CHEATING!
         Call GA_Add_Patch(-0.5d0, g_L2(ISpin), 1,     V,     1, C+A+V,
     $                      1.0d0, g_W2(ISpin), C+A+1, C+A+V, 1, C+A+V,
     $                             g_W2(ISpin), C+A+1, C+A+V, 1, C+A+V)
c$$$C
c$$$C        W2(mi) <-- L3(mi) (Note: L3(ij) = L3(ji))
c$$$C
c$$$         Call GA_Add_Patch(1.0d0, g_L34(ISpin), 1, C+A, 1,   A,
c$$$     $                     1.0d0, g_W2(ISpin),  1, C+A, C+1, C+A,
c$$$     $                            g_W2(ISpin),  1, C+A, C+1, C+A)
c$$$C
c$$$C        W2(yz) <-- L3(yz)
c$$$C
c$$$         If ( C .gt. 0)
c$$$     $        Call GA_Add_Patch(1.0d0, g_L34(ISpin), 1, C, 1, C,
c$$$     $                          1.0d0, g_W2(ISpin),  1, C, 1, C,
c$$$     $                                 g_W2(ISpin),  1, C, 1, C)
c$$$C
c$$$C        W2(om) <-- L4(om)
c$$$C
c$$$         Call GA_Add_Patch(1.0d0, g_L34(ISpin), 1, C+A, 1, C+A,
c$$$     $                     1.0d0, g_W2(ISpin),  1, C+A, 1, C+A,
c$$$     $                            g_W2(ISpin),  1, C+A, 1, C+A)
c$$$C
c$$$C        Handle minus sign and 1/2 factor
c$$$C
c$$$         Call GA_Scale_Patch(g_W2(ISpin), 1, C+A, 1, C+A, -1.0d0)
c$$4         Call GA_Scale(g_W2(ISpin), 0.5d0)
c$$$C
c$$$C        Now handle the two blocks which we know from symmetry arguments
c$$$C        but are not currently correct in W2.
c$$$C
c$$$C        W2(ma) = W2(am)
c$$$C
c$$$         Call GA_Copy_Patch('T',g_W2(ISpin), C+A+1, C+A+V, 1,     C+A,
c$$$     $                          g_W2(ISpin), 1,     C+A,   C+A+1, C+A+V)
c$$$C
c$$$C        W2(iy) = W2(yi)
c$$$C
c$$$         If ( C .gt. 0 ) Call GA_Copy_Patch('T',
c$$$     $        g_W2(ISpin), 1,   C,   C+1, C+A,
c$$$     $        g_W2(ISpin), C+1, C+A, 1,   C)
C
C
C        Print results if requested
C
         If ( Util_Print('partial w2', Print_Debug) )
     $        Call GA_Print( g_W2(ISpin) )
C     
      EndDo                     ! ISpin
C
C     Clean up memory
C
      Do ISpin = 1, TopSpin
         If ( .NOT. GA_Destroy( G_L34(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L34', ISpin)
         If ( .NOT. GA_Destroy( G_L2(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L2', ISpin)
         If ( .NOT. GA_Destroy( G_L1(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Driver_G: can''t free L1', ISpin)
      EndDo
C
      Return
      End
