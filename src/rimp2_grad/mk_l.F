C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Mk_L -- Form Lagrangian terms & their contributions elsewhere
C
C REVISION
C     $Id: mk_l.F,v 1.2 1998-07-14 21:41:49 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_Mk_L(TopSpin, NFrzO, NAct, NVir, NFit,
     $     D_Int_ai, D_Int_ij, D_Gam, g_L, g_P2, g_W2, Eig, LDEig)
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
      Integer g_P2(TopSpin) ![in]
      Integer g_W2(TopSpin) ![in]
      Integer LDEig ![in]
      Double Precision Eig(LDEig, TopSpin) ![in]
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
      Integer C, A, V, PLo, PHi, QLo, QHi, I, Y, Index, LD
      Integer Me
      Double Precision Scale
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
C        Shorthand, so we can actually read the code!
C
         C = NFrzO
         A = NAct(ISpin)
         V = NVir(ISpin)
C
         If ( NFrzO .gt. 0) then
C
C           L(iy) = -L1(yi)
C
            Call GA_Copy_Patch('T', g_L1(ISpin), 1,   C,   C+1, C+A,
     $                              g_L(ISpin),  C+1, C+A, 1,   C   )
            Call GA_Scale_Patch(g_L(ISpin), C+1, C+A, 1, C, -1.0d0 )
C
C           L(ay) = L2(ay) + L3(ay) + L4(ay)
C           Total size ap       pn       pq
C
            Write (6, *) 'Add patch in L(ay)'
            Call GA_Add_Patch( 1.0d0, g_L2(ISpin), 1,     V,     1, C,
     $                         1.0d0, g_L3(ISpin), C+A+1, C+A+V, 1, C,
     $                                g_L(ISpin),  C+A+1, C+A+V, 1, C)
            Call GA_Add_Patch( 1.0d0, g_L4(ISpin), C+A+1, C+A+V, 1, C,
     $                         1.0d0, g_L(ISpin),  C+A+1, C+A+V, 1, C,
     $                                g_L(ISpin),  C+A+1, C+A+V, 1, C)
         EndIf                  ! NFrzO .gt. 0
C
C        L(ai) = L1(ai) + L2(ai) + L3(ai) + L4(ai)
C        Total size pi       ap       pn       pq
C
         Write (6, *) 'Add patch in L(ai)'
         Call GA_Add_Patch( 1.0d0, g_L1(ISpin), C+A+1, C+A+V, 1,   A,
     $                      1.0d0, g_L2(ISpin), 1,     V,     C+1, C+A,
     $                             g_L(ISpin),  C+A+1, C+A+V, C+1, C+A)
         Call GA_Add_Patch( 1.0d0, g_L3(ISpin), C+A+1, C+A+V, C+1, C+A,
     $                      1.0d0, g_L(ISpin),  C+A+1, C+A+V, C+1, C+A,
     $                             g_L(ISpin),  C+A+1, C+A+V, C+1, C+A)
         Call GA_Add_Patch( 1.0d0, g_L4(ISpin), C+A+1, C+A+V, C+1, C+A,
     $                      1.0d0, g_L(ISpin),  C+A+1, C+A+V, C+1, C+A,
     $                             g_L(ISpin),  C+A+1, C+A+V, C+1, C+A)
C
C        Print results if requested
C
         If ( Util_Print('partial l', Print_Debug) )
     $        Call GA_Print( g_L(ISpin) )
C
      EndDo ! ISpin
C
C     ***********************
C     * L contributes to P2 *
C     ***********************
C     This is separate from the CPHF equations
C
C     P2(iy) = 1/2 L1(yi) (e(i) - e(y))^{-1}
C
C
      If ( NFrzO .gt. 0) then
         Do ISpin = 1, TopSpin
C
            C = NFrzO
            A = NAct(ISpin)
            V = NVir(ISpin)
C
            Call GA_Copy_Patch('T', g_L1(ISpin), 1,   C,   C+1, C+A,
     $                              g_P2(ISpin), C+1, C+A, 1,   C)
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
                     Write ( 6, *) 'I, Y, Scale = ', I, Y, Scale
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
C     ***********************
C     * L contributes to W2 *
C     ***********************
C
      Do ISpin = 1, TopSpin
         Call GA_Zero( g_W2(ISpin) )
C
         C = NFrzO
         A = NAct(ISpin)
         V = NVir(ISpin)
C
C        NOTE: All L contributions to W2 come with a 1/2 factor attached.
C        Also, all L contributions to the (C+A)x(C+A) region carry a minus
C        sign. These will be handled at the end, after the matrix is 
C        assembled.
C
C        W2(mi) <-- L1(mi) 
C
         Write (6, *) 'Add patch in W2(mi)'
         Call GA_Add_Patch(1.0d0, g_L1(ISpin), 1, C+A, 1,   A,
     $                     1.0d0, g_W2(ISpin), 1, C+A, C+1, C+A,
     $                            g_W2(ISpin), 1, C+A, C+1, C+A)
C
C        W2(ap) <-- L2(ap)
C
         Write (6, *) 'Add patch in W2(ap)'
         Call GA_Add_Patch(1.0d0, g_L2(ISpin), 1,     V,     1, C+A+V,
     $                     1.0d0, g_W2(ISpin), C+A+1, C+A+V, 1, C+A+V,
     $                            g_W2(ISpin), C+A+1, C+A+V, 1, C+A+V)
C
C        W2(mi) <-- L3(mi) (Note: L3(ij) = L3(ji))
C
         Write (6, *) 'Add patch in W2(mi)'
         Call GA_Add_Patch(1.0d0, g_L3(ISpin), 1, C+A, 1,   A,
     $                     1.0d0, g_W2(ISpin), 1, C+A, C+1, C+A,
     $                            g_W2(ISpin), 1, C+A, C+1, C+A)
C
C        W2(yz) <-- L3(yz)
C
         Write (6, *) 'Add patch in W2(yz)'
         If ( C .gt. 0)
     $        Call GA_Add_Patch(1.0d0, g_L3(ISpin), 1, C, 1, C,
     $                          1.0d0, g_W2(ISpin), 1, C, 1, C,
     $                                 g_W2(ISpin), 1, C, 1, C)
C
C        W2(om) <-- L4(om)
C
         Write (6, *) 'Add patch in W2(om)'
         Call GA_Add_Patch(1.0d0, g_L4(ISpin), 1, C+A, 1, C+A,
     $                     1.0d0, g_W2(ISpin), 1, C+A, 1, C+A,
     $                            g_W2(ISpin), 1, C+A, 1, C+A)
C
C        Handle minus sign and 1/2 factor
C
         Call GA_Scale_Patch(g_W2(ISpin), 1, C+A, 1, C+A, -1.0d0)
         Call GA_Scale(g_W2(ISpin), -0.5d0)
C
C        Now handle the two blocks which we know from symmetry arguments
C        but are not currently correct in W2.
C
C        W2(ma) = W2(am)
C
         Call GA_Copy_Patch('T',g_W2(ISpin), C+A+1, C+A+V, 1,     C+A,
     $                          g_W2(ISpin), 1,     C+A,   C+A+1, C+A+V)
C
C        W2(iy) = W2(yi)
C
         If ( C .gt. 0 ) Call GA_Copy_Patch('T',
     $        g_W2(ISpin), 1,   C,   C+1, C+A,
     $        g_W2(ISpin), C+1, C+A, 1,   C)
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
