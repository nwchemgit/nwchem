C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Mk_W2 -- Complete formation of W2
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Mk_W2(RTDB, Geom, AOBas, TopSpin, NFrzO, NAct,
     $     NVir, g_P2, g_W2, Eig, LDEig, g_AOXF, Tol2e_Fock)
      Implicit NONE
#include "errquit.fh"
C
      Integer RTDB          ![in]
      Integer Geom          ![in]
      Integer AOBas         ![in]
      Integer TopSpin       ![in]
      Integer NFrzO(TopSpin)![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(Topspin) ![in]
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
#include "stdio.fh"
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
      Integer g_W2ij(MaxSpin), g_Tmp
      Character*(256) String1
      Integer C, A, V, PLo, PHi, QLo, QHi, I, Y, LD, P, Q, i_Tmp
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
     $     MaxSpin, INPUT_ERR)
C
      Me = GA_NodeID()
c$$$C
c$$$      Do ISpin = 1, TopSpin
c$$$         Call GA_Zero(g_W2(ISpin) )
c$$$      EndDo
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
     $        'RIMP2_Driver_G: can''t allocate W2ij', ISpin, GA_ERR)
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
      Do ISpin = 1, TopSpin
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         Call GA_Add_Patch(1.0d0, g_W2ij(ISpin), 1, C+A, 1, C+A,
     $                     1.0d0, g_W2(ISpin),   1, C+A, 1, C+A,
     $                            g_W2(ISpin),   1, C+A, 1, C+A)
      EndDo
C
C     Print results if requested
C
      If ( Util_Print('partial w2', Print_Debug) ) then
         Write (LuOut, *) 'W2(ij) contributions from Fock build'
         Do ISpin = 1, TopSpin
            Call GA_Print( g_W2ij(ISpin) )
            Call GA_Print( g_W2(ISpin) )
         EndDo
      EndIf
C
      Do ISpin = 1, MaxSpin
         If ( .NOT. GA_Destroy( g_W2ij(ISpin) ) ) Call ErrQuit(
     $        'RIMP2_Mk_W2: can''t destroy L_tr', ISpin, GA_ERR)
      EndDo
C
C     *************************
C     * P2 contribution to W2 *
C     *************************
C     W2(mn) <-- -P2(mn) e(m)   ! CHEATING! derivation shows e(n)
C     W2(an) <-- -P2(an) e(n)    
C     W2(ab) <-- -P2(ab) e(a)   ! CHEATING! derivation shows e(b)
C
      Do ISpin = 1, TopSpin
C
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
C        For convenience, copy P2 as workspace
C
         If ( .NOT. GA_Duplicate(g_P2(ISpin), g_Tmp, 'Temp for W2') )
     $        Call ErrQuit('RIMP2_Mk_W2: can''t allocate copy of P2',
     $        ISpin, GA_ERR)
         Call GA_Copy(g_P2(ISpin), g_Tmp)
C
C        Find our patch of Tmp and access it
C
         Call GA_Distribution(g_Tmp, Me, PLo, PHi, QLo, QHi)
         If ( (PHi-PLo+1)*(QHi-QLo+1) .gt. 0) then
            Call GA_Access(g_Tmp, PLo, PHi, QLo, QHi, I_Tmp, LD)
C
C           Scale each row by the corresponding eigenvalue
C
            Do P = PLo, PHi
               Call DScal(QHi-QLo+1, -Eig(P, ISpin), Dbl_MB(I_Tmp+P-1),
     $              LD)
            EndDo               ! P
C            
            Call GA_Release_Update(g_Tmp, PLo, PHi, QLo, QHi)
         EndIf
         Call GA_Sync
c$$$         Write (6, *) 'Scaled 1-PDM'
c$$$         Call GA_Print(g_Tmp)
C
C        Since we've scaled the whole of P2, we need to patch up
C        a couple of things.
C
C        W2(an) <-- W2(na) ! CHEATING! this is transpose of what it should be!
C
         Call GA_Copy_Patch('T', g_Tmp, 1,     C+A,   C+A+1, C+A+V,
     $                           g_Tmp, C+A+1, C+A+V, 1,     C+A)
     $                           
         Call GA_Copy_Patch('T', g_Tmp, C+A+1, C+A+V, 1,     C+A,
     $                           g_Tmp, 1,     C+A,   C+A+1, C+A+V )
c$$$C
c$$$C        W2(iy) <-- W2(yi) ! CHEATING! this is transpose of what it should be
c$$$C
c$$$         If ( C .gt. 0) Call GA_Copy_Patch('T', g_Tmp, 1, C, C+1, C+A,
c$$$     $                                          g_Tmp, C+1, C+A, 1, C)
c$$$     $                                          
c$$$         If ( C .gt. 0) Call GA_Copy_Patch('T', g_Tmp, C+1, C+A, 1, C,
c$$$     $                                          g_Tmp, 1, C, C+1, C+A)
C
c$$$         Write (6, *) 'Scaled 1-PDM after patching'
c$$$         Call GA_Print(g_Tmp)
c$$$C
C        Now add into W2
C
         Call GA_Add(1.0d0, g_Tmp, 1.0d0, g_W2(ISpin), g_W2(ISpin) )
C
         If ( .NOT. GA_Destroy(g_Tmp) ) Call ErrQuit(
     $        'RIMP2_Mk_W2: can''t destroy temporary GA', g_Tmp, GA_ERR)
C
C        Print results if requested
C
         If ( Util_Print('partial w2', Print_Debug) ) then
            Write (LuOut, *) 'After RIMP2_Mk_W2...'
            Call GA_Print( g_W2(ISpin) )
         EndIf
C
      EndDo
C
      Return
      End
