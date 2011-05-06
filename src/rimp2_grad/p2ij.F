C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_P2ij -- Form P^{(2)}_{ij} for batch in memory
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_P2ij(SameSpin, IJSym, IJTri, ILo, IHi, NVirI,
     $     JLo, JHi, NVirJ, Mixed, G_T, G_P2, Occ1)
      Implicit NONE
C
      Logical SameSpin
      Logical IJSym
      Logical IJTri
      Integer ILo
      Integer IHi
      Integer NVirI
      Integer JLo
      Integer JHi
      Integer NVirJ
      Logical Mixed
      Integer G_T
      Integer G_P2
      Integer Occ1
C
C DESCRIPTION
C     Evaluate P^{(2)}_{ij} = -1/2 Sum(abk) t(ik,ab) t(jk,ab)
C
C     Non-zero spin cases:
C     i j   k a b
C     a a   a a a   t(IK,AB) t(JK,AB)
C
C     a a   b b a   t(Ik,aB) t(Jk,aB) --> t(Ik,Ba) t(Jk,Ba)
C     a a   b a b   t(Ik,Ab) t(Jk,Ab)
C
C     b b   a a b   t(iK,Ab) t(jK,Ab) --> t(Ki,Ab) t(Kj,Ab)
C     b b   a b a   t(iK,aB) t(jK,aB) --> t(Ki,Ba) t(Kj,Ba)
C
C     b b   b b b   t(ik,ab) t(jk,ab)
C
C     So we learn
C     1) Pure spin amplitudes contribute to matching P spin
C     2) Mixed-spin amplitudes contribute to both P spins
C     3) Mixed-spin contributions carry factor of 2*2, first
C        due to two contributing spin cases (above), second account
C        for the fact that t(iJ,aB) is numerically identical to 
C        t(Ij,Ab) and only latter appears in this code.
C
C PRINT CONTROLS
C     'partial p2' will print the G_P2 array at the end of the routine.
C
C INCLUDED FILES
#include "util.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, J, JTop, LT1, LT2, NI, NJ, Type, TrueI, TrueJ, K
      Double Precision Pij, Fact
      Double Precision GA_DDot_Patch
      external GA_DDot_Patch
      Integer ga_nodeid
      external ga_nodeid
C
C     Determine the size of the matrix.  Dimensions must be multiples
C     of NVirI and NVirJ respectively.
C
      Call GA_Inquire( G_T, Type, LT1, LT2)
      NI = LT1 / NVirI
      NJ = LT2 / NVirJ
C
      Call GA_Sync
C     
C     Initialize the true I and J indices based on the caller's info
C
      TrueI = ILo
      TrueJ = JLo
C
C     To understand the factor, see the discussion of spin cases, above.
C
      Fact = -1.0d0
      If ( SameSpin) Fact = -0.5d0
c$$$      If ( .NOT. IJSym ) Fact = -1.0d0
C
      If ( Mixed ) then
         Do I = 1, NJ
            Do J = 1, NJ
               Pij = 0.0d0
C
               Do K = 1, NI
                  Pij = Pij + GA_DDot_Patch(
     $               G_T, 'N', (K-1)*NVirI+1, K*NVirI,
     $                         (I-1)*NVirJ+1, I*NVirJ,
     $               G_T, 'N', (K-1)*NVirI+1, K*NVirI,
     $                         (J-1)*NVirJ+1, J*NVirJ )
               EndDo
               If ( Pij .ne. 0.0d0 .and. ga_nodeid().eq.0 ) then
                  Call GA_Acc( G_P2, Occ1+I-1, Occ1+I-1,
     $               Occ1+J-1, Occ1+J-1, Pij, 1, Fact)
               EndIf
            EndDo
         EndDo
      Else                      ! .NOT. Mixed
         Do I = 1, NI
            Do J = 1, NI
               Pij = 0.0d0
C
               Do K = 1, NJ
                  Pij = Pij + GA_DDot_Patch(
     $               G_T, 'N', (I-1)*NVirI+1, I*NVirI,
     $                         (K-1)*NVirJ+1, K*NVirJ,
     $               G_T, 'N', (J-1)*NVirI+1, J*NVirI,
     $                         (K-1)*NVirJ+1, K*NVirJ )
               EndDo
               If ( Pij .ne. 0.0d0 .and. ga_nodeid().eq.0 ) then
                  Call GA_Acc( G_P2, Occ1+I-1, Occ1+I-1,
     $               Occ1+J-1, Occ1+J-1, Pij, 1, Fact)
               EndIf
            EndDo
         EndDo
      EndIf
      Call GA_Sync
C
      If ( Util_Print('partial p2', Print_Debug)) then
         Write (LuOut, *) 'After RIMP2_P2ij...'
         Call GA_Print(g_P2)
      EndIf
C
      Return
      End
