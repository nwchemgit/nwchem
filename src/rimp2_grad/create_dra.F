C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_create_dra -- Create DRAs to previously computed sizes
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Create_DRA(MaxSpin, TopSpin, NFit, BraKetSame, 
     $     NameI, NameG, IntRange, GamRange, D_Int, D_Gam)
      Implicit NONE
#include "errquit.fh"
C
      Integer MaxSpin ![in]
      Integer TopSpin ![in]
      Integer NFit ![in]
      Logical BraKetSame ![in]
      Character*(*) NameI(MaxSpin, 2, 2) ![in]
      Character*(*) NameG(TopSpin) ![in]
      Integer IntRange(2, 2, MaxSpin, 2, 2) ![in]
      Integer GamRange(2, 2, MaxSpin) ![in]
      Integer D_Int(MaxSpin, 2, 2) ![out]
      Integer D_Gam(TopSpin) ![out]
C
C DESCRIPTION
C     Create DRAs according to names and sizes determined in
C     RIMP2_Prep_DRA.
C
C NOTES
C     The efficacy of this assumes there is no per-DRA memory
C     penalty.  This is true as of 1999-05-27.
C
C EXTERNAL ROUTINES
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "inp.fh"
#include "dra.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Spin, J, K, N1, N2
C
      If ( TopSpin .gt. MaxSpin ) Call ErrQuit(
     $     'RIMP2_Create_DRA: fatal program error: TopSpin > MaxSpin',
     $     MaxSpin, INPUT_ERR)
C
C     Initialize DRA handles to 0
C
      Call IFill(MaxSpin*2*2, 0, D_Int(1, 1, 1), 1)
      Call IFill(MaxSpin,     0, D_Gam(1),       1)
C
C     Loop over spins for both integrals and gammas
C
      Do Spin = 1, TopSpin
C
C        For integrals loop over bra/ket and class too
C
         Do J = 1, 2            ! Bra/Ket
            Do K = 1, 2         ! Class
               N1 = IntRange(2,1,Spin,J,K)-IntRange(1,1,Spin,J,K)+1
               N2 = IntRange(2,2,Spin,J,K)-IntRange(1,2,Spin,J,K)+1
               If ( Dble(N1)*Dble(N2)*Dble(NFit) .gt. 0) then
C
                  If ( DRA_Create(MT_Dbl, N1*N2, NFit,
     $                 'MO 3-center integrals',
     $                 NameI(Spin, J, K)(:Inp_StrLen(NameI(Spin,
     $                 J, K))), DRA_RW, N1, NFit, D_Int(Spin, J, K) )
C
     $                 .ne. 0) Call ErrQuit(
     $                 'RIMP2_Create_DRA: failed for int spin', Spin,
     &       INPUT_ERR)
               EndIf
            EndDo               ! Class (K)
         EndDo                  ! Bra/Ket (J)
C
C        For later convenience, we duplicate the handles when the
C        integrals are symmetric.
C
         If ( BraKetSame ) then
            D_Int(Spin, 2, 1) = D_Int(Spin, 1, 1)
         EndIf
C
C        Handle Gammas too, which have only spin cases
C
         N1 = GamRange(2, 1, Spin) - GamRange(1, 1, Spin) + 1
         N2 = GamRange(2, 2, Spin) - GamRange(1, 2, Spin) + 1
         If ( Dble(N1)*Dble(N2)*Dble(NFit) .gt. 0.0) then
C     
            If ( DRA_Create(MT_Dbl, N1*N2, NFit, 'MO 3-center gamma',
     $           NameG(Spin)(:Inp_StrLen(NameG(Spin))),
     $           DRA_RW, N1, NFit, D_Gam(Spin) )
C
     $           .ne. 0) Call ErrQuit(
     $           'RIMP2_Create_DRA: failed for gamma3 spin', Spin,
     &       INPUT_ERR)
         EndIf
      EndDo                     ! Spin
C
      Return
      End
