C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Delete_DRA -- Close and delete DRAs
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Delete_DRA(MaxSpin, TopSpin, BraKetSame,
     $     D_Int, D_Gam)
      Implicit NONE
#include "errquit.fh"
C
      Integer MaxSpin ![in]
      Integer TopSpin ![in]
      Logical BraKetSame ![in]
      Integer D_Int(MaxSpin, 2, 2) ![in]
      Integer D_Gam(TopSpin) ![in]
C
C DESCRIPTION
C     Close and delete all open DRAs.  Decides what is open by looking
C     at the handles.
C
C PRINT CONTROLS
C EXTERNAL ROUTINES
C INCLUDE FILES
#include "dra.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Spin, J, K
C
C     Sanity check
C
      If ( TopSpin .gt. MaxSpin ) Call ErrQuit(
     $     'RIMP2_Create_DRA: fatal program error: TopSpin > MaxSpin',
     $     MaxSpin, INPUT_ERR)
C
C     Loop over spins for both integrals and gammas
C
      Do Spin = 1, TopSpin
C
C        Undo our convenience from RIMP2_Create_DRA
C
         If ( BraKetSame ) D_Int(Spin, 2, 1) = 0
C
C        For integrals loop over bra/ket and class
C
         Do J = 1, 2            ! Bra/Ket
            Do K = 1, 2         ! Class
               If ( D_Int(Spin, J, K) .ne. 0 ) then
                  If ( DRA_Delete( D_Int(Spin, J, K) ) .ne. 0)
     $                 Call ErrQuit(
     $                 'RIMP2_Delete_DRA: failed for int spin', Spin,
     &       CALC_ERR)
               EndIf
            EndDo               ! Class (K)
         EndDo                  ! Bra/Ket (J)
C
C        Handle Gammas too, which have only spin cases
C
         If ( D_Gam(Spin) .ne. 0) then
            If ( DRA_Delete( D_Gam(Spin) ) .ne. 0) Call ErrQuit(
     $           'RIMP2_Delete_DRA: failed for gamma3 spin', Spin,
     &       CALC_ERR)
         EndIf
      EndDo                     ! Spin
C
      Return
      End
