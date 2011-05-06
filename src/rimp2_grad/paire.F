C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_PairE -- Calculate MP2 pair energies for IJ batch in memory
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_PairE(ILo, IHi, NVirI, JLo, JHi, NVirJ, IJSym,
     $   IJTri, G_Int, G_T, G_PairE)
      Implicit NONE
C
      Integer ILo      ! [in] Lowest occ index in batch (row)
      Integer IHi      ! [in] Highest occ index in batch (row)
      Integer NVirI    ! [in] Size of blocks (row)
      Integer JLo      ! [in] Lowest occ index in batch (col)
      Integer JHi      ! [in] Highest occ index in batch (col)
      Integer NVirJ    ! [in] Size of blocks (col)
      Logical IJSym    ! [in] Use I <--> J permutational symmetry if possible
      Logical IJTri    ! [in] G_I/G_T hold only triangular I > J
      Integer G_Int    ! [in] GA handle for integral batch
      Integer G_T      ! [in] GA handle for amplitude batch
      Integer G_PairE  ! [in] GA handle for pair energy array
C
C DESCRIPTION
C     Given a batch of integrals and amplitudes, compute pair energies
C     by dotting them together appropriately.  This is done in a
C     completely data-parallel fashion (G_Int and G_T must have the 
C     same distributions), and the only communications is accumulation
C     of pair energy contributions into G_PairE.
C
C     Integrals & amplitudes are expected to be organized as blocks of
C     all AB for a particular IJ.  A batch will contain blocks for
C     multiple Is (row dimension) and Js (column dimension) laid out
C     in rectangular fashion.  IJTri indicates that the minimum (I > J)
C     triangle is provided instead of the full rectangular matrix.
C
C     Where indicated (by IJSym and the fact that the requested ILo:IHi
C     and JLo:JHi ranges match up), only the I > J triangle of the
C     batch will be computed.
C
C PRINT CONTROLS
C     'partial pair energies' will print the G_PairE array at the end
C     of the routine.
C
C INCLUDED FILES
#include "util.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, J, JTop, LT1, LT2, NI, NJ, Type, TrueI, TrueJ
      Double Precision PairE
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
C     Use I <--> J symmetry if requested
C
      Do I = 1, NI
         If (IJSym .AND. ILo .eq. JLo .AND. IHi. eq. JHi
     $      .AND. .not. IJTri) then
            JTop = I
         Else
            JTop = NJ
         EndIf
         Do J = 1, JTop
C           
            PairE = GA_DDot_Patch(
     $         G_Int, 'N', (I-1)*NVirI+1, I*NVirI,
     $                     (J-1)*NVirJ+1, J*NVirJ,
     $         G_T,   'N', (I-1)*NVirI+1, I*NVirI,
     $                     (J-1)*NVirJ+1, J*NVirJ )
C
C           If this turns out to be too much traffic, we
C           might want to create a local array to contain
C           only our IJs and accumulate into it.  But
C           wait until perfomance says we need to.  For
C           now, we'll just avoid sending zeros across
C           the machine.
C
            If ( PairE .ne. 0.0d0 .and. ga_nodeid().eq.0 ) then
               Call GA_Acc( G_PairE, TrueI, TrueI, TrueJ, TrueJ,
     $            PairE, 1, 1.0d0)
            EndIf
C
C           Increment indices for the true I & J 
C
            TrueJ = TrueJ + 1
            If ( IJTri .AND. TrueJ .gt. TrueI ) then
               TrueI = TrueI + 1
               TrueJ = JLo
            EndIf
         EndDo
         If ( .NOT. IJTri ) then
            TrueI = TrueI + 1
            TrueJ = JLo
         EndIf
      EndDo
      Call GA_Sync
C
C     Print pair energies if requested (not necessarily complete yet)
C
      If ( Util_Print('partial pair energies', Print_Debug))
     $   Call GA_Print(g_PairE)
C
      Return
      End
