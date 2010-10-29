C***********************************************************************
C NAME
C     Mk_1d_shell_dist -- compute mapping array for distibution
C     of an array blocked by shells.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function Mk_1d_Shell_Dist(Basis, NProc, MinChunk,
     $   Dist, NBlocks)
      Implicit NONE
      Integer Basis, NProc, MinChunk, NBlocks
      Integer Dist(NProc)
C
C ARGUMENTS
C     Basis     Handle for basis to be distributed [IN]
C     NProc     Number of elements in distribution [IN]
C     MinChunk  Minimum acceptable size for blocks of the distribution [IN]
C     Dist      Distribution specification, with NBlocks valid elements [OUT]
C     NBlocks   Number of blocks in computed distribution [OUT]
C
C RETURN VALUE
C     .TRUE.    Successful completion
C     .FALSE.   Error in basis set queries
C
C DESCRIPTION
C     Returns a distribution (suitable for use in GA_Create_Irreg) for
C     the specified basis in which shell blocks are kept together.
C
C     The current algorithm is simple and runs in linear time.  It may
C     not provide the optimal distribution.
C
C INCLUDED FILES
#include "bas.fh"
#include "stdio.fh"
C***********************************************************************
C LOCAL VARIABLES
      Logical Status
      Integer NBF, NSh, NAvg, N_on_Proc, Sh, Lo, Hi, NBFSh
C
C     First check out the basis
C
      Status = .TRUE.
      Status = Status .AND. Bas_NumBF( Basis, NBF)
      Status = Status .AND. Bas_NumCont( Basis, NSh)
C      
      If ( .NOT. Status ) then
         Write (LuErr, *) 'Mk_1d_Shell_Dist: ',
     $      'Unable to get required basis info.'
         Mk_1d_Shell_Dist = .FALSE.
      EndIf
C
C     How many basis functions do we need to put on each processor
C
      NAvg = NBF / NProc
C
C     Divide up the space so that each each chunk exceeds the
C     average (or minimum)
C
 1000 N_on_Proc = 0
      NBlocks = 1
      Dist(1) = 1
C
      Do Sh = 1, NSh
C
         If ( .NOT. Bas_Cn2BFR( Basis, Sh, Lo, Hi) ) then
            Write (LuErr, *) 'Mk_1d_Shell_Dist: ',
     $         'Unable to get required basis info.'
            Mk_1d_Shell_Dist = .FALSE.
            Return
         EndIf
C
C        Size of this shell block
C
         NBFSh = Hi-Lo+1
C
C        If we haven't yet got enough on this processor, include it here.
C        Otherwise, start a new block.
C
         If ( N_on_Proc .lt. Max( NAvg, MinChunk) ) then
            N_on_Proc = N_on_Proc + NBFSh
         Else
            NBlocks = NBlocks + 1
            If ( NBlocks .le. NProc) then
               Dist( NBlocks ) = Lo
               N_On_Proc = NBFSh
            Else
               Goto 8000
            EndIf
         EndIf
      EndDo
C
 8000 If (NBlocks .gt. NProc) then
         NAvg = NAvg + 1
         Goto 1000
      EndIf
C
      Mk_1d_Shell_Dist = .TRUE.
      Return
      End
