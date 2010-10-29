C***********************************************************************
C NAME
C     ga_create_one_by_shells -- Create a global array with a single
C     block in one dimension and shell blocking in the other.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function GA_Create_one_by_shells(Type, Dim1, Basis, Name,
     $   G_A)
      Implicit NONE
      Integer Type, Dim1, Basis, G_A
      Character*(*) Name
C
C ARGUMENTS
C     Dim1   Fist dimension of matrix (not distributed) [IN]
C     Basis  Handle to fitting basis (distributed) [IN]
C     Name   Name to be given to global array [IN]
C     G_A    Handle to global array [OUT]
C
C RETURN VALUE
C     .TRUE.    Successful completion
C     .FALSE.   Error in computing distirbution or allocation
C
C DESCRIPTION
C     Creates a global array which is distributed only in the second
C     dimension, and that one is blocked so as not to break up a shell
C     block in the specified Basis.
C
C INCLUDED FILES
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
C
C ROUTINES REQUIRED
      Logical Mk_1d_Shell_Dist
      External Mk_1d_Shell_Dist
C***********************************************************************
C CONSTANTS
      Integer MinChunk
      Parameter (MinChunk = 10)
C
C LOCAL VARIABLES
      Logical Status
      Integer NBlocks, Dim2, H_Dist, Dist_Idx
C
C     Compute the distribution for the fitting basis
C
      Status = .TRUE.
      Status = MA_Push_Get( MT_Int, GA_NNodes(), 'Distribution tmp',
     $  H_Dist, Dist_Idx)
      Status = Status .AND. Mk_1d_Shell_Dist( Basis, GA_NNodes(),
     $   MinChunk, Int_MB(Dist_Idx), NBlocks)
C
C     Also need the total size of the fitting basis
C
      Status = Status .AND. Bas_NumBF( Basis, Dim2)
C
C     Now create the array
C
      Status = Status .AND. GA_Create_Irreg( Type, Dim1, Dim2,
     $   Name, 1, 1, Int_MB(Dist_Idx), NBlocks, G_A)
C
      Status = Status .AND. MA_Pop_Stack( H_Dist)
C
      GA_Create_One_By_Shells = Status
      Return
      End


