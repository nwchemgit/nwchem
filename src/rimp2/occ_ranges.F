C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     occ_ranges -- Convert an occupation vector into ranges of occupied
C     and virtual orbitals
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Occ_Ranges( N, Occ, OLo, OHi, VLo, VHi)
      Implicit NONE
      Integer N, OLo, OHi, VLo, VHi
      Double precision Occ(N)
C
C ARGUMENTS
C     N    Length of occupation vector [IN]
C     Occ  Occupation vector [IN]
C     OLo  Lower limit of occupied range [OUT]
C     OHi  Upper limit of occupied range [OUT]
C     VLo  Lower limit of virtual range [OUT]
C     VHi  Upper limit of virtual range [OUT]
C
C DESCRIPTION
C     Convert a vector of occupation numbers into ranges of occupied
C     and virtual orbitals.  If either range is empty, it will be
C     returned as 1:0, which can be fed to a Fortran DO loop and will 
C     NOT result in an iteration.
C
C     This routine is presently very simplistic, but should suffice 
C     for most purposes.  Occupied and virtual blocks are assumed to be
C     contiguous, with occupied first and virtuals second.
C
C     Any orbital with non-zero occupation is taken as occupied.
C
C INCLUDE FILES
#include "numerical_constants.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C
      Integer I
C
C     Initialize things in a way that Fortran loops will work if there
C     are no occupied or not virtual orbitals.
C
      OLo = 1
      OHi = 0
      VLo = 1
      VHi = 0
C
C     Assume that all occupied are contiguous at the bottom of the
C     array and all the virtuals are contiguous at the top.  So all
C     we have to do is find the break point.
C
      Do I = 1, N
         If ( Occ( I ) .eq. FP_Zero) then
            OHi = I - 1
            VLo = I
            VHi = N
            Return
         EndIf
      EndDo
C
      Return
      End
