!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! NAME
!     GEAXPY -- Matrix Y <-- a*X + Y
!
! REVISION
!     $Id$
!
! SYNOPSIS
      Subroutine GEAXPY(M, N, Alpha, X, LDX, Y, LDY)
      Implicit NONE
      Integer M, N, LDX, LDY
      Double precision Alpha
      Double precision X(LDX, N), Y(LDY, N)
!
! ARGUMENTS
!     M       Row dimension of X, Y [IN]
!     N       Column dimension of X, Y [IN]
!     Alpha   Scale factor for X [IN]
!     X       MxN matrix [IN]
!     LDX     Leading dimension of X [IN]
!     Y       Result matrix [INOUT]
!     LDY     Leading dimension of Y [IN]
!
! DESCRIPTION
!     A simple convenience routine to simplify the addition of two
!     matrices which may have different leading dimensions (which makes
!     use of DAXPY unsafe).
!
!     Implemented as DAXPY calls on each column of X/Y.  Easy.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! LOCAL VARIABLES
      Integer J
!
      Do J = 1, N
         Call daxpy(M, Alpha, X(1, J), 1, Y(1, J), 1)
      EndDo
!
      Return
      End
