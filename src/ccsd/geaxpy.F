!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! NAME
!     GEAXPY -- Matrix Y <-- a*X + Y
!
! REVISION
!     $Id: geaxpy.f 28747 2016-11-15 22:44:47Z jhammond $
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
      Integer J, K
!
#if USE_BLAS_NOT_LOOPS
      Do J = 1, N
         Call daxpy(M, Alpha, X(1, J), 1, Y(1, J), 1)
      EndDo
#else
      do J = 1, N
         do K = 1, M
            Y(K, J) = Y(K, J) + Alpha * X(K, J)
         enddo
      enddo
#endif
!
      Return
      End
