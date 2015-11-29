      LOGICAL FUNCTION WFN1_AD_DLAISNAN( DIN1, DIN2 )
      USE WFN1_AD1
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2010
*
*     .. Scalar Arguments ..
      TYPE(WFN1_AD_DBLE) :: DIN1, DIN2
*     ..
*
*  Purpose
*  =======
*
*  This routine is not for general use.  It exists solely to avoid
*  over-optimization in DISNAN.
*
*  DLAISNAN checks for NaNs by comparing its two arguments for
*  inequality.  NaN is the only floating-point value where NaN != NaN
*  returns .TRUE.  To check for NaNs, pass the same variable as both
*  arguments.
*
*  A compiler must assume that the two arguments are
*  not the same variable, and the test will not be optimized away.
*  Interprocedural or whole-program optimization may delete this
*  test.  The ISNAN functions will be replaced by the correct
*  Fortran 03 intrinsic once the intrinsic is widely available.
*
*  Arguments
*  =========
*
*  DIN1    (input) DOUBLE PRECISION
*
*  DIN2    (input) DOUBLE PRECISION
*          Two numbers to compare for inequality.
*
*  =====================================================================
*
*  .. Executable Statements ..
      WFN1_AD_DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
