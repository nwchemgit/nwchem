      FUNCTION WFN1_AD_DLAPY2( X, Y )
      USE WFN1_AD1
      IMPLICIT NONE
      TYPE(WFN1_AD_DBLE) :: WFN1_AD_DLAPY2
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      TYPE(WFN1_AD_DBLE) :: X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      TYPE(WFN1_AD_DBLE) :: W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
c     INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         WFN1_AD_DLAPY2 = W
      ELSE
         WFN1_AD_DLAPY2 = W*SQRT( ONE+( Z / W )**2.0d0 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
