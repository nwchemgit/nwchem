      SUBROUTINE WFN1_AD_DLASCALE( N, X, INCX, SCALE )
      USE WFN1_AD1
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      TYPE(WFN1_AD_DBLE) :: SCALE
*     ..
*     .. Array Arguments ..
      TYPE(WFN1_AD_DBLE) :: X( * )
*     ..
*
*  Purpose
*  =======
*
*  Companion routine of DLASSQ. This routine establishes SCALE as
*     scl = max( scale, abs( x( i ) ) )
*  and set SCALE = scl before returning.
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Unfortunately the approach outlined above can lead to severe numerical
*  errors using automatic differentiation to calculate the gradients.
*  These errors can be avoided but only at the cost of a second pass
*  through the vector. This problem is equivalent to the one in 
*  wfn1_ad_dnrm2, see the comments in that routine for details.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      TYPE(WFN1_AD_DBLE) :: ABSXI, scl, smsq
*     ..
*     .. Intrinsic Functions ..
c     INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
*
*        This additional loop fixes the numerical problems with the
*        gradient evaluation.
*
         do ix = 1, 1+(n-1)*incx, incx
           scale = max(scale,abs(x(ix)))
         enddo
      END IF

      RETURN
*
*     End of DLASCALE
*
      END
