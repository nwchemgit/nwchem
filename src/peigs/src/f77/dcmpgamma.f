*
* $Id$
*

      SUBROUTINE DCMPGAMMA( LAMBDA, DELTA, N, B1, BN, L, D,
     $                   LD, LLD, LPLUS, DPLUS, UMINUS, DMINUS, T,
     $                   P, K, GAMMA )
*
*  -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*
*     .. Scalar Arguments ..
	implicit none
      INTEGER            N, B1, BN, K
      DOUBLE PRECISION   DELTA, LAMBDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), L( * ), P( * ), GAMMA( * ),
     $                   DMINUS( * ), LPLUS( * ), T( * ), UMINUS( * ),
     $                   DPLUS( * ), LD( * ), LLD( * )
*     ..
*
*  Purpose
*  =======
*
*  DCMPGAMMA computes the GAMMA array, where GAMMA(I) is the
*  reciprocal of the I^{th} diagonal element of the inverse of
*  (L*D*L^T - (LAMBDA+DELTA)*I)
*
*  Arguments
*  =========
*
*  LAMBDA  (input) DOUBLE PRECISION
*          The shift.
*
*  DELTA   (input) DOUBLE PRECISION
*          Lower order bits of the shift.
*
*  N       (input) INTEGER
*          The order of the matrix L * D * L^T.
*
*  B1      (input) INTEGER
*          Starting index of the submatrix (of L * D * L^T).
*
*  BN      (input) INTEGER
*          Last index of the submatrix (of L * D * L^T).
*
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of the bidiagonal matrix
*          L, in elements 1 to N-1.  L(N) need not be set.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the diagonal matrix D.
*
*  LD      (input) DOUBLE PRECISION array, dimension (N-1)
*          The n-1 elements L(i)*D(i).
*
*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)
*          The n-1 elements L(i)*L(i)*D(i).
*
*  LPLUS   (output) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) diagonal elements of L+.
*
*  DPLUS   (output) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of D+.
*
*  UMINUS  (output) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) diagonal elements of U-.
*
*  DMINUS  (output) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of D-.
*
*  T       (output) DOUBLE PRECISION array, dimension (N)
*          Intermediate results of the dstqds algorithm.
*
*  P       (output) DOUBLE PRECISION array, dimension (N)
*          Intermediate results of the dqds algorithm.
*
*  K       (output) INTEGER
*          The k^{th} column of the inverse of (L*D*L^T - (LAMBDA+DELTA)*I).
*
*  GAMMA   (output) DOUBLE PRECISION array, dimension (N)
*          GAMMA(i) is the reciprocal of the i^{th} diagonal element
*          of the inverse of (L*D*L^T - (LAMBDA+DELTA)*I).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d+0, ONE = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   DIF, EPS, MINDIF
*     ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
      integer            doprnt1,doprnt2
      common             doprnt1,doprnt2

*     ..
*     .. Executable Statements ..
*

      EPS = 0.111022302462515654E-15
      K = B1
      IF( B1.EQ.1 ) THEN
         MINDIF = P( B1 ) - DELTA
      ELSE
         MINDIF = ( LLD( B1-1 ) + P( B1 ) ) - DELTA
      END IF
      IF( MINDIF.EQ.ZERO ) THEN
         MINDIF = EPS * P( B1 )
      END IF
      GAMMA( B1 ) = MINDIF
      DIF = ( D( BN ) + T( BN ) ) - DELTA
      IF( DIF.EQ.ZERO ) THEN
         DIF = EPS * T( BN )
      END IF
      GAMMA( BN ) = DIF
      IF( ABS( DIF ).LT.ABS( MINDIF ) ) THEN
          MINDIF = DIF
          K = BN
      END IF

      DO I = B1+1, BN-1
         DIF = ( ( P( I ) + T( I ) ) + LAMBDA ) - DELTA
         IF( DIF.EQ.ZERO ) THEN
            DIF = EPS * P( I )
         END IF
         GAMMA( I ) = DIF
         IF( ABS( DIF ).LT.ABS( MINDIF ) ) THEN
             MINDIF = DIF
             K = I
         END IF
      END DO


      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*
*     End of DCMPGAMMA
*
      END
