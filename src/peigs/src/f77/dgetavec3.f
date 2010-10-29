*
* $Id$
*
      SUBROUTINE DGETAVEC3( LAMBDA, DELTA, N, B1, BN, L, D, LD,
     $     LLD, LPLUS,
     $     DPLUS, UMINUS, DMINUS, T, P, GAMMA, Z, K,
     $     ZTZ, ZBEGIN, ZEND, vecno, index)
*     
*  -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*
*     .. Scalar Arguments ..
      INTEGER            N, B1, BN, K, ZBEGIN, ZEND,
     $     vecno, index(*)
      DOUBLE PRECISION   DELTA, LAMBDA, ZTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ),
     $                   P( * ), GAMMA( * ),
     $                   DMINUS( * ), LPLUS( * ),
     $                   T( * ), UMINUS( * ), DPLUS( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETAVEC3 computes eigenvectors using inverse iteration
*	by factoring up and down.
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
*  GAMMA   (output) DOUBLE PRECISION array, dimension (N)
*          GAMMA(i) is the reciprocal of the i^{th} diagonal element
*          of the inverse of (L*D*L^T - (LAMBDA+DELTA)*I).
*
*  Z       (output) DOUBLE PRECISION array, dimension (N)
*          The FP vector. Z(k) is returned to be 1.
*
*  K       (output) INTEGER
*          The k^{th} column of the inverse of (L*D*L^T - (LAMBDA+DELTA)*I).
*
*  ZTZ     (output) DOUBLE PRECISION
*          The square of the norm of the FP vector.
*
*  ZBEGIN  (output) INTEGER
*          For i < ZBEGIN, Z(i) < EPS. ZBEGIN >= B1.
*
*  ZEND    (output) INTEGER
*          For i > ZEND, Z(i) < EPS. ZEND <= BN.
*
*  vecno   (input) INTEGER --- ( used for printing purposes )
*          Index of the eigenvector desired.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d+0, ONE = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
c
      double precision dnrm2
      external dnrm2
c      
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT, ABS
      
      integer            doprnt1,doprnt2
      common             doprnt1,doprnt2
*     ..
*     .. Executable Statements ..
*     
      ZBEGIN = B1
      ZEND = BN
c
c     inverse iterate
c     
c     
c     
      CALL DDSTQDS( LAMBDA, DELTA, N, B1, BN, L, D, LD,
     $     LPLUS, DPLUS, T )
      CALL DDQDS( LAMBDA, DELTA, N, B1, BN, L, D, LLD,
     $     UMINUS, DMINUS, P )
c     
c     CALL DCMPGAMMA( LAMBDA, DELTA, N, B1, BN, L, D, LD, LLD,
c     $     LPLUS, DPLUS, UMINUS, DMINUS, T, P, K, GAMMA )
c     
c     
c     should probably lower threshold this
c     
      k = bn
c     
      DO I = K-1, B1, -1
         IF( Z( I+1 ).NE.ZERO ) THEN
            Z( I ) = - ( LPLUS( I ) * Z( I+1 ) )
         ELSE
            Z( I ) = - ( LD( I+1 ) / LD( I ) ) * Z( I+2 )
         END IF
      END DO
      DO I = K, BN-1
         IF( Z( I ).NE.ZERO ) THEN
            Z( I+1 ) = - ( UMINUS( I ) * Z( I ) )
         ELSE
            Z( I+1 ) = - ( LD( I-1 ) / LD( I ) ) * Z( I-1 )
         END IF
      END DO
c     
c     
c     ZTZ = ZERO
c     DO I = B1, BN
c     ZTZ = ZTZ + Z( I ) * Z( I )
c     END DO
c     
      ztz = dnrm2( bn-b1+1, z(b1), 1)
c     
      return
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*
*     End of DGETAVEC
*
      END






