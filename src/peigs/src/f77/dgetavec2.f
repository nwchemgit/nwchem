*
* $Id$
*
      SUBROUTINE DGETAVEC2( iii, LAMBDA, DELTA, N,
     $     B1, BN, L, D, LD,
     $     LLD, LPLUS,
     $     DPLUS, UMINUS, DMINUS, T, P, GAMMA, Z, K,
     $     ZTZ, ZBEGIN, ZEND, vecno, index)
*     
*  -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*
*     .. Scalar Arguments ..
      integer iii
      INTEGER            N, B1, BN, K, ZBEGIN, ZEND, vecno, index(*)
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
*  DGETAVEC computes the FP vector of the submatrix indexed from
*  B1 to BN of (L*D*L^T - (LAMBDA+DELTA)*I) using the qd algorithms.
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
c assume that vector is already there
c     
      ztz = dnrm2( bn-b1+1, z(b1), 1)
      write(*,*) ' in dgetavec2 ', ztz

c     
      CALL DDSTQDS( LAMBDA, DELTA, N, B1, BN, L, D, LD,
     $     LPLUS, DPLUS, T )
      CALL DDQDS( LAMBDA, DELTA, N, B1, BN, L, D, LLD,
     $     UMINUS, DMINUS, P )
c
      CALL DCMPGAMMA( LAMBDA, DELTA, N, B1, BN, L, D,
     $     LD, LLD,
     $     LPLUS, DPLUS, UMINUS, DMINUS, T, P, K, GAMMA )
c     
c     DO I = B1, bn
c     write(*,*) ' gamma ', i, gamma(i), k
c     enddo
c     
      z(k) = 1.0d0
      if ( iii .eq. 0 ) then
         z(k) = 1.0
      endif
c      if ( iii .gt. 0 ) then
c         k = 0
c         do i = b1, bn
c            index(i) = i
c            gamma(i) = abs(gamma(i))
c         enddo
c         j = bn-b1+1
c         call dshellsort2(j, gamma(b1), index(b1))
c         DO I = B1, bn
c            write(*,*)  ' after sorting gamma ', i,
c     $           gamma(i), index(i)
c         enddo
c         do i = b1, bn
c            if ( gamma(i) .lt. 10.*dlamch('u') ) then
c               j = index(i)
c               Z(j) = ONE
c               k = max(k, j)
c               write(*,*) 'fill  z(j) = 1 ', j
c            endif
c         enddo
c      endif
c     
c     should probably lower threshold this
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
      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*
*     End of DGETAVEC
*
      END




