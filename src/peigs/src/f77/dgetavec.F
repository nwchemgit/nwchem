*
* $Id$
*
#include "blas_lapackf.h"
c      
      SUBROUTINE DDSTQDS( LAMBDA, DELTA, N, B1, BN, L, D,
     $     LD, LPLUS, DPLUS, T )
*      
*     -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*      
*     .. Scalar Arguments ..
      implicit none
      INTEGER            N, B1, BN
      DOUBLE PRECISION   DELTA, LAMBDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), L( * ), LPLUS( * ),
     $                   LD( * ), T(*), DPLUS( * )
*     ..
*      
*  Purpose
*  =======
*      
*     DDSTQDS computes the L+ D+ L+^T decomposition of the submatrix
*     indexed from B1 to BN of (L*D*L^T - (LAMBDA+DELTA)*I)
*     by the differential form of the stationary qd algorithm.
*      
*     Arguments
*     =========
*      
*     LAMBDA  (input) DOUBLE PRECISION
*     The shift.
*      
*     DELTA   (input) DOUBLE PRECISION
*     Lower order bits of the shift.
*      
*     N       (input) INTEGER
*     The order of the matrix L * D * L^T.
*      
*     B1      (input) INTEGER
*     Starting index of the submatrix (of L * D * L^T).
*      
*     BN      (input) INTEGER
*     Last index of the submatrix (of L * D * L^T).
*      
*     L       (input) DOUBLE PRECISION array, dimension (N-1)
*     The (n-1) subdiagonal elements of the bidiagonal matrix
*     L, in elements 1 to N-1.  L(N) need not be set.
*      
*     D       (input) DOUBLE PRECISION array, dimension (N)
*     The n diagonal elements of the diagonal matrix D.
*      
*     LD      (input) DOUBLE PRECISION array, dimension (N-1)
*     The n-1 elements L(i)*D(i).
*      
*     LPLUS   (output) DOUBLE PRECISION array, dimension (N-1)
*     The (n-1) subdiagonal elements of L+.
*      
*     DPLUS   (output) DOUBLE PRECISION array, dimension (N)
*     The n diagonal elements of D+.
*      
*     T       (output) DOUBLE PRECISION array, dimension (N)
*     Intermediate results of the algorithm.
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
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
#include "blas_lapack.data"
*     ..
*     .. Executable Statements ..
*      


 310  CONTINUE
      IF( B1.EQ.1 ) THEN
         T( B1 ) = -LAMBDA
      ELSE
         T( B1 ) = LD( B1-1 ) * L( B1-1 ) - LAMBDA
      END IF
      
      DO I = B1, BN-1
         DPLUS( I ) = ( D( I ) + T( I ) ) - DELTA
c     LPLUS( I ) = LD( I ) / DPLUS( I )
         if ( dplus(i) .lt. DLAMCHU ) then
         LPLUS( I ) = LD( I ) /(DPLUS( I ) + DLAMCHU)
         else
         LPLUS( I ) = LD( I ) /DPLUS( I )
         endif

         T ( I+1 ) =
     $        ( T( I )-DELTA ) * LPLUS( I ) * L( I ) - LAMBDA
      END DO
      DPLUS( BN	) = ( D( BN ) + T( BN ) ) - DELTA
      
      IF(.NOT.( ( DPLUS( BN ).GT.ZERO) .OR.
     $     ( DPLUS( BN ).LT.ONE ) ) ) THEN
         IF( B1.EQ.1 ) THEN
            T( B1 ) = -LAMBDA
         ELSE
            T( B1 ) = LD( B1-1 ) * L( B1-1 ) - LAMBDA
         END IF
         
         DO I = B1, BN-1
            DPLUS( I ) = ( D( I ) + T( I ) ) - DELTA
         if ( dplus(i) .lt. DLAMCHU ) then
         LPLUS( I ) = LD( I ) /(DPLUS( I ) + DLAMCHU)
         else
         LPLUS( I ) = LD( I ) /DPLUS( I )
         endif
c            LPLUS( I ) = LD( I ) /DPLUS( I )

*     
*     Need to check the next few lines
*     
         IF(LPLUS( I ).EQ.0 
     $           .AND. ONE/(T( I )-DELTA).EQ.ZERO ) THEN
               T ( I+1 ) = L( I ) * LD( I ) - LAMBDA
            ELSE IF((T(I)-DELTA).EQ.ZERO
     $              .AND. ONE/LPLUS(I).EQ.0 ) THEN
               T ( I+1 ) = -ONE / (DPLUS( I )+ DLAMCHU)
            ELSE
               T ( I+1 ) =
     $              ( T( I )-DELTA )*LPLUS( I )*L( I ) - LAMBDA
            END IF
         END DO
         DPLUS( BN ) = ( D( BN ) + T( BN ) ) - DELTA
      END IF
       
      IF(.NOT.( ( DPLUS( BN ).GT.ZERO) .OR.
     $          ( DPLUS( BN ).LT.ONE ) ) ) THEN
         print *,"DDSTQDS : NaN detected!"
         GO TO 310
      END IF
       
       
      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*      
*     End of DDSTQDS
*      
      END
       
      SUBROUTINE DDQDS( LAMBDA, DELTA, N, B1, BN, L, D,
     $                  LLD, UMINUS, DMINUS, P )
*      
*  -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*      
*     .. Scalar Arguments ..
      implicit none
      INTEGER            N, B1, BN
      DOUBLE PRECISION   DELTA, LAMBDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( N ), L( N-1 ), P(N ), DMINUS( N ),
     $                   LLD( N-1 ), UMINUS( N )
*     ..
*      
*  Purpose
*  =======
*      
*  DDQDS computes the U- D- U-^T decomposition of the submatrix
*  indexed from B1 to BN of (L*D*L^T - (LAMBDA+DELTA)*I)
*  by the differential form of the progressive qd algorithm.
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
*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)
*          The n-1 elements L(i)*L(i)*D(i).
*      
*  UMINUS  (output) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) superdiagonal elements of U-.
*      
*  DMINUS  (output) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of D-.
*      
*  P       (output) DOUBLE PRECISION array, dimension (N)
*          Intermediate results of the algorithm.
*      
*  =====================================================================
*      
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d+0, ONE = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   tmp
*     ..
*     .. External Functions ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
#include "blas_lapack.data"
*     ..
*     .. Executable Statements ..
*      


 210  CONTINUE
      P( BN ) = D( BN ) - LAMBDA

      DO I = BN-1, B1, -1
         DMINUS( I+1 ) = ( LLD( I ) + P( I+1 ) ) - DELTA
         if ( dminus(i+1) .lt. DLAMCHU ) then
            TMP = D( I )/(DMINUS( I+1 )+DLAMCHS)
         else
            TMP = D( I ) / DMINUS( I+1 )
         endif
         UMINUS( I ) = L( I ) * TMP
         P ( I ) = ( P( I+1 ) - DELTA ) * TMP - LAMBDA
      END DO
      IF( B1.EQ.1 ) THEN
         DMINUS( B1 ) = P( B1 ) - DELTA
      ELSE
         DMINUS( B1 ) = ( LLD( B1-1 ) + P( B1 ) ) - DELTA
      END IF
       
      IF(.NOT.( ( DMINUS( B1 ).GT.ZERO) .OR.
     $          ( DMINUS( B1 ).LT.ONE ) ) ) THEN
         P( BN ) = D( BN ) - LAMBDA

         DO I = BN-1, B1, -1
            DMINUS( I+1 ) = ( LLD( I ) + P( I+1 ) ) - DELTA
            if ( dminus(i+1) .lt. DLAMCHS ) then
               TMP = D( I ) / (DMINUS( I+1 )+DLAMCHS)
            else
               TMP = D( I ) / DMINUS( I+1 )
            endif
            UMINUS( I ) = L( I ) * TMP
*      
*  Need to check the next few lines
*      
            IF( TMP.EQ.0 .AND. ONE/(P( I+1) - DELTA).EQ.ZERO ) THEN
               P ( I ) = ( L( I ) ) - LAMBDA
            ELSE IF((P( I+1)-DELTA).EQ.ZERO .AND. ONE/TMP.EQ.0) THEN
               P ( I ) = - ONE / DMINUS( I+1 )
            ELSE
               P ( I ) = ( P( I+1 ) - DELTA ) * TMP - LAMBDA
            END IF
         END DO
         IF( B1.EQ.1 ) THEN
            DMINUS( B1 ) = P( B1 ) - DELTA
         ELSE
            DMINUS( B1 ) = ( LLD( B1-1 ) + P( B1 ) ) - DELTA
         END IF
      END IF
       
      IF(.NOT.( ( DMINUS( B1 ).GT.ZERO) .OR.
     $     ( DMINUS( B1 ).LT.ONE ) ) ) THEN
         print *," peigs DDQDS : NaN detected!"
         GO TO 210
      END IF
      
       
      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*      
*     End of DDQDS
*      
      END
       
      SUBROUTINE DCMPGAMMA( LAMBDA, DELTA, N, B1,
     $     BN, L, D,
     $     LD, LLD, LPLUS, DPLUS, UMINUS,
     $     DMINUS, T,
     $     P, K, GAMMA )
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
      DOUBLE PRECISION   D( N ), L( N-1 ), P( N ), GAMMA( N ),
     $     DMINUS( N-1 ), LPLUS( N-1 ), T( N ),
     $     UMINUS( N-1 ),
     $     DPLUS( N ), LD( N-1 ), LLD( N-1 )
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
#include "blas_lapack.data"
       
*     ..
*     .. Executable Statements ..
*      


       
      EPS = DLAMCHP
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
*        DIF = P( I ) + ( (L( I-1 ) * L( I-1 )) * / DPLUS( I-1 ) ) * T(
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
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 11/ 3/97 at 12:05:10
C****************************************************************
C
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2
C   Version 6.1C1 on 10/14/97 at 11:50:45
C****************************************************************
C      
       
      SUBROUTINE DGETAVEC( iii, LAMBDA, DELTA, N, B1, BN, L, D, LD,
     $     LLD, LPLUS,
     $     DPLUS, UMINUS, DMINUS, T, P, GAMMA, Z, K,
     $     ZTZ, ZBEGIN, ZEND, vecno, index )
*      
*  -- LAPACK routine (version 0.0) -- in progress --
*     September 1995
*      
*     .. Scalar Arguments ..
      implicit none
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
      DOUBLE PRECISION   dummy, dummy1, dummy2
*     ..
*     .. External Functions ..
      INTEGER            idamax
      EXTERNAL           idamax
c      
      DOUBLE PRECISION dnrm2
      EXTERNAL dnrm2, dscal
c      
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT, ABS
#include "blas_lapack.data"
*     ..
*     .. Executable Statements ..
*      
      INTEGER I1, I2, kk
      INTEGER I3, I4
      ZBEGIN = B1
      ZEND = BN


      CALL DDSTQDS( LAMBDA, DELTA, N, B1, BN, L, D, LD,
     $     LPLUS, DPLUS, T )
      CALL DDQDS( LAMBDA, DELTA, N, B1, BN, L, D, LLD,
     $     UMINUS, DMINUS, P )
      CALL DCMPGAMMA( LAMBDA, DELTA, N, B1, BN, L, D, LD, LLD,
     $     LPLUS, DPLUS, UMINUS, DMINUS, T, P, K, GAMMA )
c      
c      
c      
c this is for an occasionally breakage on cray mpp
c can''t tell if it is the compiler or what
c
c      if ( abs(gamma(k)) .lt. DLAMCHS) then
c         Z( K ) = DLAMCHS
c      else	
c         Z( K ) = gamma(k)* ONE
c      endif
	z(k) = ONE
c     
c     DO I = B1, bn
c     write(*,*) ' gamma ', i, gamma(i), k
c     enddo
c      
c      if ( iii .gt. 0 ) then
c         goto 111
c         k = 0
c         do i = b1, bn
c            index(i) = i
c            gamma(i) = abs(gamma(i))
c         enddo
c         j = bn-b1+1
c         call dshellsort2(j, gamma(b1), index(b1))
c         DO I = B1, bn
c            write(*,*)  ' after sorting gamma ', i, gamma(i), index(i)
c         enddo
c         do i = b1, b1+iii
c            j = index(i)
c            Z(j) = ONE
c            k = max(k, j)
c            write(*,*) 'fill  z(j) = 1 ', j
c         enddo
c 111     continue
c         k = bn
c      else
c         Z( K ) = ONE
c      endif
c      write(*,*) 'final  k ', k, ' iii ', iii
c      
c      
      I2 = IAND(MAX0(K - B1,0),3)
      I4 = IAND(MAX0(I2,0),3)
      DO I = K - 1, K - I4, -1
         IF (Z(I+1) .NE. 0.0D+0) THEN
            Z(I) = -LPLUS(I)*Z(I+1)
         ELSE
            Z(I) = -(LD(I+1)/LD(I))*Z(I+2)
         ENDIF
      END DO

      DO I = K - I4 - 1, K - I2, -4
         kk = i-3 
         dummy = z(kk)
         dummy1 = lplus(kk)
         dummy2 = ld(kk)
         IF (Z(I+1) .NE. 0.0D+0) THEN
            Z(I) = -LPLUS(I)*Z(I+1)
         ELSE
            Z(I) = -(LD(I+1)/LD(I))*Z(I+2)
         ENDIF
         IF (Z(I) .NE. 0.0D+0) THEN
            Z(I-1) = -LPLUS(I-1)*Z(I)
         ELSE
            Z(I-1) = -(LD(I)/LD(I-1))*Z(I+1)
         ENDIF
         IF (Z(I-1) .NE. 0.0D+0) THEN
            Z(I-2) = -LPLUS(I-2)*Z(I-1)
         ELSE
            Z(I-2) = -(LD(I-1)/LD(I-2))*Z(I)
         ENDIF
         IF (Z(I-2) .NE. 0.0D+0) THEN
            Z(I-3) = -LPLUS(I-3)*Z(I-2)
         ELSE
            Z(I-3) = -(LD(I-2)/LD(I-3))*Z(I-1)
         ENDIF
      END DO

      DO I = K - I2 - 1, B1, -4
         kk = i-3 
         dummy = z(kk)
         dummy1 = lplus(kk)
         dummy2 = ld(kk)
         IF (Z(I+1) .NE. 0.0D+0) THEN
            Z(I) = -LPLUS(I)*Z(I+1)
         ELSE
            Z(I) = -(LD(I+1)/LD(I))*Z(I+2)
         ENDIF
         IF (Z(I) .NE. 0.0D+0) THEN
            Z(I-1) = -LPLUS(I-1)*Z(I)
         ELSE
            Z(I-1) = -(LD(I)/LD(I-1))*Z(I+1)
         ENDIF
         IF (Z(I-1) .NE. 0.0D+0) THEN
            Z(I-2) = -LPLUS(I-2)*Z(I-1)
         ELSE
            Z(I-2) = -(LD(I-1)/LD(I-2))*Z(I)
         ENDIF
         IF (Z(I-2) .NE. 0.0D+0) THEN
            Z(I-3) = -LPLUS(I-3)*Z(I-2)
         ELSE
            Z(I-3) = -(LD(I-2)/LD(I-3))*Z(I-1)
         ENDIF
      END DO
      I1 = IAND(MAX0(BN - K,0),3)
      I3 = IAND(MAX0(I1,0),3)
      DO I = K, I3 - 1 + K
         IF (Z(I) .NE. 0.0D+0) THEN
            Z(I+1) = -UMINUS(I)*Z(I)
         ELSE
            Z(I+1) = -(LD(I-1)/LD(I))*Z(I-1)
         ENDIF
      END DO
      DO I = I3 + K, I1 - 1 + K, 4
         IF (Z(I) .NE. 0.0D+0) THEN
            Z(I+1) = -UMINUS(I)*Z(I)
         ELSE
            Z(I+1) = -(LD(I-1)/LD(I))*Z(I-1)
         ENDIF
         IF (Z(I+1) .NE. 0.0D+0) THEN
            Z(I+2) = -UMINUS(I+1)*Z(I+1)
         ELSE
            Z(I+2) = -(LD(I)/LD(I+1))*Z(I)
         ENDIF
         IF (Z(I+2) .NE. 0.0D+0) THEN
            Z(I+3) = -UMINUS(I+2)*Z(I+2)
         ELSE
            Z(I+3) = -(LD(I+1)/LD(I+2))*Z(I+1)
         ENDIF
         IF (Z(I+3) .NE. 0.0D+0) THEN
            Z(I+4) = -UMINUS(I+3)*Z(I+3)
         ELSE
            Z(I+4) = -(LD(I+2)/LD(I+3))*Z(I+2)
         ENDIF
      END DO
      DO I = I1 + K, BN - 1, 4
         IF (Z(I) .NE. 0.0D+0) THEN
            Z(I+1) = -UMINUS(I)*Z(I)
         ELSE
            Z(I+1) = -(LD(I-1)/LD(I))*Z(I-1)
         ENDIF
         IF (Z(I+1) .NE. 0.0D+0) THEN
            Z(I+2) = -UMINUS(I+1)*Z(I+1)
         ELSE
            Z(I+2) = -(LD(I)/LD(I+1))*Z(I)
         ENDIF
         IF (Z(I+2) .NE. 0.0D+0) THEN
            Z(I+3) = -UMINUS(I+2)*Z(I+2)
         ELSE
            Z(I+3) = -(LD(I+1)/LD(I+2))*Z(I+1)
         ENDIF
         IF (Z(I+3) .NE. 0.0D+0) THEN
            Z(I+4) = -UMINUS(I+3)*Z(I+3)
         ELSE
            Z(I+4) = -(LD(I+2)/LD(I+3))*Z(I+2)
         ENDIF
      END DO
c      
c      
c     
      i = bn-b1+1
      ztz = dnrm2( i, z(b1), 1)
      ztz = 1.0e0/ztz
      call dscal(i, ztz, z(b1), 1)
c     
c + DLAMCHU
c      	
      RETURN
101   format (E22.14)
102   format (E22.14,E22.14)
103   format (E22.14,E22.14,E22.14)
*      
*     End of DGETAVEC
*      
      END
