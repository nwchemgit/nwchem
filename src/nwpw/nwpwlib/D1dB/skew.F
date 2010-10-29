*
* $Id$
*


*     ****************************************
*     *                                      *
*     *            Rotate_Skew               *
*     *                                      *
*     ****************************************

*  This routine returns an NxN Rotation matrix
*
*                 R = Exp(t*K) 
*
*  given  K = U*Sigma*U^H, where U = (V+i*W)
*
* where K is a skew matrix that is factored into 
* V, W, and Sigma using Factor_Skew.
*
*     Entry - N: matrix dimension
*            V: real factor matrix
*            W: imaginary factor matrix
*            Sigma: eigenvalues
*            t: rotational parameter
*            tmpv,tmpw: NxN tmp space
*
*     Exit - R: rotation matrix
*
      subroutine Rotate_Skew(N,V,W,Sigma,t,R,tmpv,tmpw)
      implicit none
      integer N
      real*8 V(N,N),W(N,N)
      real*8 Sigma(N)
      real*8 t
      real*8 R(N,N)
      real*8 tmpv(N,N),tmpw(N,N)

*     **** local variables ****
      integer I,J,K
      real*8 sum,csa,sna

      do I=1,N
       csa = cos(Sigma(I)*t)
       sna = sin(Sigma(I)*t) 
       do J=1,N
        tmpv(I,J) = csa*V(J,I)+sna*W(J,I)
        tmpw(I,J) = csa*W(J,I)-sna*V(J,I)
      end do
      end do

      do J=1,N
      do I=1,N
        sum = 0.0d0
        do K=1,N
          sum = sum + V(I,K)*tmpv(K,J) + W(I,K)*tmpw(K,J)
        end do
        R(I,J) = sum
      end do
      end do

      return
      end



*     ****************************************
*     *                                      *
*     *             Factor_Skew              *
*     *                                      *
*     ****************************************

*  This routine factors a NxN skew matrix K such that
*
*  K = U*Sigma*U^H, where U = (V+i*W)
*
*     Entry - N: matrix dimension
*             K: NxN skew matrix - destroyed on output
*     Exit - V: real factor matrix
*            W: imaginary factor matrix
*            Sigma: eigenvalues
*

      subroutine Factor_Skew(N,K,V,W,Sigma)
      implicit none
      integer N
      real*8 K(N,N)
      real*8 V(N,N),W(N,N)
      real*8 Sigma(N)

*     **** local variables ****                                     
      REAL*8   CON
      INTEGER   I,J,JP1,IERR
      LOGICAL MATZ,SKEW
      SKEW = .TRUE.
      MATZ = .TRUE.

*     **** COMPUTE EIGENVALUES AND EIGENVECTORS ****
      CALL TRIZD (N,N,K,Sigma)
      CALL IMZD (N,Sigma,MATZ,SKEW,N,V,IERR)
      IF (IERR .NE. 0) WRITE (6,*) "Factor Skew: IERR=",IERR
      CALL TBAKZD (N,N,K,N,N,V)
      call dcopy(N*N,V,1,K,1)


      J = 0
   15    J = J + 1
         IF (Sigma(J) .EQ. 0.0d0) GO TO 20
         JP1 = J+1
         DO I=1,N
           V(I,J) = K(I,J)  ! set the eigenvector
           W(I,J) = K(I,JP1)
         end do

         do I=1,N
            CON = -K(I,JP1)
            V(I,J+1) = K(I,J)   ! set the eigenvector
            W(I,J+1) = CON
         end do
         J = J + 1
         GO TO 30
   20 do I=1,N
        V(I,J) = K(I,J)
        W(I,J) = 0.0d0
      end do      
   30 IF (J .LT. N) GO TO 15

      call dscal(N*N,1.0d0/dsqrt(2.0d0),V,1)
      call dscal(N*N,1.0d0/dsqrt(2.0d0),W,1)

      return
      end 



      SUBROUTINE TRIZD ( NA, N, A, E )
C
C ****
C
C  FUNCTION  -  REDUCES A REAL*8 SKEW-SYMMETRIC MATRIX TO A SKEW-SYMMETRIC
C                 TRIDIAGONAL MATRIX USING ORTHOGONAL SIMILARITY
C                 TRANSFORMATIONS
C
C  PARAMETERS
C
C     NA       - INPUT INTEGER SPECIFYING THE ROW DIMENSION OF A AS
C                  DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT
C
C     N        - INPUT INTEGER SPECIFYING THE ORDER OF A
C
C     A(NA,N)  - ON INPUT, A CONTAINS THE REAL*8 SKEW-SYMMETRIC MATRIX.
C                  ONLY THE STRICT LOWER TRIANGLE OF THE MATRIX NEED
C                  BE SUPPLIED.
C                ON OUTPUT, A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C                  TRANSFORMATIONS USED IN THE REDUCTION IN ITS FULL
C                  LOWER TRIANGLE. THE STRICT UPPER TRIANGLE OF A IS
C                  UNALTERED.
C
C     E(N)     - OUTPUT ARRAY CONTAINING THE LOWER SUBDIAGONAL ELEMENTS
C                  OF THE TRIDIAGONAL MATRIX IN ITS LAST N-1 POSITIONS.
C                  E(1) IS SET TO ZERO.
C
C  REQUIRED FUNCTIONS - ABS,SIGN,SQRT
C
C ****
C
      implicit none
      integer NA,N
      REAL*8   A(NA,N),E(N)
      REAL*8   F,G,H,SCALE
c      REAL*8   ABS,SIGN,SQRT
      INTEGER   I,J,K,L,II,JM1,JP1
C
      IF (N .EQ. 1) GO TO 230
C
C *** MAIN DO LOOP  I=N STEP -1 UNTIL 2
C
      DO 220 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0d0
         SCALE = 0.0d0
C
C *** NORMALIZE ROW
C
         DO 100 K = 1, L
            SCALE = SCALE + DABS(A(I,K))
  100    CONTINUE
C
         IF (SCALE .NE. 0.0d0) GO TO 120
         E(I) = 0.0d0
         GO TO 215
C
C *** COMPUTE ELEMENTS OF U VECTOR
C
  120    DO 130 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  130    CONTINUE
C
         F = A(I,L)
         G = -SIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 200
C
C *** COMPUTE ELEMENTS OF A*U/H
C
         DO 180 J = 1, L
            G = 0.0d0
            IF (J .EQ. 1) GO TO 150
            JM1 = J - 1
C
            DO 140 K = 1, JM1
               G = G + A(J,K) * A(I,K)
  140       CONTINUE
C
  150       IF (J .EQ. L) GO TO 170
            JP1 = J + 1
C
            DO 160 K = JP1, L
               G = G - A(K,J) * A(I,K)
  160       CONTINUE
C
  170       E(J) = G / H
  180    CONTINUE
C
C *** COMPUTE REDUCED A
C
         DO 190 J = 2, L
            F = A(I,J)
            G = E(J)
            JM1 = J - 1
C
            DO 190 K = 1, JM1
               A(J,K) = A(J,K) + F * E(K) - G * A(I,K)
  190    CONTINUE
C
  200    DO 210 K = 1, L
            A(I,K) = SCALE * A(I,K)
  210    CONTINUE
C
  215    A(I,I) = SCALE * DSQRT(H)
  220 CONTINUE
C
  230 E(1) = 0.0d0
      RETURN
      END
      SUBROUTINE IMZD ( N, E, MATZ, SKEW, NZ, Z, IERR )
C
C ****
C
C  FUNCTION  -  COMPUTE THE EIGENVALUES AND OPTIONALLY THE EIGENVECTORS
C                 OF A SYMMETRIC TRIDIAGONAL MATRIX WITH ZERO DIAGONALS
C                 OR A SKEW-SYMMETRIC TRIDIAGONAL MATRIX USING AN
C                 IMPLICIT QR-TYPE ITERATION
C
C  PARAMETERS
C
C     N        - INPUT INTEGER SPECIFYING THE ORDER OF THE TRIDIAGONAL
C                  MATRIX
C
C     E(N)     - ON INPUT, ARRAY CONTAINING THE LOWER SUBDIAGONAL
C                  ELEMENTS OF THE TRIDIAGONAL MATRIX IN ITS LAST N-1
C                  POSITIONS. E(1) IS ARBITRARY.
C                ON OUTPUT, ARRAY CONTAINS THE EIGENVALUES. THE NON-ZERO
C                  EIGENVALUES OCCUR IN PAIRS WITH OPPOSITE SIGNS AND
C                  ARE FOUND IN ADJACENT LOCATIONS IN E. THE EIGENVALUES
C                  OF SYMMETRIC MATRICES ARE REAL*8 AND THE EIGENVALUES
C                  OF SKEW-SYMMETRIC MATRICES ARE PURELY IMAGINARY
C                  COMPLEX NUMBERS. IF AN ERROR EXIT IS MADE, THE
C                  EIGENVALUES ARE CORRECT FOR INDICES IERR+1,IERR+2...N
C
C     MATZ     - INPUT LOGICAL VARIABLE SPECIFYING THE EIGENVECTOR
C                  OPTION
C                  = .TRUE.   EIGENVECTORS ARE TO BE COMPUTED
C                  = .FALSE.  EIGENVECTORS ARE NOT TO BE COMPUTED
C
C     SKEW     - INPUT LOGICAL VARIABLE SPECIFYING TYPE OF INPUT MATRIX
C                  = .TRUE.   INPUT TRIDIAGONAL MATRIX IS SKEW-SYMMETRIC
C                  = .FALSE.  INPUT TRIDIAGONAL MATRIX IS SYMMETRIC WITH
C                               ZERO DIAGONALS
C                  SKEW IS NOT REFERENCED IF MATZ = .FALSE.
C
C     NZ       - INPUT INTEGER SPECIFYING THE ROW DIMENSION OF Z AS
C                  DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT
C
C     Z(NZ,N)  - OUTPUT ARRAY CONTAINING THE ORTHOGONAL EIGENVECTORS
C                  OF THE INPUT TRIDIAGONAL MATRIX. EIGENVECTORS CORRE-
C                  SPONDING TO ZERO EIGENVALUES ARE NORMALIZE TO UNIT
C                  2-NORM (LENGTH) AND THOSE CORRESPONDING TO NON-ZERO
C                  EIGENVALUES HAVE 2-NORM OF SQUARE ROOT 2. IF THE J-TH
C                  EIGENVALUE IS ZERO OR REAL*8 (I.E. E(J)), ITS EIGEN-
C                  VECTOR IS FOUND IN THE J-TH COLUMN OF Z. IF THE J-TH
C                  EIGENVALUE IS IMAGINARY (I.E. E(J)*I) WITH E(J+1) =
C                  -E(J), THE REAL*8 PART OF ITS EIGENVECTOR IS FOUND IN
C                  THE J-TH COLUMN OF Z AND ITS IMAGINARY PART FOUND IN
C                  THE (J+1)-TH COLUMN. IF AN ERROR EXIT IS MADE, Z
C                  CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C                  EIGENVALUES.
C                  Z IS NOT REFERENCED IF MATZ = .FALSE.
C
C     IERR     - OUTPUT ERROR CODE
C                  = 0   NORMAL RETURN (ALL EIGENVALUES/VECTORS FOUND)
C                  = J   IF THE J-TH EIGENVALUE HAS NOT BEEN DETERMINED
C                          AFTER 30 ITERATIONS
C
C  REQUIRED FUNCTIONS - ABS,SIGN,SQRT,MOD
C
C ****
C
      implicit none
      INTEGER  N,NZ
      REAL*8   E(N),Z(NZ,N)
      REAL*8   F,G,Q,C,S,R,P,TEST,TMAG
c      REAL*8   ABS,SIGN,SQRT
      INTEGER   I,J,K,L,M,L0,L0M1,M0,LS,IM1,JP1,KM1,KP1,LM1,MM1,
     1          IP3,IERR,IEO,ITS,IP1
c      INTEGER   MOD
      LOGICAL MATZ,SKEW
C
      IF (.NOT. MATZ) GO TO 115
C
C *** PLACE IDENTITY MATRIX IN Z
C
      DO 110 I = 1, N
         DO 100 J = 1, N
            Z(I,J) = 0.0d0
  100    CONTINUE
         Z(I,I) = 1.0d0
  110 CONTINUE
C
  115 IERR = 0
      M = N
      MM1 = M - 1
      E(1) = 0.0d0
      ITS = 0
C
  120 IF (M .LT. 2) GO TO 370
      M0 = M
C
C *** SEARCH FOR NEXT SUBMATRIX TO SOLVE  (MATRIX SPLITTING)
C
      F = 0.0d0
      DO 130 I = 1, MM1
         J = M - I
         JP1 = J + 1
         G = DABS(E(JP1))
         TMAG = DABS(E(J)) + F
         TEST = TMAG + G
         IF (TEST .EQ. TMAG) GO TO 140
         F = G
  130 CONTINUE
      JP1 = 1
C
  140 L0 = JP1 + 1
      L0M1 = JP1
      IF (L0M1 .EQ. M) GO TO 290
      IF (.NOT. MATZ) GO TO 160
      IF (.NOT. SKEW) GO TO 160
C
C *** PLACE CORRECT SIGN ON IDENTITY DIAGONALS
C
      DO 150 I = L0M1, M, 4
         Z(I,I) = -Z(I,I)
         IP3 = I + 3
         IF (IP3 .GT. M) GO TO 160
         Z(IP3,IP3) = -Z(IP3,IP3)
  150 CONTINUE
C
  160 IF (L0 .EQ. M) GO TO 300
      IEO = M - L0
      IEO = MOD(IEO,2)
      L = L0
      IF (IEO .EQ. 0) GO TO 230
C
C *** FIND ZERO EIGENVALUE OF ODD ORDERED SUBMATRICES
C
      C = 0.0d0
      S = -1.0d0
      DO 190 I = L0, MM1, 2
         K = MM1 + L0 - I
         KP1 = K + 1
         Q = -S * E(KP1)
         E(KP1) = C * E(KP1)
         IF (DABS(E(K)) .GT. ABS(Q)) GO TO 170
         C = E(K) / Q
         R = DSQRT(C*C + 1.0d0)
         E(K) = Q * R
         S = 1.0d0 / R
         C = C * S
         GO TO 180
  170    S = Q / E(K)
         R = DSQRT(1.0d0 + S*S)
         E(K) = E(K) * R
         C = 1.0d0 / R
         S = S * C
  180    IF (.NOT. MATZ) GO TO 190
C
C *** ACCUMULATE TRANSFORMATIONS FOR EIGENVECTORS
C
         KM1 = K - 1
         Z(KM1,M) = -S * Z(KM1,KM1)
         Z(KM1,KM1) = C * Z(KM1,KM1)
         DO 185 J = KP1, M, 2
            Z(J,KM1) = S * Z(J,M)
            Z(J,M) = C * Z(J,M)
  185    CONTINUE
C
  190 CONTINUE
      M = MM1
      MM1 = M - 1
      IF (L0 .EQ. M) GO TO 300
C
C *** CHECK FOR CONVERGENCE OR SMALL SUBDIAGONAL ELEMENT
C
  200 DO 210 I = L0, MM1, 2
         K = MM1 + L0 - I
         L = K + 1
         TMAG = DABS(E(L)) + DABS(E(K-1))
         TEST = TMAG + E(K)
         IF (TEST .EQ. TMAG) GO TO 220
  210 CONTINUE
      L = L0
  220 IF (L .EQ. M) GO TO 300
C
C *** FORM SHIFT
C
  230 ITS = ITS + 1
      IF (ITS .GT. 30) GO TO 360
      F = E(M-3)
      G = E(M-2)
      C = E(MM1)
      S = E(M)
      P = ((C-F) * (C+F) + (S-G) * (S+G)) / (2. * G * C)
      R = DSQRT(P*P + 1.0d0)
      Q = (G / (P + SIGN(R,P))) - C
      F = E(L)
      LM1 = L - 1
      E(LM1) = ((F-S) * (F+S) + C * Q) / F
C
C *** PERFORM ONE IMPLICIT QR ITERATION ON CHOLESKY FACTOR
C
      LS = L0M1
      C = 1.0d0
      S = 1.0d0
      DO 280 I = L, MM1
         IP1 = I + 1
         IM1 = I - 1
         Q = S * E(IP1)
         E(IP1) = C * E(IP1)
         IF (DABS(E(IM1)) .GT. DABS(Q)) GO TO 240
         C = E(IM1) / Q
         R = DSQRT(C*C + 1.0d0)
         E(IM1) = Q * R
         S = 1.0d0 / R
         C = C * S
         GO TO 250
  240    S = Q / E(IM1)
         R = DSQRT(1.0d0 + S*S)
         E(IM1) = E(IM1) * R
         C = 1.0d0 / R
         S = S * C
  250    F = E(IP1)
         E(IP1) = -S * E(I) + C * F
         E(I) = C * E(I) + S * F
         IF (.NOT. MATZ) GO TO 280
C
C *** ACCUMULATE TRANSFORMATIONS FOR EIGENVECTORS
C
         DO 260 J = LS, M0, 2
            F = Z(J,IP1)
            Z(J,IP1) = -S * Z(J,IM1) + C * F
            Z(J,IM1) = C * Z(J,IM1) + S * F
  260    CONTINUE
         IF (LS .EQ. L0M1) GO TO 270
         LS = L0M1
         GO TO 280
  270    LS = L0
  280 CONTINUE
      E(LM1) = 0.0d0
      GO TO 200
C
C *** ITERATION CONVERGED TO ONE ZERO EIGENVALUE
C
  290 E(M) = 0.0d0
      M = MM1
      GO TO 310
C
C *** ITERATION CONVERGED TO EIGENVALUE PAIR
C
  300 E(MM1) = E(M)
      E(M) = -E(M)
      M = M - 2
C
  310 ITS = 0
      MM1 = M - 1
      IF (M .GT. L0) GO TO 200
      IF (M .EQ. L0) GO TO 300
      IF (.NOT. MATZ) GO TO 120
      IF (SKEW) GO TO 120
C
C *** COMPUTE EIGENVECTORS FROM ORTHONORMAL COLUMNS OF Z IF NOT SKEW
C
  320 K = M0
  330 IF (E(K) .EQ. 0.0d0) GO TO 350
      KM1 = K - 1
      DO 340 J = L0M1, M0, 2
         Z(J,K) = Z(J,KM1)
         F = Z(J+1,K)
         Z(J+1,KM1) = F
         Z(J+1,K) = -F
  340 CONTINUE
      K = KM1
  350 K = K-1
      IF (K .GT. L0M1) GO TO 330
      IF (IERR .NE. 0) GO TO 370
      GO TO 120
C
C *** ERROR EXIT
C
  360 IERR = M
      IF (.NOT. MATZ) GO TO 370
      IF (.NOT. SKEW) GO TO 320
C
  370 RETURN
      END
      SUBROUTINE TBAKZD ( NA, N, A, M, NZ, Z )
C
C ****
C
C  FUNCTION  -  FORMS THE EIGENVECTORS OF A REAL*8 SKEW-SYMMETRIC MATRIX
C                 BY BACK TRANSFORMING THOSE OF THE CORRESPONDING SKEW-
C                 SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRIZD
C
C  PARAMETERS
C
C     NA       - INPUT INTEGER SPECIFYING THE ROW DIMENSION OF A
C                  AS DECLARED IN CALLING PROGRAM DIMENSION STATEMENT
C
C     N        - INPUT INTEGER SPECIFYING THE ORDER OF A
C
C     A(NA,N)  - INPUT ARRAY CONTAINING INFORMATION ABOUT THE ORTHOGONAL
C                  TRANSFORMATIONS USED IN THE REDUCTION BY  TRIZD  IN
C                  ITS FULL LOWER TRIANGLE
C
C     M        - INPUT INTEGER SPECIFYING THE NUMBER OF EIGENVECTORS TO
C                  BE BACK TRANSFORMED
C
C     NZ       - INPUT INTEGER SPECIFYING THE ROW DIMENSION OF Z AS
C                  DECLARED IN CALLING PROGRAM DIMENSION STATEMENT
C
C     Z(NZ,M)  - ON INPUT, Z CONTAINS THE REAL*8 AND IMAGINARY (IF
C                  COMPLEX) PARTS OF THE EIGENVECTORS TO BE BACK TRANS-
C                  FORMED IN ITS FIRST M COLUMNS
C                ON OUTPUT, Z CONTAINS THE REAL*8 AND IMAGINARY (IF
C                  COMPLEX) PARTS OF THE TRANSFORMED EIGENVECTORS IN
C                  ITS FIRST M COLUMNS
C
C ****
C
      implicit none
      INTEGER  NA,N,M,NZ
      REAL*8   A(NA,N),Z(NZ,M)
      REAL*8   H,S
      INTEGER   I,J,K,L
C
      IF (M .EQ. 0) GO TO 140
      IF (N .EQ. 1) GO TO 140
C
      DO 130 I = 2, N
         L = I - 1
         H = A(I,I)
         IF (H .EQ. 0.0d0) GO TO 130
C
         DO 120 J = 1, M
            S = 0.0d0
C
            DO 100 K = 1, L
               S = S + A(I,K) * Z(K,J)
  100       CONTINUE
C
            S = (S / H) / H
C
            DO 110 K = 1, L
               Z(K,J) = Z(K,J) - S * A(I,K)
  110       CONTINUE
C
  120    CONTINUE
C
  130 CONTINUE
C
  140 RETURN
      END
