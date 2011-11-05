      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
c
c  njtj
c  ###  Cray conversions
c  ###    1)Switch double precision to real.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
Cray      REAL D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
Cray      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,PFIVE=0.5D0)
Cray      PARAMETER(ZERO=0.0,ONE=1.0,TWO=2.0,PFIVE=0.5)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY,
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES,
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT-
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE,
C
C        D AND E ARE UNALTERED,
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO,
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES,
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER,
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP = TWO**(-40)
C
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = ZERO
C     ********** LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES **********
      DO 40 I = 1, N
         X1 = U
         U = ZERO
         IF (I .NE. N) U = ABS(E(I+1))
         XU = MIN(D(I)-(X1+U),XU)
         X0 = MAX(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         IF (ABS(E(I)) .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1))))
     X      GO TO 40
   20    E2(I) = ZERO
   40 CONTINUE
C
      X1 = MAX(ABS(XU),ABS(X0)) * MACHEP * REAL(N)
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     ********** DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES **********
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS **********
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = ZERO
C
      DO 120 Q = P, N
         X1 = U
         U = ZERO
         V = ZERO
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = MIN(D(Q)-(X1+U),XU)
         X0 = MAX(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0) GO TO 140
  120 CONTINUE
C
  140 X1 = MAX(ABS(XU),ABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     ********** CHECK FOR ISOLATED ROOT WITHIN INTERVAL **********
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * REAL(Q-P+1)
      LB = MAX(T1,XU-X1)
      UB = MIN(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     ********** FIND ROOTS BY BISECTION **********
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     ********** LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) **********
      K = M2
  250    XU = LB
C     ********** FOR I=K STEP -1 UNTIL M1 DO -- **********
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     ********** NEXT BISECTION STEP **********
  300    X1 = (XU + X0) * PFIVE
         IF ((X0 - XU) .LE. (TWO * MACHEP *
     X      (ABS(XU) + ABS(X0)) + ABS(EPS1))) GO TO 420
C     ********** IN-LINE PROCEDURE FOR STURM SEQUENCE **********
  320    S = P - 1
         U = ONE
C
         DO 340 I = P, Q
            IF (U .NE. ZERO) GO TO 325
            V = ABS(E(I)) / MACHEP
            IF (E2(I) .EQ. ZERO) V = ZERO
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. ZERO) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     ********** REFINE INTERVALS **********
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     ********** K-TH EIGENVALUE FOUND **********
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     ********** ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS **********
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     ********** SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES **********
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
C     ********** LAST CARD OF TRIDIB **********
      END
c $Id$
