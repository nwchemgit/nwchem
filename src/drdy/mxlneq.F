C $Id$
      Subroutine drdy_mxlneq (a,nn,ida,dett,jrank,eps,in,mm)
C
C PROGRAMMED BY R. HOTCHKISS, U. COMP. CTR., U. OF MINN., REVISED OCT.73
C      modified for the VAX by Bruce Garrett, Nov. 1980
C*
C*    Dimensions changed by Tom Joseph Jan. 1988.
C*
C*    Solves for the matrix X in C*X = B. C must be a sq. matrix
C*    A(NN,NN+abs(MM)) contains matrices C and B, A(i,j) = C(i,j)
C*    A(i,j+NN) = B(i,j)
C*    IDA - first dim. of A, NN - order of matrix C (ie. # of eqns.)
C*    IN(NN) - work array
C*    M = 0 is used to calculate inverse of C
C*
C     CALLED BY:
C                FIVPT,NEWT,PROJCT,TREPT
C
      Implicit Double Precision (A-H,O-Z)
#include "errquit.fh"
C
      Double Precision MACHEP
      DIMENSION A(IDA,NN),IN(IDA)
C
C     DATA ZERO,ONE/0.0D0,1.0D0/,MACHEP/1.D35/
C     REVISED BY N.ABUSALBI TO ALLOW A DET OF THE ORDER OF 10**37
C
      DATA ZERO,ONE / 0.0D0,1.0D0 /,MACHEP / 1.D35 /
      DATA INTMX / 400000 /
      Integer fu6
      Data fu6/6/
C
C INITIATE SOME LOCAL AND OUTPUT VARIABLES
C
      JRANK = NN
      N = NN
      ID = IDA
      APS = EPS
      OPS = EPS
C
C error checking
C
      IF (N.GT.0.AND.N.LE.ID) GO TO 10
      WRITE(FU6,1700)
      WRITE(FU6,1000) N,ID
*     CALL EXIT
      call errquit('drdy_mxlneq: fatal error',911,
     &       UNKNOWN_ERR)
   10 IF (ID.GE.N.AND.ID.LT.INTMX) GO TO 20
      WRITE(FU6,1700)
      WRITE(FU6,1100) ID,N,INTMX
*     CALL EXIT
      call errquit('drdy_mxlneq: fatal error',911, UNKNOWN_ERR)
   20 IF ((ONE+OPS).NE.ONE.AND.OPS.GE.ZERO) GO TO 30
      WRITE(FU6,1700)
      WRITE(FU6,1200) OPS
*     CALL EXIT
      call errquit('drdy_mxlneq: fatal error',911, UNKNOWN_ERR)
   30 MOD = MM
      NM = IABS(MM)
      IF (NM.LT.INTMX) GO TO 40
      WRITE(FU6,1700)
      WRITE(FU6,1300) MM,INTMX
*     CALL EXIT
      call errquit('drdy_mxlneq: fatal error',911, UNKNOWN_ERR)
   40 NM = NN+NM
      NMIDA = NM*ID*2
      IF (NMIDA.LT.INTMX) GO TO 50
      WRITE(FU6,1700)
      WRITE(FU6,1400) NMIDA,INTMX
*     CALL EXIT
      call errquit('drdy_mxlneq: fatal error',911, UNKNOWN_ERR)
   50 CONTINUE
C
C  end of error checking
C CONTINUE TO INITIALIZE
C
      K1 = 1
      NFLAG = 0
      DET = ONE
C
C MOD IS -1 DET ONLY (this version does not have this option)
C        = or > 0 FOR INV, DET AND 0 OR MORE SETS OF LIN EQNS
C        <-1 FOR LIN EQNS ONLY AND DET
C MAIN GAUSS JORDAN LOOP BEGINS
C
      DO 150 K = 1, N
C
C SEARCH FOR LARGEST PIVOT CANDIDATE IN REMAINING LOWER RIGHT SQUARE
C  MATRIX
C
         PIV = ZERO
         L = K
         DO 60 I = K, N
            P = ABS(A(I,K))
            IF (PIV.GE.P) GO TO 60
            PIV = P
            L = I
   60    CONTINUE
C
C PIVOT WITH ABS VALUE PIV AND SUBSCRIPTS L AND M HAS BEEN FOUND
C
         PIVOT = A(L,K)
C
C CONTINUE IF PIV LARGER THAN USER EPS
C
         IF (PIV.GT.OPS) GO TO 80
         IF (EPS.EQ.OPS) JRANK = K-1
C
C EPS TEST FAILED, CHECK FOR ZERO PIVOT
C
         IF (PIV.GT.ZERO) GO TO 70
C
C PIVOT IS ZERO, TERMINATE PROGRAM UNLESS MOD=-1,IE, DET ONLY CASE
C      this version does not have the det only mode
C ZERO PIVOT MEANS ZERO DET AND EXIT IF DET ONLY MODE
C
         WRITE(FU6,1700)
         WRITE(FU6,1500) K
         jrank = 0
         return
C
C ISSUE NON-FATAL MESSAGE, PIV .LE. EPS
C
   70    WRITE(FU6,1700)
         WRITE(FU6,1600) K,EPS
C
C SET OPS TO 0 SO SOLUTION MAY CONTINUE AFTER ERROR MESSAGE
C
         OPS = ZERO
         PIV = PIVOT
C
C CALCULATE DETERMINANT AND CHECK FOR OVERFLOW
C
   80    DET = PIVOT*DET
         IF (ABS(DET).LT.MACHEP) GO TO 90
         WRITE(FU6,1700)
         WRITE(FU6,1800) K,DET
         jrank = 0
         return
C
C RESET LEADING ROW DO INDEX FOR DET ONLY AND LIN EQN ONLY CASE
C
   90    IF (MOD.LT.0) K1 = K
C
C SAVE PIVOT INDEX
C
         IN(K) = L
C
C CHECK FOR ROW INTERCHANGE
C
         IF (L.EQ.K) GO TO 110
         DET = -DET
C
C INTERCHANGE ROW CONTAINING PIVOT AND CURRENT ROW
C ONLY PARTIAL ROWS NEED BE EXCHANGED FOR DET ONLY OR LIN EQN ONLY
C  SINCE LOWER LEFT PARTIALLY FORMED TRIANGLE IS NOT NEEDED
C
         DO 100 J = K1, NM
            Z = A(L,J)
            A(L,J) = A(K,J)
            A(K,J) = Z
  100    CONTINUE
C
C PIVOT ELEMENT IS NOW ON DIAGONAL
C SAVE DIVISION TIME BY USING RECIPROCAL OF PIVOT
C
  110    PIVOT = ONE/PIVOT
C
C PRE-DIVIDE NECESSARY PORTION OF PIVOT ROW
C
         DO 120 J = K1, NM
            A(K,J) = A(K,J)*PIVOT
  120    CONTINUE
C
C SET PIVOT ELEMENT TO ZERO SO MAIN REDUCTION STEP DOESNT OPERATE ON
C  PIVOT ROW
C
         A(K,K) = ZERO
C
C SWEEP THROUGH ALL OR PART OF MATRIX USING KTH ROW, PIVOT ROW, TO
C  REDUCE THE MATRIX
C
         DO 140 I = K1, N
            Z = A(I,K)
            IF (Z.EQ.ZERO) GO TO 140
C
C THIS CHECK NOT ONLY PREVENTS OPERATING ON PIVOT ROW BUT CATCHES
C  OTHER ZEROES IN PIVOT COLUMNS. THESE OTHER ZEROES WOULD LEAVE JTH
C  ROW UNCHANGED IN FOLLOWING LOOP SO CONSIDERABLE TIME MAY BE SAVED BY
C  SKIPPING OPERATION
C
            DO 130 J = K1, NM
               A(I,J) = A(I,J)-Z*A(K,J)
  130       CONTINUE
C
C THE INVERSE IS CREATED IN PLACE BY SUBSTITUTING AN IDENTITY MATRIX
C  COL BY COL, SINCE WE ARE SUBT. THE PIVOT ROW FROM OFF DIAGONAL 0
C  ELEMENTS AT THIS POINT, WE NOW PLACE -A(I,K)/A(K,K) AT THIS POINT IN
C  THE PIVOT COL
C
            A(I,K) = -Z*PIVOT
  140    CONTINUE
C
C SIMILARLY DIVIDING PIVOT ROW BY THE PIVOT IS EQUIVALENT TO PLACING
C  ONE/A(K,K) AT THE PIVOT POINT FOR THE INVERSE
C
         A(K,K) = PIVOT
  150 CONTINUE
      IF (N.EQ.1) GO TO 210
      IF (MOD.GE.0) GO TO 180
C
C BACK SUBSTITUTION FOR LIN EQN ONLY CASE
C
      K1 = K1+1
      DO 170 K = K1, NM
         I = N
         DO 170 L = 2, N
            I1 = I
            I = I-1
            Z = ZERO
            DO 160 J = I1, N
               Z = Z+A(I,J)*A(J,K)
  160       CONTINUE
            A(I,K) = A(I,K)-Z
  170 CONTINUE
      GO TO 210
C
C FINAL REORDERING OF MATRIX
C
  180 K = N
C
C SKIP LAST STEP SINCE NO INTERCHANGE COULD OCCUR THERE
C
      DO 200 J = 2, N
C
C PERFORM INTERCHANGES IN EXACT REVERSE ORDER OF PREVIOUS EXECUTION
C
         K = K-1
C
C ROW INTERCHANGE DURING INVERSION IMPLIES COL INTERCHANGE HERE
C
         M = IN(K)
         IF (M.EQ.K) GO TO 200
C
C COL INTERCHANGE
C
         DO 190 I = 1, N
            Z = A(I,K)
            A(I,K) = A(I,M)
            A(I,M) = Z
  190    CONTINUE
  200 CONTINUE
  210 DETT = DET
      RETURN
C
 1000 FORMAT(' arg 2, n =',I10,/,7X,' id=',I10,
     &       ' n must be .ge. 1 and .le. id')
 1100 FORMAT(' arg 3, id =',I10,/,9x,' n=',I10,
     &       ' id must be .ge. N and .le.',I10)
 1200 FORMAT(' arg 6, eps =',1PE13.5,
     &       ' eps must be .ge. zero and finite')
 1300 FORMAT(' arg 8, m =',I10,
     &       ' abs(m) must be .le.',I10)
 1400 FORMAT(' size =',I10,/,
     &       ' size = id*(n+abs(m))*2 must be .le.',I10)
 1500 FORMAT(' k =',I10,/,
     &       ' at step k a gauss-jordan pivot value was zero')
 1600 FORMAT(' k =',I10,/,' eps =',1PE13.5,/,
     &   ' at step k a gauss-jordan pivot value was .le. eps ')
 1700 FORMAT(' *** mxlneq ***' )
 1800 FORMAT(' k =',I10,/,' det =' ,1PE13.5,/,
     &       ' at step k det is too large' )
*holerith_replaced: 1000 FORMAT(11H arg 2, n = ,I10/7X,4H id=,I10,
*holerith_replaced:     &       29H n must be .ge. 1 and .le. id)
*holerith_replaced: 1100 FORMAT(12H arg 3, id = ,I10 /9x,3h n=,I10,
*holerith_replaced:     &       27H id must be .ge. N and .le. I10)
*holerith_replaced: 1200 FORMAT(13H arg 6, eps = ,1PE13.5,
*holerith_replaced:     &       33H eps must be .ge. zero and finite)
*holerith_replaced: 1300 FORMAT(11H arg 8, m = ,I10,
*holerith_replaced:     &   20H abs(m) must be .le. I10 )
*holerith_replaced: 1400 FORMAT(7H size = ,I10/
*holerith_replaced:     &       36H size = id*(n+abs(m))*2 must be .le. I10)
*holerith_replaced: 1500 FORMAT(4H k =,I10/
*holerith_replaced:     &       46H at step k a gauss-jordan pivot value was zero)
*holerith_replaced: 1600 FORMAT(4H k = ,I10/6H eps = ,1PE13.5/
*holerith_replaced:     &   50H at step k a gauss-jordan pivot value was .le. eps )
*holerith_replaced: 1700 FORMAT(15H *** mxlneq *** )
*holerith_replaced: 1800 FORMAT(4H k = ,I10/6H det = ,1PE13.5/
*holerith_replaced:     &        27H at step k det is too large )
C
      END
