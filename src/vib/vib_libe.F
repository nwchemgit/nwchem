      SUBROUTINE vib_LIBE(EQVAL,NOINT,NLIBE,I,J,K,C,B,NDIM)
C*MODULE VIB  *DECK LIBE
* $Id$
C
      IMPLICIT NONE ! DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
C
      integer ndim, noint, nlibe, i, j, k
      double precision B(NDIM,*),C(3,*), eqval
      DOUBLE PRECISION ZERO, ONE, TOL1, TOL2
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (TOL1=1.0D-04, TOL2=5.0D-04)
      INTEGER MAXBE, NNAM
      PARAMETER (MAXBE=5, NNAM=1)
C
      double precision A(3),RJI(3),RJK(3),UN(3),APTS(3*MAXBE),RJA(3),
     *          UNIT(3),UP(3),EJI(3),EJK(3)
c
      DOUBLE PRECISION PI, FACT, DJKSQ, DJISQ, DAJSQ, DOT, RJASQ
      DOUBLE PRECISION DJI, DJK, DAJ, DOTJ, DOTP, TEST, B1, B2
      INTEGER IPT, M, LOC, NOCOL1, NOCOL2, NOCOL3 
C
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR ONE OR
C           BOTH OF 2 PERPENDICULAR LINEAR BENDING COORDINATES.
C           SEE FERIGLE AND MEISTER, J CHEM PHYS 19,982(1951) -----
C
C     ----- I AND K ARE THE END ATOMS AND J IS THE CENTRAL ATOM. A
C           GIVES THE CARTESIAN COORDINATES OF A POINT IN SPACE SUCH
C           THAT THE VECTOR FROM ATOM J TO POINT A IS PERPENDICULAR TO
C           THE LINE I-J-K AND SERVES TO ORIENT THE COORDS IN SPACE.
C           THE FIRST OF THE TWO INTERNAL COORDINATES IS IN THE I-A-K
C           PLANE AND THE SECOND IN A PLANE PERPENDICULAR TO THE FIRST,
C           THROUGH POINTS I,J, AND K -----
C
      PI = ACOS(-ONE)
C
C     ----- READ THE PLANE CONTAINING THE POINTS A -----
C
       IPT = 1 + (NLIBE - 1)*3
       READ(5,*,ERR=100,END=110)(APTS(I),I=IPT,IPT+2)
       GOTO 115
  100  WRITE(6,9040)
       CALL errquit('vib_libe error',911, UNKNOWN_ERR)
  110  WRITE(6,9050)
       CALL errquit('vib_libe error',911, UNKNOWN_ERR)
C
C     ----- TAKE THE POINT A LYING IN THE PLANE, AND GENERATE
C           THE POINT A SUCH THAT A-J IS PERP. TO I-J-K -----
C
  115 DJKSQ=ZERO
      DOT=ZERO
      LOC=3*NLIBE-3
      DO 120 M=1,3
         A(M) = APTS(LOC+M)
         RJK(M) = C(M,K) - C(M,J)
         DJKSQ = DJKSQ + RJK(M)*RJK(M)
         RJA(M) = A(M) - C(M,J)
         DOT = DOT + RJK(M)*RJA(M)
  120 CONTINUE
      FACT=DOT/DJKSQ
      RJASQ = ZERO
      DO 140 M=1,3
         RJA(M) = RJA(M) - FACT*RJK(M)
         RJASQ = RJASQ + RJA(M)*RJA(M)
         A(M) = C(M,J) + RJA(M)
  140 CONTINUE
      IF(ABS(RJASQ).LT.TOL1) THEN
         WRITE(6,9010) (APTS(LOC+M),M=1,3)
         CALL errquit('vib_libe error',911, UNKNOWN_ERR)
      END IF
C
C     ----- A IS NOW CORRECTLY LOCATED -----
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DAJSQ = ZERO
      DO 210 M = 1,3
         RJI(M) = C(M,I)-C(M,J)
         DJISQ = DJISQ+RJI(M)*RJI(M)
         RJK(M) = C(M,K)-C(M,J)
         DJKSQ = DJKSQ+RJK(M)*RJK(M)
         UN(M) = A(M)-C(M,J)
         DAJSQ = DAJSQ+UN(M)*UN(M)
  210 CONTINUE
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DAJ = SQRT(DAJSQ)
      DOTJ = ZERO
      DOTP = ZERO
      DO 230 M = 1,3
         EJI(M) = RJI(M)/DJI
         EJK(M) = RJK(M)/DJK
         UNIT(M) = UN(M)/DAJ
         DOTJ = DOTJ+EJI(M)*EJK(M)
         DOTP = DOTP+EJI(M)*UNIT(M)
  230 CONTINUE
      TEST = (ABS(DOTJ)-ONE)
      IF ((ABS(TEST).GT.TOL1) .OR.
     *    (ABS(DOTP).GT.TOL2)) THEN
         WRITE (6,9020)  TEST,DOTP
         CALL errquit('vib_libe error',911, UNKNOWN_ERR)
      END IF
      IF (ABS(DOTJ)-ONE .GT. ZERO) DOTJ=SIGN(ONE,DOTJ)
      EQVAL = ACOS(DOTJ)
      IF(ABS(EQVAL)   .LT.TOL1) EQVAL=ZERO
      IF(ABS(EQVAL-PI).LT.TOL1) EQVAL=PI
C
C     ----- FIRST COMPONENT OF THE BEND -----
C
      NOCOL1 = 3*(I-1)
      NOCOL2 = 3*(K-1)
      NOCOL3 = 3*(J-1)
      DO 300 M = 1,3
         B(NOINT,NOCOL1+M) = -UNIT(M)/DJI
         B(NOINT,NOCOL2+M) = -UNIT(M)/DJK
         B(NOINT,NOCOL3+M) = (ONE/DJI+ONE/DJK)*UNIT(M)
  300 CONTINUE
C
C     ----- SECOND COMPONENT OF THE BEND -----
C
      UP(1) = EJK(2)*UNIT(3)-EJK(3)*UNIT(2)
      UP(2) = EJK(3)*UNIT(1)-EJK(1)*UNIT(3)
      UP(3) = EJK(1)*UNIT(2)-EJK(2)*UNIT(1)
      DO 400 M = 1,3
         B1 = -UP(M)/DJI
         B2 = -UP(M)/DJK
         B(NOINT+1,NOCOL1+M) = B1
         B(NOINT+1,NOCOL2+M) = B2
         B(NOINT+1,NOCOL3+M) = -(B1+B2)
  400 CONTINUE
      RETURN
C
 9000 FORMAT(1X,'NO $LIBE INPUT DATA FOUND...THIS IS REQUIRED')
 9010 FORMAT(1X,'THE POINT A IS COLINEAR WITH I,J,K'/
     *    1X,'A=',3F15.8)
 9020 FORMAT(1X,'*** ROUNDOFF ERROR IN LIBE - TEST,DOTP=',2E20.10)
 9030 FORMAT(1X,'ONLY 5 LINEAR BENDS MAY BE USED --- ERROR')
 9040 FORMAT(1X,'Error reading input for linear bend info.',/,
     *       1X,'Call Bomb and Exit.')
 9050 FORMAT(1X,'End of file encountered during read of linear ',
     *          'bend infor.   Call Bomb and Exit.')
      END
