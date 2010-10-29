      SUBROUTINE vib_BEND(EQVAL,NOINT,I,J,K,C,B,NDIM)
* $Id$
C*MODULE VIB  *DECK BEND
C
      IMPLICIT NONE
#include "errquit.fh"
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS OF A
C          VALENCE ANGLE BENDING COORDINATE AS DEFINED BY WILSON.
C          SEE PP. 56-58 OF "MOLECULAR VIBRATIONS" BY E.B. WILSON,
C          J.C. DECIUS AND P.C. CROSS, MCGRAW-HILL, NY, 1955 ----
C
C     -----I AND K ARE THE NUMBERS OF THE END ATOMS.  J IS THE
C          NUMBER OF THE CENTRAL ATOM -----
C
C
      integer noint, i, j, k, ndim
      double precision eqval
      double precision C(3,*),B(NDIM,*)
c:old      DIMENSION C(3,1),B(NDIM,1)
c:old      DIMENSION RJI(3),RJK(3),EJI(3),EJK(3)
C
      DOUBLE PRECISION ZERO, ONE, TOL
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00, TOL=1.0D-07)
c
      double precision RJI(3),RJK(3),EJI(3),EJK(3)
      double precision pi
      double precision djisq, djksq
      double precision dji  , djk
      double precision dotj, sinj, b1, b2, dotval
      integer m, nocol1, nocol2, nocol3
c
C
C
      PI = ACOS(-ONE)
      DJISQ = ZERO
      DJKSQ = ZERO
C
      DO 100 M = 1,3
         RJI(M) = C(M,I)-C(M,J)
         RJK(M) = C(M,K)-C(M,J)
         DJISQ = DJISQ+RJI(M)*RJI(M)
         DJKSQ = DJKSQ+RJK(M)*RJK(M)
  100 CONTINUE
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DOTJ = ZERO
      DO 120 M = 1,3
         EJI(M) = RJI(M)/DJI
         EJK(M) = RJK(M)/DJK
         DOTJ = DOTJ+EJI(M)*EJK(M)
  120 CONTINUE
      IF (ONE.LE.ABS(DOTJ)) GO TO 200
      SINJ = SQRT(ONE-DOTJ*DOTJ)
      NOCOL1 = 3*(I-1)
      NOCOL2 = 3*(K-1)
      NOCOL3 = 3*(J-1)
      DO 160 M = 1,3
         B1 = ((DOTJ*EJI(M)-EJK(M)))/(DJI*SINJ)
         B2 = ((DOTJ*EJK(M)-EJI(M)))/(DJK*SINJ)
         B(NOINT,NOCOL1+M) = B1
         B(NOINT,NOCOL2+M) = B2
         B(NOINT,NOCOL3+M) = -(B1+B2)
  160 CONTINUE
      DOTVAL = ZERO
      DO 180 M = 1,3
  180 DOTVAL = DOTVAL+EJI(M)*EJK(M)
      EQVAL =   ACOS(DOTVAL)
      IF(ABS(EQVAL).LT.TOL) EQVAL=ZERO
      IF(ABS(EQVAL-PI).LT.TOL) EQVAL=PI
      RETURN
C
  200 CONTINUE
      WRITE (6,9008)
      CALL errquit(' vib_bend',911, UNKNOWN_ERR)
 9008 FORMAT(31H ROUNDOFF ERROR IN BEND - STOP )
      END
