      SUBROUTINE vib_TORS(EQVAL,NOINT,I,J,K,L,C,B,NDIM)
C*MODULE VIB  *DECK TORS
* $Id$
C
      IMPLICIT NONE ! DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR THE
C          TORSION AS DEFINED BY WILSON.
C          SEE PP. 60-61 OF "MOLECULAR VIBRATIONS" BY E.B. WILSON,
C          J.C. DECIUS AND P.C. CROSS, MCGRAW-HILL, NY, 1955.
C          NOTE THAT THE SIGN OF THE SECOND TERM ON THE RIGHT HAND
C          SIDE IN EQUATION (22) IS INCORRECT AS PRINTED ON P61 OF
C          WDC -----
C
C     -----I AND L ARE THE END ATOMS.  I IS NEARER THE OBSERVER.
C          J AND K ARE THE CENTRAL ATOM. J IS NEARER THE OBSERVER -----
C
      integer noint, i, j, k, l, ndim
      double precision C(3,*),B(NDIM,*), EQVAL
c
#include "util.fh"      
c
      double precision zero, one, tol, tolrd
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00, TOL=1.0D-07)
      PARAMETER (TOLRD=1.0001D+00)
c
      double precision RIJ(3), RJK(3), RKL(3)
      double precision EIJ(3), EJK(3), EKL(3)
      double precision CR1(3), CR2(3), PI
      double precision dijsq, djksq, dklsq, dij, djk, dkl
      double precision dotpj, dotpk, sinpj, sinpk
      double precision smi, smj, sml, sense, f1, f2, dot
      integer m, nocol1, nocol2, nocol3, nocol4
C
C
C
      PI = ACOS(-ONE)
      DIJSQ = ZERO
      DJKSQ = ZERO
      DKLSQ=ZERO
      DO 120 M = 1,3
         RIJ(M) = C(M,J)-C(M,I)
         DIJSQ = DIJSQ+RIJ(M)*RIJ(M)
         RJK(M) = C(M,K)-C(M,J)
         DJKSQ = DJKSQ+RJK(M)*RJK(M)
         RKL(M) = C(M,L)-C(M,K)
         DKLSQ = DKLSQ+RKL(M)*RKL(M)
  120 CONTINUE
      DIJ = SQRT(DIJSQ)
      DJK = SQRT(DJKSQ)
      DKL = SQRT(DKLSQ)
      DO 180 M = 1,3
         EIJ(M) = RIJ(M)/DIJ
         EJK(M) = RJK(M)/DJK
         EKL(M) = RKL(M)/DKL
  180 CONTINUE
C
      CR1(1) = EIJ(2)*EJK(3)-EIJ(3)*EJK(2)
      CR1(2) = EIJ(3)*EJK(1)-EIJ(1)*EJK(3)
      CR1(3) = EIJ(1)*EJK(2)-EIJ(2)*EJK(1)
C
      CR2(1) = EJK(2)*EKL(3)-EJK(3)*EKL(2)
      CR2(2) = EJK(3)*EKL(1)-EJK(1)*EKL(3)
      CR2(3) = EJK(1)*EKL(2)-EJK(2)*EKL(1)
C
      DOTPJ = ZERO
      DOTPK = ZERO
      DO 220 M = 1,3
         DOTPJ = DOTPJ-EIJ(M)*EJK(M)
         DOTPK = DOTPK-EJK(M)*EKL(M)
  220 CONTINUE
      IF (ONE.LE.ABS(DOTPJ)) GO TO 320
      IF (ONE.LE.ABS(DOTPK)) GO TO 320
      SINPJ = SQRT(ONE-DOTPJ*DOTPJ)
      SINPK = SQRT(ONE-DOTPK*DOTPK)
      DO 280 M = 1,3
         SMI = -CR1(M)/(DIJ*SINPJ*SINPJ)
         NOCOL1 = 3*(I-1)+M
         B(NOINT,NOCOL1) = SMI
         F1 = (CR1(M)*(DJK-DIJ*DOTPJ))/(DJK*DIJ*SINPJ*SINPJ)
         F2 = (DOTPK*CR2(M))/(DJK*SINPK*SINPK)
         SMJ = F1-F2
         NOCOL2 = 3*(J-1)+M
         B(NOINT,NOCOL2) = SMJ
         SML = CR2(M)/(DKL*SINPK*SINPK)
         NOCOL3 = 3*(L-1)+M
         B(NOINT,NOCOL3) = SML
         NOCOL4 = 3*(K-1)+M
         B(NOINT,NOCOL4) = -(SMI+SMJ+SML)
  280 CONTINUE
      DOT = ZERO
      DO 300 M = 1,3
  300 DOT = DOT+(CR1(M)*CR2(M))/(SINPJ*SINPK)
      IF (ABS(DOT) .GT. TOLRD) GO TO 320
         IF (ABS(DOT) .GT. ONE) DOT = SIGN(ONE,DOT)
         EQVAL =   ACOS(DOT)
         IF(ABS(EQVAL).LT.TOL) EQVAL=ZERO
         IF(ABS(EQVAL-PI).LT.TOL) EQVAL=PI
         SENSE = ddot(3,CR2,1,RIJ,1)
         IF(SENSE.LT.ZERO) EQVAL = -EQVAL
         RETURN
C
  320 CONTINUE
      WRITE (6,9008)
      CALL errquit('tors: error',911, UNKNOWN_ERR)
      STOP
 9008 FORMAT(31H ROUNDOFF ERROR IN TORS - STOP )
      END
