      SUBROUTINE vib_BSTR(EQVAL,NOINT,I,J,C,B,NDIM)
C*MODULE VIB  *DECK BSTR
* $Id$
C
      IMPLICIT NONE ! DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C          BOND STRETCH AS DEFINED BY WILSON (SEE WDC P.55) -----
C
      integer i, j, noint, ndim
      double precision C(3,*),B(NDIM,*)
      double precision RIJ(3)
      double precision eqval
C
      double precision zero
      PARAMETER (ZERO=0.0D+00)
C
      double precision dijsq
      integer m, nocol1, nocol2
C
      DIJSQ = ZERO
      DO 100 M = 1,3
         RIJ(M) = C(M,J)-C(M,I)
         DIJSQ = DIJSQ+RIJ(M)*RIJ(M)
  100 CONTINUE
      DO 120 M = 1,3
         NOCOL1 = 3*(I-1)+M
         NOCOL2 = 3*(J-1)+M
         B(NOINT,NOCOL1) = -RIJ(M)/SQRT(DIJSQ)
         B(NOINT,NOCOL2) = -B(NOINT,NOCOL1)
  120 CONTINUE
      EQVAL = SQRT(DIJSQ)
      RETURN
      END
