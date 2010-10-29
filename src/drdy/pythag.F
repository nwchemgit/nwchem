C $Id$
      Double Precision Function drdy_pythag(a,b)
C***BEGIN PROLOGUE  PYTHAG
C***REFER TO  EISDOC
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
C***ROUTINES CALLED    (NONE)
C***END PROLOGUE  PYTHAG
      Implicit None
      DOUBLE PRECISION A,B
C
      DOUBLE PRECISION P,Q,R,S,T
C***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 drdy_PYTHAG = P
      RETURN
      END
