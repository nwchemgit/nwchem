c $Id$
      SUBROUTINE hess_hssout(V,M,N,NDIM)
      IMPLICIT none
#include "stdio.fh"
c
      double precision v
      integer m, n, ndim
      DIMENSION V(NDIM,1)
c
      integer max, imax, imin, i, j
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C     This is the hnd_hssout routine -> nwchemized
C
      MAX=7
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (LuOut,9008)
      WRITE (LuOut,9128) (I,I = IMIN,IMAX)
      WRITE (LuOut,9008)
      DO 120 J = 1,N
      WRITE (LuOut,9148) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
 9008 FORMAT(/)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.7)
      END
