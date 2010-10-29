      SUBROUTINE vib_CLEAN(A,LEN,THR)
* $Id$
C
C  This routine "cleans" or rids a vecter of the elements smaller than
C       the supplied threshold
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
      DOUBLE PRECISION D0
      PARAMETER (D0=0.0D00)
      INTEGER LEN
      DOUBLE PRECISION A(LEN) ! Vector to be cleaned up
      DOUBLE PRECISION THR
C
      DOUBLE PRECISION XDUM
      INTEGER II
C
      DO 00100 II = 1,LEN             ! length of vector
         XDUM = ABS(A(II))            ! absolute value
         IF (XDUM.LT.THR) A(II) = D0  ! if element small zero
00100 CONTINUE
      RETURN
      END
