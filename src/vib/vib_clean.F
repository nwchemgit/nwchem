      SUBROUTINE vib_CLEAN(A,LEN,THR)
* $Id: vib_clean.F,v 1.1 1995-07-29 23:19:24 d3e129 Exp $
C
C  This routine "cleans" or rids a vecter of the elements smaller than
C       the supplied threshold
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (D0=0.0D00)
      DIMENSION A(LEN)   ! Vector to be cleaned up
C
      DO 00100 II = 1,LEN             ! length of vector
         XDUM = ABS(A(II))            ! absolute value
         IF (XDUM.LT.THR) A(II) = D0  ! if element small zero
00100 CONTINUE
      RETURN
      END
