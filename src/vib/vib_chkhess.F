      SUBROUTINE vib_CHKHESS(HESS,NAT3)
* $Id: vib_chkhess.F,v 1.1 1995-07-29 23:19:23 d3e129 Exp $
C
C  This routine checks to make sure the analytic Hessian is symmetric
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "stdio.fh"
      LOGICAL LSYM
      PARAMETER(THRESH = 1.0D-5)
      DIMENSION HESS(NAT3,NAT3) ! full square Hessian
C
      DO 00100 III = 2,NAT3  ! loop over size-1
         JJJEND = III-1
         LSYM = .TRUE.
         DO 00100 JJJ = 1,JJJEND
            XDIFF = HESS(III,JJJ)-HESS(JJJ,III)  ! difference of ij and ji
            XDIFF = ABS(XDIFF)
            IF(XDIFF.GT.THRESH)THEN ! is it symmetric??
               WRITE(luout,*)
               WRITE(luout,*)' **** HESSIAN NOT SYMMETRIC ****'
               WRITE(luout,*)' HESS(',III,',',JJJ,') =',HESS(III,JJJ)
               WRITE(luout,*)' HESS(',JJJ,',',III,') =',HESS(JJJ,III)
               WRITE(luout,*)
               LSYM = .FALSE.
            ENDIF
00100 CONTINUE
C
      IF (LSYM) THEN
        WRITE(luout,*)
        WRITE(luout,*)' Analytic Hessian passed symmetry test '
        WRITE(luout,*)
      END IF
      RETURN
      END
