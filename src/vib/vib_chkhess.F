      SUBROUTINE vib_CHKHESS(HESS,NAT3)
* $Id: vib_chkhess.F,v 1.2 1995-10-21 09:16:04 d3e129 Exp $
C
C  This routine checks to make sure the analytic Hessian is symmetric
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
#include "stdio.fh"
      LOGICAL LSYM
      double precision thresh
      PARAMETER(THRESH = 1.0D-5)
      integer nat3
      double precision HESS(NAT3,NAT3) ! full square Hessian
c
      double precision xdiff
      integer iii,  jjj, jjjend
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
