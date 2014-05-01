C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_chkhess(HESS)
c $Id$
C
C     This routine checks to make sure the analytic Hessian
C     is symmetric.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      LOGICAL LSYM
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION HESS(NAT3,NAT3)
      DATA SYMCHK/1.D-6/
C
      DO 00100 I = 2,NAT3
         LSYM = .TRUE.
         DO 00100 J = 1,I-1
            XDIFF = HESS(I,J)-HESS(J,I)
            XDIFF = DABS(XDIFF)
            IF(XDIFF.GT.SYMCHK)THEN
               WRITE(6,*)
               WRITE(6,*)'  **** Hessian not symmetric ****'
               WRITE(6,*)'  HESS(',I,',',J,') =',HESS(I,J)
               WRITE(6,*)'  HESS(',J,',',I,') =',HESS(J,I)
               WRITE(6,*)
               LSYM = .FALSE.
            ENDIF
00100 CONTINUE
C
      IF (LSYM) THEN
        WRITE(6,*)'  Analytic Hessian passed symmetry test. '
      else
        call errquit('stepper:chkhess fatal error',911, UNKNOWN_ERR)
      END IF
      RETURN
      END
C> @}
