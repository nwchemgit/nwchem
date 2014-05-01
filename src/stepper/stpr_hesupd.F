C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_hesupd ( STEP, GRAD, HESIAN, SCR )
c $Id$
C
C     This subroutine updates the hessian using the hessian from the
C     old geometry and the difference in the coordinates and
C     gradients between the current and old geometries.
C
C     For finding a minimum on the surface (OMIN or STD =.TRUE.)
C     the BFGS  update is used; for TS=.TRUE. the Powell method
C     is chosen.
C     Note that hessian is a symmetric packed vector.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TOLER=1.0D-25)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
#include "chow.fh"
      DIMENSION STEP(NAT3), GRAD(NAT3), HESIAN(NAT3TR), SCR(NAT3)
      IF (OMIN.OR.STD) THEN
C
C       BFGS update
C       first, form vector hessian matrix times step vector
C       also, form dot products grad dot step and step dot
C       (hessian*step)
C
        DM1 = ZERO
        DM2 = ZERO
        DO 20 I = 1,NAT3
          SCR(I) = ZERO
          DO 10 J = 1,NAT3
            IJ = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
            SCR(I) = SCR(I) + HESIAN(IJ)*STEP(J)
   10     CONTINUE
          DM1 = DM1 + STEP(I)*GRAD(I)
          DM2 = DM2 + STEP(I)*SCR(I)
   20   CONTINUE
        IF (ABS(DM1) .LE. TOLER) THEN
          WRITE (ITAP06,1000) DM1
 1000     FORMAT ('  * ERROR * Dot product of STEP with',
     &               '  GRAD = ',1PD12.4)
          call errquit('stepper:hesupd fatal error',911, UNKNOWN_ERR)
        ENDIF
        IF (ABS(DM2) .LE. TOLER) THEN
          WRITE (ITAP06,2000) DM2
 2000     FORMAT ('  * ERROR * STEP*HESIAN*STEP = ',1PD12.4)
          call errquit('stepper:hesupd fatal error',911, UNKNOWN_ERR)
        ENDIF
        DM1 = ONE/DM1
        DM2 = ONE/DM2
C
C       update
C
        DO 40 I = 1,NAT3
          DO 30 J = 1,I
            IJ = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
            HESIAN(IJ) = HESIAN(IJ) +
     &                   DM1*GRAD(I)*GRAD(J) - DM2*SCR(I)*SCR(J)
   30     CONTINUE
   40   CONTINUE
      ELSE IF (TS) THEN
C
C       Powell update
C       First, form T vector = grad - hessian times step
C       Also form dot products T dot step and step dot step
C
        DM1 = ZERO
        DM2 = ZERO
        DO 60 I = 1,NAT3
          SCR(I) = ZERO
          DO 50 J = 1,NAT3
            IJ = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
            SCR(I) = SCR(I) + HESIAN(IJ)*STEP(J)
   50     CONTINUE
          SCR(I) = GRAD(I) - SCR(I)
          DM1 = DM1 + STEP(I)*STEP(I)
          DM2 = DM2 + STEP(I)*SCR(I)
   60   CONTINUE
        IF (ABS(DM1) .LE. TOLER) THEN
          WRITE (ITAP06,3000) DM1
 3000     FORMAT ('  * ERROR * norm of step vector = ',
     &    1PD12.4)
          call errquit('stepper:hesupd fatal error',911, UNKNOWN_ERR)
        ENDIF
        DM1 = ONE/DM1
C
C       update
C
        DO 80 I = 1,NAT3
          DO 70 J = 1,I
            IJ = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
            HESIAN(IJ) = HESIAN(IJ) + DM1*(SCR(I)*STEP(J) +
     &                   STEP(I)*SCR(J) - STEP(I)*STEP(J)*DM1*DM2)
   70     CONTINUE
   80   CONTINUE
      ELSE
        WRITE (ITAP06,*)' UNKNOWN TYPE OF CALCULATION IN HESUPD'
        call errquit('stepper:hesupd fatal error',911, INPUT_ERR)
      ENDIF
      RETURN
      END
C> @}
