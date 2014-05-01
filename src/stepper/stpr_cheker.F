C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_cheker(
     &       ENERGY,EOLD,E2OLD,EDIFF,
     &       GMXMUM,GLNGTH,ERATIO,
     &       CONVGE, CONVGG, CONVGGM,
     &       TRUSTI,TRUST,TRUSTO,SOLD,NNEG,NWALK,
     &       BCKSTP,CNVGRD)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      LOGICAL BCKSTP,CNVGRD
#include "chow.fh"
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
C
      IF (MWC) THEN
        CONVGG=1.0D-08
      ENDIF
C
C     Set energy ratio backstep criterion.
C
      IF(ICALC.EQ.2)THEN
        ERATM = 0.5D0
      ELSEIF(ICALC.EQ.1)THEN
        ERATM = 0.0D0
      ELSE
        WRITE(6,*)'  Unknown calculation type. Program aborts. '
        call errquit('stepper:cheker fatal error',911, INPUT_ERR)
      ENDIF
C
C     Retrieve previous step parameters (if not first step).
C
      IF(IWCTR.GT.1)THEN
        IF(BCKSTP)THEN
          CALL stpr_gtnpo(EOLD,E2OLD,TRUST,SOLD,NWALK)
        ELSE
          CALL stpr_getnp(EOLD,E2OLD,TRUST,SOLD,NWALK)
        ENDIF
      ELSE	
         E2OLD=0.d0
         EOLD=0.d0
      ENDIF
C
C     Initialize current step parameters.
C
      NWALK = NWALK + 1
C
C     Calculate actual energy difference ... EDIFF.
C
      EDIFF = ENERGY - EOLD
C
C     Calculate actual to predicted energy ratio.
C
      ERATIO = MIN(ABS(EDIFF),ABS(E2OLD))/
     &         MAX(ABS(EDIFF),ABS(E2OLD))
C
C     Check for convergence ... look at:
C       1. Structure of HESS (proper number of negative eigenvalues),
C       2. Energy convergence, and
C       3. Gradient convergence.
C
      write(6,9876)
      write(6,9875) ediff,glngth,gmxmum
      write(6,9874) convge,convgg,convggm
 9876 format('               e-diff       ','  g-norm  ','   g-max   ')
 9875 format(' value',f20.10,f10.6,f10.6)
 9874 format(' thres',f20.10,f10.6,f10.6)
      IF((NNEG.EQ.0.AND..NOT.TS).OR.(NNEG.EQ.1.AND.TS))THEN
        IF(ABS(EDIFF).LT.CONVGE.AND.GLNGTH.LT.CONVGG
     &                         .AND.GMXMUM.LT.CONVGGM)THEN
          WRITE(6,*)'  ENERGY = ',ENERGY
          WRITE(6,*)'  The hessian exhibits the "proper" '
          WRITE(6,*)'  structure, and the gradient and energy '
          WRITE(6,*)'  convergence criterion has been met.  '
          WRITE(6,*)'  This geometry walk is HISTORY.'
          CNVGRD = .TRUE.
          TRUST = TRUSTI
          TRUSTO = TRUSTI
          RETURN
        ELSE
          CNVGRD = .FALSE.
        ENDIF
      ELSE
        IF(ABS(EDIFF).LT.CONVGE.AND.GLNGTH.LT.CONVGG
     &                         .AND.GMXMUM.LT.CONVGGM)THEN
          WRITE(6,*)'  ENERGY = ',ENERGY
          WRITE(6,*)'  The gradient and energy convergence '
          WRITE(6,*)'  criterion has been met, BUT the hessian '
          WRITE(6,*)'  does not exhibit proper structure. '
          WRITE(6,*)'  There are: ',NNEG,' negative hessian '
          WRITE(6,*)'  eigenvalues. '
          CNVGRD = .TRUE.
          TRUST = TRUSTI
          TRUSTO = TRUSTI
          RETURN
        ELSE
          CNVGRD = .FALSE.
        ENDIF
      ENDIF
C
C     Set TRUST radius and check for backstepping.  Use algorithm
C     suitable for analytical versus update hessian.
C
      IF (ICALC.GT.0) THEN
        IF(NWALK.EQ.1)THEN
          TRUST = TRUSTI
          TRUSTO = TRUSTI
        ELSE
          IF(BCKSTP)THEN
            WRITE(6,*)'  Backstep in progress; trust radius ',
     $            'set to one-quarter of previous step length.'
            TRUSTO = TRUST
            TRUST = 0.25d0*SOLD
          ELSE
            TRUST = TRUSTI
            TRUSTO = TRUSTI
          ENDIF
            write(6,'(a,f8.2)') '  trust radius = ',trust
        ENDIF
      ELSE
        WRITE(6,*)'  Unknown calculation type. Program aborts. '
        call errquit('stepper:cheker fatal error',911, INPUT_ERR)
      ENDIF
C
C     Special Case: If on initial iteration there are no
C     significant forces then stop.
C
      IF(NWALK.EQ.1.AND.GLNGTH.LT.CONVGG
     &             .AND.GMXMUM.LT.CONVGGM)THEN
        WRITE(6,*)'  >>> W A R N I N G     W A R N I N G <<< '
        WRITE(6,*)'  >>> Negligible gradient on first',
     &            ' iteration. <<< '
        WRITE(6,*)'  >>> Arbitrary step could be taken. <<< '
        WRITE(6,*)'  >>> Nudge coordinates away from',
     &            ' minimum and restart. <<< '
        call errquit('stepper:cheker fatal error',911, INPUT_ERR)
      ENDIF
      RETURN
      END
C> @}
