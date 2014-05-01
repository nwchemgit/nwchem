C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_wlkchk(ENEW,CONVGE,BCKSTP,NMODE)
c $Id$
C
C     Determination of walk information.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "util.fh"
      LOGICAL BCKSTP, LPRINT
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / OLDES / ELAST,E2PRED
      BCKSTP = .FALSE.
      LPRINT = util_print('line information', print_low)
      IF(ICALC.GT.0)THEN
        IF(IWCTR.LE.1)THEN
c          IF(LPRINT)THEN
c            WRITE(6,*)'  There is not yet enough information '
c            WRITE(6,*)'  available to determine geometry walk '
c            WRITE(6,*)'  characteristics. '
c          ENDIF
          RETURN
        ELSE
c          IF(LPRINT)THEN
c            WRITE(6,*)'  Checking geometry walk information. '
c          ENDIF
        ENDIF
      ELSE
        RETURN
      ENDIF
C
C     Calculate actual energy difference ... EDIFF
C
      EDIFF = ENEW - ELAST
C
C     If actual or 2nd order predicted energies are less then
C     an energy tolerance e.g. ETOL then DO NOT worry about
C     backstepping.
C
      ETOL = CONVGE
      IF(ABS(EDIFF).LT.ETOL.OR.ABS(E2PRED).LT.ETOL)THEN
        IF(LPRINT)THEN
          WRITE(6,*)'  Actual energy difference is : ',EDIFF
          WRITE(6,*)'  2nd order predicted energy difference ',
     &              ' was : ',E2PRED
          WRITE(6,*)'  One of these is below the energy ',
     &              ' tolerance: ',ETOL
          WRITE(6,*)'  Backstepping possibilities are NO '
          WRITE(6,*)'  longer considered appropriate. '
        ENDIF
C
C       Set backstep flag .FALSE.
C
        BCKSTP = .FALSE.
        ERATIO = 0.0D0
        RETURN
      ELSE
C
C       Calculate actual to predicted energy ratio.
C
        ENUM = MIN(ABS(EDIFF),ABS(E2PRED))
        EDEN = MAX(ABS(EDIFF),ABS(E2PRED))
        ERATIO = 0.0D0
        IF(ABS(EDEN).GT.1.D-12)ERATIO = ENUM/EDEN
      ENDIF
C
C     Check sign of actual to predicted energies.
C
      ISEDIF = DNINT(ABS(EDIFF)/EDIFF)
      IF(E2PRED.NE.0)THEN
        ISE2OL = DNINT(ABS(E2PRED)/E2PRED)
      ELSE
        ISE2OL = ISEDIF
      ENDIF
      IF(ICALC.EQ.2)THEN
        IF(ISEDIF.NE.ISE2OL)THEN
          IF(LPRINT)THEN
            WRITE(6,*)'  Signs : ',ISEDIF,ISE2OL
            WRITE(6,*)'  Actual energy difference is : ',EDIFF
            WRITE(6,*)'  2nd order predicted energy ',
     &                ' difference was : ',E2PRED
            WRITE(6,*)'  Actual and predicted energies have '
            WRITE(6,*)'  sign disagreement.  Back-step is in '
            WRITE(6,*)'  order here. '
          ENDIF
          BCKSTP = .TRUE.
        ENDIF
      ELSEIF(ICALC.EQ.1.AND.IWCTR.GT.2)THEN
        IF(ISEDIF.NE.ISE2OL)THEN
          IF(LPRINT)THEN
            WRITE(6,*)'  Signs : ',ISEDIF,ISE2OL
            WRITE(6,*)'  Actual energy difference is : ',EDIFF
            WRITE(6,*)'  2nd order predicted energy ',
     &                ' difference was : ',E2PRED
            WRITE(6,*)'  Actual and predicted energies have '
            WRITE(6,*)'  sign disagreement.  Back-step is in '
            WRITE(6,*)'  order here. '
          ENDIF
          BCKSTP = .TRUE.
        ENDIF
      ENDIF
      IF(BCKSTP)RETURN
C
C     Set energy ratio backstep criterion.
C
      IF(ICALC.EQ.2)THEN
        ERATM = 0.0D0
      ELSEIF(ICALC.EQ.1)THEN
        ERATM = 0.0D0
      ENDIF
C
C     Check ratio of 2nd order energy prediction to actual energy.
C
      IF(LPRINT)THEN
        WRITE(6,'(a,f10.6)')
     $        '  Energy ratio, actual to predicted = ',
     &            ERATIO
      ENDIF
      IF(ERATIO.LT.ERATM)THEN
          IF(LPRINT)THEN
            WRITE(6,*)'  Actual energy difference is : ',EDIFF
            WRITE(6,*)'  2nd order predicted energy ',
     &                ' difference was : ',E2PRED
          ENDIF
        IF(NMODE.GT.0)THEN
          IF(LPRINT)THEN
            WRITE(6,*)'  There is poor agreement between '
            WRITE(6,*)'  actual and predicted energies ... '
            WRITE(6,*)'  less than: ',ERATM,' percent.  '
            WRITE(6,*)'  Conservatism in transition state '
            WRITE(6,*)'  walks --> Back-step is in order here. '
          ENDIF
          BCKSTP = .TRUE.
        ELSE
          IF(LPRINT)THEN
            WRITE(6,*)'  There is poor agreement between '
            WRITE(6,*)'  actual and predicted energies ... '
            WRITE(6,*)'  less than: ',ERATM,' percent.  '
            WRITE(6,*)'  Conservatism suggests that a  '
            WRITE(6,*)'  Back-step is in order here. '
          ENDIF
          BCKSTP = .TRUE.
        ENDIF
      ENDIF
C
C     Write walk backstep flag.
C
      IF(BCKSTP)THEN
          IF(LPRINT)THEN
            WRITE(6,*)' >>>> Geometry backstep in progress <<<<'
          ENDIF
      ELSE
c          IF(LPRINT)THEN
c            WRITE(6,*)'  Predicted versus actual energy in '
c            WRITE(6,*)'  geometry walk is JUST FINE. '
c          ENDIF
      ENDIF
      RETURN
      END
C> @}
