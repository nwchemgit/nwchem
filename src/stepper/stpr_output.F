      SUBROUTINE stpr_output(STEP,COORD,BCKSTP)
c $Id: stpr_output.F,v 1.2 1994-10-26 23:12:38 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL BCKSTP
      LOGICAL SECO, SWSEC, SWFIR
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / SWITCH / SECO, SWSEC, SWFIR
      DIMENSION STEP(NAT3),
     &          COORD(3,NATOM)
      CALL stpr_getp(GLNGTH,SLNGTH,EOLD,ENERGY,
     &           EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &           E2NEW,NWALK,NNEG)
      IF ((ICALC.EQ.2).AND.(SECO.AND.(IWCTR.EQ.1))) THEN
        WRITE(6,*)' STEPPER| Stepping switched to 1-st derivatives'
      ENDIF
      IF ((ICALC.EQ.2).AND.SWFIR) THEN
        WRITE(6,*)' STEPPER| Using Analytical 1st derivatives. '
      ENDIF
      IF ((ICALC.EQ.1).AND.SWSEC) THEN
        WRITE(6,*)' STEPPER| Using Analytical 1st and 2nd derivatives. '
      ENDIF
      WRITE(6,1000)NWALK
      WRITE(6,1006)EOLD
      WRITE(6,1007)ENERGY
      WRITE(6,1008)EDIFF
      WRITE(6,1009)E2OLD
      WRITE(6,1010)ERATIO
      IF(BCKSTP)THEN
        WRITE(6,1012)TRUST
      ELSE
        WRITE(6,1011)TRUSTO
      ENDIF
      WRITE(6,1001)NNEG
      WRITE(6,1002)GLNGTH
#if defined(KSR) || defined(CRAY)
      GLNGTH = GLNGTH/SQRT(real(NATOM))
#else
      GLNGTH = GLNGTH/SQRT(DBLE(NATOM))
#endif
      WRITE(6,1022)GLNGTH
      WRITE(6,1005)SLNGTH
      WRITE(6,1014)E2NEW
C
C     Write old coordinates.
C
      WRITE(6,1015)
      WRITE(6,1016)
      DO 5 I = 1,NATOM
        WRITE(6,1017)I,(COORD(J,I),J=1,3)
    5 CONTINUE
C
C     Write step.
C
      WRITE(6,1018)
      WRITE(6,1019)
      DO 6 I = 1,NATOM
        JB = (I-1)*3 + 1
        JE = I*3
        WRITE(6,1017)I,(STEP(J),J=JB,JE)
    6 CONTINUE
      RETURN
 1000 FORMAT(' STEPPER|        Walk statistics for step: ',I3)
 1001 FORMAT(' STEPPER|  Number of negative eigenvalues: ',I3)
 1002 FORMAT(' STEPPER|                   Gradient norm: ',D20.12)
 1022 FORMAT(' STEPPER|               RMS Gradient norm: ',D20.12)
 1003 FORMAT(' STEPPER|                 Value of lambda: ',D20.12)
 1004 FORMAT(' STEPPER|         Value of scale factor A: ',D20.12)
 1005 FORMAT(' STEPPER| Largest component of step v-tor: ',D20.12)
 1006 FORMAT(' STEPPER|    Previous actual total energy: ',D20.12)
 1007 FORMAT(' STEPPER|     Current actual total energy: ',D20.12)
 1008 FORMAT(' STEPPER|        Actual energy difference: ',D20.12)
 1009 FORMAT(' STEPPER|   Second order energy predicted: ',D20.12)
 1010 FORMAT(' STEPPER|            Ratio of differences: ',D20.12)
 1011 FORMAT(' STEPPER|              Fixed trust radius: ',D20.12)
 1012 FORMAT(' STEPPER|  Modified trust radius for step: ',D20.12)
 1014 FORMAT(' STEPPER| New 2nd order energy prediction: ',D20.12)
 1015 FORMAT(' STEPPER| ',18('-'),' Input coordinates ',18('-'))
 1016 FORMAT(' STEPPER| ATOM:',12x,'X',14x,'Y',14x,'Z')
 1017 FORMAT(' STEPPER| ',I3,6X,3(1X,F14.10))
 1018 FORMAT(' STEPPER| ',20('-'),' Current step ',21('-'))
 1019 FORMAT(' STEPPER| STEP:',12x,'X',14x,'Y',14x,'Z')
 1020 FORMAT(' STEPPER| ',19('-'),' New coordinates ',19('-'))
      END
