      SUBROUTINE stpr_SAVE(HESS,GRAD,COORD,ENERGY,IERR)
c $Id: stpr_save.F,v 1.1 1994-06-24 20:42:42 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL SECO, SWSEC, SWFIR
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / PARAM / PSTEPP(20)
      COMMON / SWITCH / SECO, SWSEC, SWFIR
      DIMENSION HESS(NAT3TR),
     &          GRAD(NAT3),
     &          COORD(3,NATOM)
      IERR = 0
      REWIND ITAP41
C
C     Write previous step parameters PSTEPP.
C
      WRITE(ITAP41, ERR=10)PSTEPP
C
C     Write HESS.
C
      CALL stpr_wrtsq (ITAP41,HESS,NAT3TR,IERR)
      IF(IERR.NE.0)GOTO 10
C
C     Write GRAD.
C
      CALL stpr_wrtsq (ITAP41,GRAD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 10
C
C     Write COORD.
C
      CALL stpr_wrtsq (ITAP41,COORD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 10
C
C     Write ENERGY.
C
      WRITE(ITAP41, ERR=10)ENERGY
      RETURN
   10 CONTINUE
      IERR = 1
      RETURN
      ENTRY stpr_RETRVE(HESS,GRAD,COORD,ENERGY,IERR)
      IERR = 0
      REWIND ITAP41
C
C     Read previous step parameters PSTEPP.
C
      READ(ITAP41, END=20, ERR=20)PSTEPP
C
C     Read HESS.
C
      CALL stpr_readsq (ITAP41,HESS,NAT3TR,IERR)
      IF(IERR.NE.0)GOTO 20
C
C     Read GRAD.
C
      CALL stpr_readsq (ITAP41,GRAD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 20
C
C     Read COORD.
C
      CALL stpr_readsq (ITAP41,COORD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 20
C
C     Read ENERGY.
C
      READ(ITAP41, END=20, ERR=20)ENERGY
      RETURN
   20 CONTINUE
      IERR = 1
      RETURN
      ENTRY stpr_PUTP(GLNGTH,SLNGTH,EOLD,ENERGY,
     &           EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &           E2NEW,NWALK,NNEG)
C
C     PSTEPP(2) and PSTEPP(3) currently not used.
C
      PSTEPP(1) = GLNGTH
      PSTEPP(4) = SLNGTH
      PSTEPP(5) = EOLD
      PSTEPP(6) = ENERGY
      PSTEPP(7) = EDIFF
      PSTEPP(8) = E2OLD
      PSTEPP(9) = ERATIO
      PSTEPP(10) = TRUSTO
      PSTEPP(11) = TRUST
      PSTEPP(12) = TRUSTN
      PSTEPP(13) = E2NEW
#if defined(CRAY) || defined(KSR)
      PSTEPP(14) = real(NWALK)
      PSTEPP(15) = real(NNEG)
#else
      PSTEPP(14) = dble(NWALK)
      PSTEPP(15) = dble(NNEG)
#endif
C
C     Write previous step parameters PSTEPP.
C
      REWIND ITAP41
      WRITE(ITAP41)PSTEPP
      RETURN
      ENTRY stpr_GETP(GLNGTH,SLNGTH,EOLD,ENERGY,
     &           EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &           E2NEW,NWALK,NNEG)
C
C     If not first step read previous step parameters PSTEPP from
C     disk.
C
      IF(NWALK.GT.1)THEN
        REWIND ITAP41
        READ(ITAP41)PSTEPP
      ENDIF
C
C     PSTEPP(2) and PSTEPP(3) currently not used.
C
      GLNGTH = PSTEPP(1)
      SLNGTH = PSTEPP(4)
      EOLD = PSTEPP(5)
      ENERGY = PSTEPP(6)
      EDIFF = PSTEPP(7)
      E2OLD = PSTEPP(8)
      ERATIO = PSTEPP(9)
      TRUSTO = PSTEPP(10)
      TRUST = PSTEPP(11)
      TRUSTN = PSTEPP(12)
      E2NEW = PSTEPP(13)
#if defined(CRAY) || defined(KSR)
      NWALK = NINT(PSTEPP(14))
      NNEG  = NINT(PSTEPP(15))
#else
      NWALK = IDNINT(PSTEPP(14))
      NNEG  = IDNINT(PSTEPP(15))
#endif
      RETURN
      ENTRY stpr_GTNPO(EOLD,E2OLD,TRUST,SOLD,NWALK)
C
C     Read previous step parameters PSTEPP.
C
      REWIND ITAP41
      READ(ITAP41)PSTEPP
      EOLD = PSTEPP(5)
      E2OLD = PSTEPP(8)
      TRUST = PSTEPP(11)
      SOLD = PSTEPP(4)
#if defined(CRAY) || defined(KSR)
      NWALK = NINT(PSTEPP(14))
#else
      NWALK = IDNINT(PSTEPP(14))
#endif
      RETURN
      ENTRY stpr_GETNP(EOLD,E2OLD,TRUST,SOLD,NWALK)
C
C     Read previous step parameters PSTEPP.
C
      REWIND ITAP41
      READ(ITAP41)PSTEPP
      EOLD = PSTEPP(6)
      E2OLD = PSTEPP(13)
      TRUST = PSTEPP(12)
      SOLD = PSTEPP(4)
#if defined(CRAY) || defined(KSR)
      NWALK = NINT(PSTEPP(14))
#else
      NWALK = IDNINT(PSTEPP(14))
#endif
      RETURN
      END
