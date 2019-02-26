C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_SAVE(HESS,GRAD,COORD,ENERGY,IERR)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL SECO, SWSEC, SWFIR
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / PARAM / PSTEPP(20)
      COMMON / SWITCH / SECO, SWSEC, SWFIR
#include "cstprfiles.fh"
      DIMENSION HESS(NAT3TR),
     &          GRAD(NAT3),
     &          COORD(3,NATOM)
      IERR = 0
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN', err=10)
      REWIND (UNIT=ITAP41)
#else
      REWIND ITAP41
#endif
C
C     Write previous step parameters PSTEPP.
C
      WRITE(ITAP41, ERR=10)PSTEPP
C
C     Write HESS.
C
*debug:      write(6,*)'hess stpr_save',nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
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
#if defined(OPENCLOSE)
      close(unit=itap41,status='keep')
#endif
      RETURN
   10 CONTINUE
      IERR = 1
      RETURN
      end
      subroutine stpr_RETRVE(HESS,GRAD,COORD,ENERGY,IERR)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL SECO, SWSEC, SWFIR
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / PARAM / PSTEPP(20)
      COMMON / SWITCH / SECO, SWSEC, SWFIR
#include "cstprfiles.fh"
      DIMENSION HESS(NAT3TR),
     &          GRAD(NAT3),
     &          COORD(3,NATOM)
c
      IERR = 0
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN', err=20)
      REWIND (UNIT=ITAP41)
#else
      REWIND ITAP41
#endif
C
C     Read previous step parameters PSTEPP.
C
      READ(ITAP41, END=20, ERR=20)PSTEPP
C
C     Read HESS.
C
      CALL stpr_readsq (ITAP41,HESS,NAT3TR,IERR)
      IF(IERR.NE.0)GOTO 20
*debug:      write(6,*)'hess stpr_retrve',nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
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
#if defined(OPENCLOSE)
      close(unit=itap41,status='keep')
#endif
      RETURN
   20 CONTINUE
      IERR = 1
      RETURN
      end
      subroutine stpr_PUTP(GMXMUM,GLNGTH,SLNGTH,EOLD,ENERGY,
     &           EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &           E2NEW,NWALK,NNEG)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / PARAM / PSTEPP(20)
#include "cstprfiles.fh"
C
C     PSTEPP(3) currently not used.
C
      PSTEPP(1)  = GLNGTH
      PSTEPP(2)  = GMXMUM
      PSTEPP(4)  = SLNGTH
      PSTEPP(5)  = EOLD
      PSTEPP(6)  = ENERGY
      PSTEPP(7)  = EDIFF
      PSTEPP(8)  = E2OLD
      PSTEPP(9)  = ERATIO
      PSTEPP(10) = TRUSTO
      PSTEPP(11) = TRUST
      PSTEPP(12) = TRUSTN
      PSTEPP(13) = E2NEW
      PSTEPP(14) = dble(NWALK)
      PSTEPP(15) = dble(NNEG)
C
C     Write previous step parameters PSTEPP.
C
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE41,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN', err=30)
      REWIND (UNIT=ITAP41)
      WRITE(ITAP41)PSTEPP
00001 continue
      read(itap41,end=00030)
      goto 00001
00030 continue
      close(unit=itap41,status='keep')
#else
      REWIND ITAP41
      WRITE(ITAP41)PSTEPP
#endif
      RETURN
      end
      subroutine stpr_GETP(GMXMUM,GLNGTH,SLNGTH,EOLD,ENERGY,
     &           EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &           E2NEW,NWALK,NNEG)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / PARAM / PSTEPP(20)
#include "cstprfiles.fh"
C
C     If not first step read previous step parameters PSTEPP from
C     disk.
C
      IF(NWALK.GT.1)THEN
#if defined(OPENCLOSE)
        OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        REWIND (UNIT=ITAP41)
        READ(ITAP41)PSTEPP
        close(unit=itap41,status='keep')
#else
        REWIND ITAP41
        READ(ITAP41)PSTEPP
#endif
      ENDIF
C
C     PSTEPP(3) currently not used.
C
      GLNGTH = PSTEPP(1)
      GMXMUM = PSTEPP(2)
      SLNGTH = PSTEPP(4)
      EOLD   = PSTEPP(5)
      ENERGY = PSTEPP(6)
      EDIFF  = PSTEPP(7)
      E2OLD  = PSTEPP(8)
      ERATIO = PSTEPP(9)
      TRUSTO = PSTEPP(10)
      TRUST  = PSTEPP(11)
      TRUSTN = PSTEPP(12)
      E2NEW  = PSTEPP(13)
#if defined(CRAY) || defined(KSR)
      NWALK  = NINT(PSTEPP(14))
      NNEG   = NINT(PSTEPP(15))
#else
      NWALK = IDNINT(PSTEPP(14))
      NNEG  = IDNINT(PSTEPP(15))
#endif
      RETURN
      end
      subroutine stpr_GTNPO(EOLD,E2OLD,TRUST,SOLD,NWALK)
C
C     Read previous step parameters PSTEPP.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / PARAM / PSTEPP(20)
#include "cstprfiles.fh"
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &   ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP41)
      READ(ITAP41)PSTEPP
      close(unit=itap41,status='keep')
#else
      REWIND ITAP41
      READ(ITAP41)PSTEPP
#endif
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
      end
      subroutine stpr_GETNP(EOLD,E2OLD,TRUST,SOLD,NWALK)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / PARAM / PSTEPP(20)
#include "cstprfiles.fh"
C
C     Read previous step parameters PSTEPP.
C
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &    ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP41)
      READ(ITAP41)PSTEPP
      close(unit=itap41,status='keep')
#else
      REWIND ITAP41
      READ(ITAP41)PSTEPP
#endif
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
C> @}
