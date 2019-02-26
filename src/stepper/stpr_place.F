C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_place(HESS,GRAD,COORD,ENERGY,IERR)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
#include "cstprfiles.fh"
      DIMENSION HESS(NAT3TR),
     &          GRAD(NAT3),
     &          COORD(3,NATOM)
      IERR = 0

#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP45, FORM='UNFORMATTED', FILE=FILE45,
     &   ACCESS='SEQUENTIAL', STATUS='UNKNOWN',err=10)
      REWIND (UNIT=ITAP45)
#else
      REWIND ITAP45
#endif
C
C     Write HESS.
C
*debug:      write(6,*)'hess stpr_place',nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
      CALL stpr_wrtsq (ITAP45,HESS,NAT3TR,IERR)
      IF(IERR.NE.0)GOTO 10
C
C     Write GRAD.
C
      CALL stpr_wrtsq (ITAP45,GRAD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 10
C
C     Write COORD.
C
      CALL stpr_wrtsq (ITAP45,COORD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 10
C
C     Write ENERGY.
C
      WRITE(ITAP45, ERR=10)ENERGY
#if defined(OPENCLOSE)
      close(unit=itap45,status='keep')
#endif
      RETURN
   10 CONTINUE
      IERR = 1
      RETURN
      end
      subroutine stpr_recall(HESS,GRAD,COORD,ENERGY,IERR)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
#include "cstprfiles.fh"
      DIMENSION HESS(NAT3TR),
     &          GRAD(NAT3),
     &          COORD(3,NATOM)
      IERR = 0
#if defined(OPENCLOSE)
      OPEN (UNIT=ITAP45, FORM='UNFORMATTED', FILE=FILE45,
     &   ACCESS='SEQUENTIAL', STATUS='UNKNOWN',err=20)
      REWIND (UNIT=ITAP45)
#else
      REWIND ITAP45
#endif
C
C     Read HESS.
C
      CALL stpr_readsq (ITAP45,HESS,NAT3TR,IERR)
      IF(IERR.NE.0)GOTO 20
*debug:      write(6,*)'hess stpr_recall', nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
C
C     Read GRAD.
C
      CALL stpr_readsq (ITAP45,GRAD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 20
C
C     Read COORD.
C
      CALL stpr_readsq (ITAP45,COORD,NAT3,IERR)
      IF(IERR.NE.0)GOTO 20
C
C     Read ENERGY.
C
      READ(ITAP45, END=20, ERR=20)ENERGY
#if defined(OPENCLOSE)
      close(unit=itap45,status='keep')
#endif
      RETURN
   20 CONTINUE
      IERR = 1
      RETURN
      END
C> @}
