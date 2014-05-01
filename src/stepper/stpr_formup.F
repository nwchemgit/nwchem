C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_formup(HESS,GRADC,GRADD,COORDC,COORDD,SCR,rtdb)
c $Id$
C
C     HESS = hessian to be updated
C     GRADC = current gradient
C     GRADD = difference of gradients
C     COORDC = current coordinates
C     COORDD = difference of coordinates (STEP)
C     SCR = scratch
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / PARAM / PSTEPP(20)
#include "cstprfiles.fh"
      integer rtdb
      DIMENSION HESS(NAT3TR),
     &          GRADC(NAT3),GRADD(NAT3),
     &          SCR(NAT3),
     &          COORDC(3,NATOM),COORDD(3,NATOM)
      logical does_it_exist
*debug:      write(6,*)' stpr_formup: iwctr =',iwctr
      IF(IWCTR.EQ.1)THEN
C
C       Attempt to retrieve any "old" GRAD, HESS, etc.
C
        IERR = 0
        CALL stpr_retrve(HESS,GRADD,COORDD,ENERGY,IERR)
        IF(IERR.EQ.0)THEN
          WRITE(6,*)'  Nuclear Hessian retrieved from',
     &              ' previous run. '
          GOTO 9
        ENDIF
C
C       Try a user defined ASCII file.
C
        does_it_exist = .false.
        inquire(file=FILEHESS,exist=does_it_exist)
        if (does_it_exist) then
          OPEN (UNIT=46, FILE=FILEHESS, FORM='FORMATTED',
     &        ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
          REWIND (UNIT=46)
          READ(46,*,END=6,ERR=6)(HESS(I),I=1,NAT3TR)
          CLOSE (UNIT=46,STATUS='KEEP')
          WRITE(6,*)
     &        '  Nuclear Hessian retrieved from ASCII file.'
          RETURN
        endif
    6   CONTINUE
C
C       Last resort, use unit hessian.
C
        CALL stpr_inihes(HESS,rtdb)
        RETURN
      ELSE
        IERR = 0
        CALL stpr_retrve(HESS,GRADD,COORDD,ENERGY,IERR)
        IF(IERR.NE.0)
     &         call errquit('stepper:formup fatal error',911,
     &       UNKNOWN_ERR)
      ENDIF
*debug:      write(6,*)'hess retrieved', nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
    9 CONTINUE
      CC = 0.0D0
      DO 20 J = 1,NATOM
        DO 10 I = 1,3
          COORDD(I,J) = COORDC(I,J) - COORDD(I,J)
          CC = CC + DABS(COORDD(I,J))
   10   CONTINUE
   20 CONTINUE
      DO 30 K = 1,NAT3
        GRADD(K) = GRADC(K) - GRADD(K)
   30 CONTINUE
C
C     Check to see if COORDD .ne. 0 to update
C
      IF(CC.GT.1.0D-10)THEN
        CALL stpr_hesupd ( COORDD, GRADD, HESS, SCR )
      ELSE
        WRITE(6,*)'  Negligible coordinate change. No Update. '
      ENDIF
*debug:      write(6,*)'hess updated', nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
      RETURN
      END
C> @}
