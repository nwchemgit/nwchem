C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_driver(
     &       HESS, GRAD, STEP, VC,
     &       COORD, COOROLD, HIEIGS, GISR,
     &       CHARGE, SELECT, SCRSYM, SCRSQ,
     &       VECSQ, ATMASS, P, HMP,
     &       CMASS, TENIN, ENERGY, E2NEW, TRUST,
     &       CONVGE, CONVGG, CONVGGM,
     &       NMODE, TRACK, RPATH, BCKSTP, CNVGRD, rtdb )
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      LOGICAL TRACK,BCKSTP,CNVGRD,RPATH
#include "chow.fh"
      integer rtdb
      LOGICAL SECO, SWSEC, SWFIR
      integer idamax
      external idamax
      PARAMETER( ZERO=0.D0, ONE=1.D0 )
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / SWITCH / SECO, SWSEC, SWFIR
      DIMENSION HESS(NAT3TR),GRAD(NAT3),STEP(NAT3),VC(NAT3,NAT3),
     &          COORD(3,NATOM),COOROLD(3,NATOM),HIEIGS(NAT3),
     &          GISR(NAT3),CHARGE(NATOM),
     &          SELECT(NAT3),SCRSYM(NAT3TR),SCRSQ(NAT3,NAT3),
     &          VECSQ(NAT3,NAT3),ATMASS(NATOM),P(NAT3TR),
     &          HMP(NAT3TR),CMASS(3),TENIN(3,3)
      IF(BCKSTP)THEN
        WRITE(6,*)'  Retrieving old hessian, gradient '
        WRITE(6,*)'  coordinates, and energy. '
        IERR = 0
        CALL stpr_retrve(HESS,GRAD,COORD,ENERGY,IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'  Unsuccessful "RETRVE" attempt. '
          STOP
        ENDIF
      ENDIF

*      grad_max = dabsmax(nat3,grad)
      grad_max = abs(grad(idamax(nat3,grad,1)))

      NINTER = NAT3
C
C     Zero all arrays with the exception of the
C     incoming hessian; HESS, and gradient; GRAD.
C
      call dfill(NAT3,0.0d00,STEP,1)
      call dfill(NAT3,0.0d00,GISR,1)
      call dfill(NAT3,0.0d00,HIEIGS,1)
      call dfill(NAT3TR,0.0d00,SCRSYM,1)
      call dfill(NAT3SQ,0.0d00,SCRSQ,1)
C
C     Place necessary current iteration information: HESS, GRAD,
C     COORD, and ENERGY in cold storage (scratch disk) to recall
C     when needed.
C
      IERR = 0
      CALL stpr_place(HESS,GRAD,COORD,ENERGY,IERR)
      IF(IERR.NE.0)THEN
        WRITE(6,*)'  Unsuccessful "PLACE" attempt. '
        STOP
      ENDIF
C
C     Partition hessian.
C
      CALL stpr_partit(rtdb,
     $     HESS, GRAD, COORD, HIEIGS, GISR, VC, SCRSYM,
     &            SCRSQ, VECSQ, STEP, SELECT, P, HMP,
     &            ATMASS, CMASS, TENIN, ENERGY,
     &            CONVGE, CONVGG, CONVGGM,
     &            NINTER, NMODE, RPATH, TRACK, BCKSTP)
*      write(6,*)' step after partit '
*      do jjjj = 1,nat3
*      write(6,'(1x,i5,1x,f20.10)')jjjj,step(jjjj)
*      enddo
      IF(NINTER.NE.0)THEN
        CALL stpr_mkstep(HIEIGS,GISR,VC,STEP,SELECT,SCRSQ,ATMASS,ENERGY,
     &              CONVGE, CONVGG, CONVGGM, GRAD_MAX,
     &              E2NEW,NINTER,NAT3,NMODE,TRUST,BCKSTP,CNVGRD,
     &              rtdb)
*      write(6,*)' step after mkstep '
*      do jjjj = 1,nat3
*      write(6,'(1x,i5,1x,f20.10)')jjjj,step(jjjj)
*      enddo
      ELSE
       CNVGRD = .TRUE.
       WRITE(6,*)'  ENERGY = ',ENERGY
      ENDIF
      IF (.NOT.CNVGRD) THEN
        CALL stpr_output(STEP,COORD,BCKSTP, grad_max,
     &                   CONVGGM, CONVGG, CONVGE )
      ELSE
c
c     Write out useful hessian info:
c

      ENDIF
C
C     Restore currently used information and write to backstep disk.
C
      IERR = 0
      CALL stpr_recall(HESS,GRAD,COORD,ENERGY,IERR)
      IF(IERR.NE.0)THEN
        WRITE(6,*)'  Unsuccessful "RECALL" attempt. '
        STOP
      ENDIF
      IERR = 0
      CALL stpr_save(HESS,GRAD,COORD,ENERGY,IERR)
      IF(IERR.NE.0)THEN
        WRITE(6,*)'  Unsuccessful "SAVE" attempt. '
        STOP
      ENDIF
      CALL stpr_sumstc(STEP,COORD,ATMASS,CMASS,TENIN,CNVGRD)
      RETURN
      END
C> @}
