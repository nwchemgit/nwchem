C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_partit(rtdb,
     &       HESS, GRAD, COORD, HIEIGS, GISR, VC, SCRSYM,
     &       SCRSQ, VECSQ, STEP, SELECT, P, HMP,
     &       ATMASS, CMASS, TENIN, ENERGY,
     &       CONVGE, CONVGG, CONVGGM,
     &       NINTER, NMODE, RPATH, TRACK, BCKSTP)
c $Id$
c
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      LOGICAL TRACK, BCKSTP, RPATH, TROUB
      LOGICAL COOK, LPRINT, oprint
      PARAMETER( ZERO=0.D0, ONE=1.D0 )
#include "rtdb.fh"
#include "util.fh"
#include "chow.fh"
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / THRESH / RCONVG,RSTD
      COMMON / MASS / TOTM, NUMAS
      DIMENSION HESS(NAT3TR), GRAD(NAT3), COORD(3,NATOM),
     &          HIEIGS(NAT3), GISR(NAT3), VC(NAT3,NAT3),
     &          SCRSYM(NAT3TR), SCRSQ(NAT3,NAT3),
     &          VECSQ(NAT3,NAT3), STEP(NAT3), SELECT(NAT3),
     &          P(NAT3TR), HMP(NAT3TR), ATMASS(NATOM),
     &          CMASS(3), TENIN(3,3)
      DIMENSION UNIVEC(3), SQMAS(3)
      character*26 date
      integer rtdb
C
C     Internal hessian/gradient symmetry tolerance.
C
      DATA SYMCHK/1.D-8/
c
      ISYM2(I,J)=MAX(I,J)*((MAX(I,J))-1)/2 + MIN(I,J)
c
      LPRINT = util_print('partition', print_debug)
      oprint = util_print('information', print_default)
c
*debug:      call stpr_datestp(6,'partit: start')
      AUCM = 2.194746D+05
*debug:      LPRINT = .TRUE.
      IF (MWC) THEN
        SYMCHK=1.0D-08
      ENDIF
C
C     RCONVG can not be more restrictive than SYMCHK.
C
      IF (DABS(RCONVG).LT.DABS(SYMCHK)) RCONVG = SYMCHK
C
C     Initialize TROUB.
C
      TROUB = .FALSE.
C
C     Initialize COOK.
C
      COOK = .FALSE.
C
C     Tricky stuff.
C
      IF (MWC) THEN
        COOK = .TRUE.
        MWC = .FALSE.
      ENDIF
 999  CONTINUE
C
C     Zero incoming scratch arrays, and the output arrays:
C     HIEIGS, GISR, and the VC array.
C
      call dfill(nat3sq,0.0d00,vecsq,1)
      call dfill(NAT3SQ,0.0d00,SCRSQ,1)
      call dfill(NAT3TR,0.0d00,SCRSYM,1)
      call dfill(NAT3,0.0d00,HIEIGS,1)
      call dfill(NAT3,0.0d00,GISR,1)
      call dfill(NAT3SQ,0.0d00,VC,1)
      call dfill(3,0.0d00,UNIVEC,1)
      call dfill(3,0.0d00,SQMAS,1)
C
C     Construct the translation unit vectors;  these are stored in
C     the first three columns of array VC.
C
      IF (MWC) THEN
      VNORM = ONE/SQRT(TOTM)
      DO 30 I=1,3
        DO 20 IAT=1,NATOM
          DO 10 IAXIS=1,3
            IF( IAXIS.NE.I ) THEN
              VC(3*(IAT-1)+IAXIS,I) = ZERO
            ELSE
              VC(3*(IAT-1)+IAXIS,I) = VNORM*SQRT(ATMASS(IAT))
            END IF
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      ELSE
      xnatom = natom
      VNORM = ONE/SQRT(xNATOM)
      DO 35 I=1,3
        DO 25 IATOM=1,NATOM
          DO 15 IAXIS=1,3
            IF( IAXIS.NE.I ) THEN
              VC(3*(IATOM-1)+IAXIS,I) = ZERO
            ELSE
              VC(3*(IATOM-1)+IAXIS,I) = VNORM
            END IF
   15     CONTINUE
   25   CONTINUE
   35 CONTINUE
      ENDIF
*debug:      call stpr_datestp(6,'partit: after tran-vec form')
C
C     Construct rotation unit vectors; store in next 3 columns of VC.
C
      IF (MWC) THEN
      DO 60 IAXIS=1,3
        DO 50 IAT=1,NATOM
           DO 40 J=1,3
             IF (J.EQ.IAXIS) THEN
               SQMAS(J) = ONE*SQRT(ATMASS(IAT))
             ELSE
               SQMAS(J) = ZERO
             END IF
   40      CONTINUE
           CALL stpr_cross(COORD(1,IAT),SQMAS,
     &      VC(3*(IAT-1)+1,IAXIS+3))
   50   CONTINUE
   60 CONTINUE
      ELSE
      DO 65 IAXIS=1,3
        DO 45 J=1,3
          IF (J.EQ.IAXIS) THEN
            UNIVEC(J) = ONE
          ELSE
            UNIVEC(J) = ZERO
          END IF
   45   CONTINUE
        DO 55 IATOM=1,NATOM
          CALL stpr_cross(COORD(1,IATOM),UNIVEC,
     &      VC(3*(IATOM-1)+1,IAXIS+3))
   55   CONTINUE
   65 CONTINUE
      ENDIF
*debug:      call stpr_datestp(6,'partit: after rot-vec form')
C
C     Normalize initial vectors
C
      DO 61 I=1,6
*        RNORM = ZERO
*        DO 62 J=1,NAT3
*          RNORM = RNORM + VC(J,I)**2
*   62   CONTINUE
        Rnorm = ddot(nat3,VC(1,I),1,VC(1,I),1)
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
          call dscal(nat3,rnorm,Vc(1,I),1)
*          DO 63 K=1,NAT3
*            VC(K,I) = VC(K,I)*RNORM
*   63     CONTINUE
        ELSE
          call dfill(nat3,0.0d00,Vc(1,I),1)
*          DO 64 K=1,NAT3
*            VC(K,I) = ZERO
*   64     CONTINUE
        ENDIF
   61 CONTINUE
C      WRITE(6,*)'  Translational & Rotational vectors'
C        DO 66 I=1,6
C          WRITE(6,*)(VC(K,I),K=1,NAT3)
C   66   CONTINUE
C
C     Schmidt orthogonalize the constraint vectors.
C
*debug:      call stpr_datestp(6,'partit: after rot/tran-vec norm')
      NEXTER = 1
      IRED = 0
      DO 130 I=2,6
        DO 90 J=1,I-1
*          DOT = ZERO
*          DO 70 K=1,NAT3
*            DOT = DOT+VC(K,J)*VC(K,I)
*   70     CONTINUE
          Dot = ddot(nat3,VC(1,J),1,VC(1,I),1)
C          WRITE(6,*)'  Overlap between vectors ',I,J,
C     &    ' is ',DOT
*          DO 80 K=1,NAT3
*            VC(K,I) = VC(K,I)-DOT*VC(K,J)
*   80     CONTINUE
          dot=-dot
          call daxpy(nat3,dot,Vc(1,J),1,Vc(1,I),1)
   90   CONTINUE
*        RNORM = ZERO
*        DO 100 K=1,NAT3
*          RNORM = RNORM+VC(K,I)*VC(K,I)
*  100   CONTINUE
         rnorm = ddot(nat3,Vc(1,i),1,Vc(1,i),1)
C        WRITE(6,*)'  Squared norm of vector ',I,' is ',RNORM
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
*          DO 110 K=1,NAT3
*            VC(K,I) = VC(K,I)*RNORM
*  110     CONTINUE
          call dscal(nat3,rnorm,Vc(1,I),1)
          NEXTER = NEXTER+1
        ELSE
          call dfill(nat3,ZERO,Vc(1,I),1)
*          DO 120 K=1,NAT3
*            VC(K,I) = ZERO
*  120     CONTINUE
          ITEMP = I-3
          IRED = I
          WRITE(6,*)'  Dependent rotation vector no.',ITEMP
          WRITE(6,*)'  found in ECKART; linear geometry ',
     &              ' assumed. '
        END IF
  130 CONTINUE


      if (rtdb_get_info(rtdb,'geometry:actlist', ma_type, 
     $     nactive,date)) then
c
c     DO NOT project out rotn/trans if have frozen some cartesians
c
**         write(6,*) ' NOT PROJECTING !!!! '
         nexter = 0
         ired = 0
      endif


      IF(IRED.NE.0.AND.IRED.NE.6)THEN
        call dcopy(nat3,Vc(1,6),1,Vc(1,IRED),1)
*        DO 140 I = 1,NAT3
*          VC(I,IRED) = VC(I,6)
*  140   CONTINUE
      ENDIF
C
C     Find orthogonal vectors which span the internal space.
C
      IINT = 0
      IONE = 0
      DO 210 I=1,NAT3
        IINT = IINT + 1
        IONE = IONE + 1
        call dfill(nat3,zero,Vc(1,(nexter+iint)),1)
*        DO 150 J = 1,NAT3
*          VC(J,NEXTER+IINT) = 0.0D0
*  150   CONTINUE
        VC(IONE,NEXTER+IINT) = ONE
        DO 180 J=1,NEXTER+IINT-1
*          DOT = ZERO
*          DO 160 K=1,NAT3
*            DOT = DOT+VC(K,J)*VC(K,NEXTER+IINT)
*  160     CONTINUE
          dot = ddot(nat3,Vc(1,j),1,Vc(1,(NEXTER+IINT)),1)
          dot = -dot
          call daxpy(nat3,dot,VC(1,J),1,VC(1,(NEXTER+IINT)),1)
*          DO 170 K=1,NAT3
*            VC(K,NEXTER+IINT) = VC(K,NEXTER+IINT)-DOT*VC(K,J)
*  170     CONTINUE
  180   CONTINUE
*        RNORM = ZERO
*        DO 190 K=1,NAT3
*          RNORM = RNORM+VC(K,NEXTER+IINT)*VC(K,NEXTER+IINT)
*  190   CONTINUE
         rnorm = ddot(nat3,VC(1,(NEXTER+IINT)),1,VC(1,(NEXTER+IINT)),1)
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
          call dscal(nat3,rnorm,VC(1,(NEXTER+IINT)),1)
*          DO 200 K=1,NAT3
*            VC(K,NEXTER+IINT) = VC(K,NEXTER+IINT)*RNORM
*  200     CONTINUE
          IF(NEXTER+IINT.EQ.NAT3)GOTO 220
        ELSE
          IINT = IINT - 1
        END IF
  210 CONTINUE
  220 CONTINUE
*debug:      call stpr_datestp(6,'partit: after ortho/norm vecs ')
      NINTER = NAT3 - NEXTER
C
C     Form the internal hessian matrix first in a scratch array
C     and than move to HESS.
C
*debug:      call stpr_datestp(6,'partit: b4 internal hess form')
#define PHASE2
#if defined(OLDIFORM)      
      DO 260 M=1,NINTER
        MV = NEXTER + M
        DO 250 N=1,M
          NV = NEXTER + N
          SCRSYM(ISYM2(M,N)) = ZERO
          DO 240 I=1,NAT3
            DO 230 J=1,NAT3
              SCRSYM(ISYM2(M,N)) = SCRSYM(ISYM2(M,N)) +
     &                       VC(I,MV)*VC(J,NV)*HESS(ISYM2(I,J))
  230       CONTINUE
  240     CONTINUE
  250   CONTINUE
  260 CONTINUE
#endif
#if defined(PHASE1)      
      is2_mn = 0
      DO 260 M=1,NINTER
        MV = NEXTER + M
        DO 250 N=1,M
          NV = NEXTER + N
          is2_mn = is2_mn + 1
          is2_ij = 0
          DO 240 I=1,NAT3
            DO 230 J=1,I
              is2_ij = is2_ij +1 
              if (i.eq.j) then
                SCRSYM(is2_mn) = SCRSYM(is2_mn) +
     &                 VC(I,MV)*VC(J,NV)*HESS(is2_ij)
              else
                SCRSYM(is2_mn) = SCRSYM(is2_mn) +
     &                 (VC(I,MV)*VC(J,NV) + VC(J,MV)*VC(I,NV))*
     &                 HESS(is2_ij)
              endif
230         CONTINUE
240       CONTINUE
250     CONTINUE
260   CONTINUE
#endif
#if defined(PHASE2)
      is2_mn = 0
      DO 260 M=1,NINTER
        MV = NEXTER + M
        DO 250 N=1,M
          NV = NEXTER + N
          is2_mn = is2_mn + 1
          DO 240 I=1,NAT3
            is2_i = i*(i-1)/2
            DO 230 J=1,(I-1)
              is2_ij = is2_i + j
              SCRSYM(is2_mn) = SCRSYM(is2_mn) +
     &               (VC(I,MV)*VC(J,NV) + VC(J,MV)*VC(I,NV))*
     &               HESS(is2_ij)
230         CONTINUE
240       CONTINUE
          do 00241 i=1,nat3
            is2_ij = i*(i+1)/2
            SCRSYM(is2_mn) = SCRSYM(is2_mn) +
     &             VC(I,MV)*VC(i,NV)*HESS(is2_ij)
00241     continue
250     CONTINUE
260   CONTINUE
#endif
*debug:      call stpr_datestp(6,'partit: after internal hess form')
      call dfill(NAT3TR,0.0d00,HESS,1)
      DO 270 I = 1,NAT3TR
C
C       Zap any elements with absolute value less than 1.D-12
C
        IF(DABS(SCRSYM(I)).GT.1.D-12) HESS(I) = SCRSYM(I)
  270 CONTINUE
*debug:      call stpr_datestp(6,'partit: after internal hess zap')
*debug:      write(6,*)'hess stpr_partit after zap', nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
      call dfill(NAT3TR,0.0d00,SCRSYM,1)
C
C     Form the internal gradient first in a scratch array and then
C     move to GRAD.
C
      DO 290 M=1,NINTER
        MV = NEXTER + M
        SCRSYM(M) = ZERO
        DO 280 I=1,NAT3
          SCRSYM(M) = SCRSYM(M) + VC(I,MV)*GRAD(I)
  280   CONTINUE
  290 CONTINUE
      call dfill(NAT3,0.0d00,GRAD,1)
*debug:      call stpr_datestp(6,'partit: after internal grad form')
      DO 300 I = 1,NINTER
C
C       Zap any elements with absolute value less than 1.D-12
C
        IF(DABS(SCRSYM(I)).GT.1.D-12) GRAD(I) = SCRSYM(I)
  300 CONTINUE
      call dfill(NAT3TR,0.0d00,SCRSYM,1)
      call dfill(NAT3SQ,0.0d00,SCRSQ,1)
C
C     Construct a reduced internal hessian and gradient to remove
C     any eigenmodes with zero force.
C
*debug:      write(6,*)'hess stpr_partit before jacobi',nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
*debug:      write(6,*)' '
*debug:      write(6,*)' hessian full print'
*debug:      call stpr_prntpd(hess,nat3tr,nat3,6)
*debug:      write(6,*)' '
*debug:      write(6,*)' hessian ninternal print'
*debug:      call stpr_prntpd(hess,(ninter*(ninter+1)/2),ninter,6)
*debug:      call stpr_datestp(6,'partit: before jacobi')
      CALL stpr_sjacobi(NAT3,NINTER,HESS,HIEIGS,SCRSQ,IERR,
     .     .true.)
*debug:      call stpr_datestp(6,'partit: after  jacobi')
*debug:      write(6,*)'hess stpr_partit after jacobi',nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hess(i)).gt.1.0d-5)
*debug:     &        write(6,*)' hess (',i,') =',hess(i)
*debug:      enddo
      IF(LPRINT)THEN
        WRITE(6,*)'  The eigenvalues of the internal hessian. '
        CALL stpr_matout(HIEIGS,NINTER,1,NINTER)
c        WRITE(6,*)(HIEIGS(I),I=1,NINTER)
      ENDIF
      IF (MWC) THEN
        WRITE(6,*)'  The frequencies in cm-1 '
        DO 306 I=1,NINTER
          IF (HIEIGS(I).GT.0.D0) THEN
            FREQ=SQRT(HIEIGS(I))*AUCM
            WRITE(6,*) FREQ
          ELSE
            FREQ=SQRT(-HIEIGS(I))*AUCM
            WRITE(6,*)'  i',FREQ
          ENDIF
  306   CONTINUE
      ENDIF
C
C     Form the internal gradient in the spectral representation.
C
      GLNGTH = 0.0D0
      GMXMUM = 0.0D0
      DO 320 I = 1,NINTER
        GLNGTH = GLNGTH + GRAD(I)**2
        IF(ABS(GRAD(I)).GT.GMXMUM) GMXMUM=ABS(GRAD(I))
        GISR(I) = 0.0D0
        DO 310 J = 1,NINTER
          GISR(I) = GISR(I) + SCRSQ(J,I)*GRAD(J)
  310   CONTINUE
  320 CONTINUE
      GLNGTH = SQRT(GLNGTH)
      IF(LPRINT)THEN
        WRITE(6,*)'  The corresponding internal gradient '
        WRITE(6,*)'  in the spectral representaion. '
        CALL stpr_matout(GISR,NINTER,1,NINTER)
c        WRITE(6,*)(GISR(I),I=1,NINTER)
      ENDIF
C
C    Determine all eigenvectors of internal hessian.
C
      call dfill(NAT3SQ,0.0d00,VECSQ,1)
      DO 330 N = 1, NINTER
        DO 340 ICART = 1,NAT3
          VECSQ(ICART,N) = 0.0D0
          DO 350 J = 1,NINTER
            VECSQ(ICART,N) = VECSQ(ICART,N) +
     &                           VC(ICART,NEXTER+J)*SCRSQ(J,N)
  350     CONTINUE
  340   CONTINUE
  330 CONTINUE
C
C     Write out the eigenvectors.
C
      IF(LPRINT)THEN
        WRITE(6,1001)
        WRITE(6,*)'  The corresponding eigenvectors '
        CALL stpr_matout(VECSQ,NAT3,NINTER,NAT3)
c        WRITE(6,1002)
c        WRITE(6,1001)
c        DO 355 J = 1,NINTER
c          WRITE(6,1000)J, (VECSQ(M,J),M=1,NAT3)
c  355   CONTINUE
      ENDIF
C
C     Tricky stuff.
C
      IF (COOK) THEN
        COOK = .FALSE.
        MWC = .TRUE.
        call dfill(NAT3TR,0.0d00,HESS,1)
        call dfill(NAT3,0.0d00,GRAD,1)
        call dfill(NAT3,0.0d00,COORD,1)
        IERR = 0
        CALL stpr_recall(HESS,GRAD,COORD,ENERGY,IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'  Unsuccessful "RECALL" attempt. '
          STOP
        ENDIF
C
C       Calculate vector of the center of mass and the inertia tensor.
C
        CALL stpr_mgeom(COORD,ATMASS,CMASS,TENIN)
C
C       Mass weight hessian.
C
        CALL stpr_hmass(HESS,ATMASS)
C
C       Mass weight gradient.
C
        CALL stpr_gmass(GRAD,ATMASS)
        GOTO 999
      ENDIF
C      ENDIF
C
C     Checkup for pathological situations in minimization.
C
      IF (OMIN) THEN
        DO 361 I = 1,NINTER
          IF ((HIEIGS(I).LT.0.D0).AND.(DABS(GISR(I)).LE.SYMCHK))
     1    THEN
            WRITE(6,1001)
            WRITE(6,*)'  W A R N I N G '
            WRITE(6,*)'  In mode ',I,' you are on a ridge. '
          ENDIF
 361    CONTINUE
      ENDIF
C
C     Find NMODE
C     The following 2 lines are to satisfy compiler warnings.
C
      NMIN = 0
      NPLUS = 2
      IF (TS.OR.STD) THEN
        ICOUNT = 0
        DO 362 I = 1,NINTER
          IF (DABS(GISR(I)).GT.SYMCHK) THEN
            ICOUNT = ICOUNT + 1
          ENDIF
          IF (ICOUNT.EQ.NMODE) THEN
            NMODE = I
            GOTO 363
          ENDIF
 362    CONTINUE
 363    CONTINUE
        IF(TRACK) THEN
          CALL stpr_fndmde(STEP,SELECT,VECSQ,NINTER,NAT3,NMODE)
          call dfill(NAT3,0.0d00,STEP,1)
          TRACK = .FALSE.
        ENDIF
        NMIN = NMODE - 1
        NPLUS = NMODE + 1
      ENDIF
C
C     Checkup for pathological situations in TS or STD walk.
C
C     Checkup for negative curvatures and ridges along
C     modes transverse to NMODE.
C
      IF (TS.OR.STD) THEN
        DO 364 I = 1,NMIN
          IF (HIEIGS(I).LT.0.D0) THEN
            WRITE(6,1001)
            WRITE(6,*)'  W A R N I N G '
            WRITE(6,*)'  Negative hessian eigenvalue for'
            WRITE(6,*)'  tranverse mode # ',I
            IF (DABS(GISR(I)).LE.SYMCHK) THEN
              WRITE(6,*)'  In mode ',I,' you are on a ridge. '
            ENDIF
          ENDIF
  364   CONTINUE
        DO 365 I = NPLUS,NINTER
          IF (HIEIGS(I).LT.0.D0) THEN
            WRITE(6,1001)
            WRITE(6,*)'  W A R N I N G '
            WRITE(6,*)'  Negative hessian eigenvalue for'
            WRITE(6,*)'  tranverse mode # ',I
            IF (DABS(GISR(I)).LE.SYMCHK) THEN
              WRITE(6,*)'  In mode ',I,' you are on a ridge. '
            ENDIF
          ENDIF
  365   CONTINUE
      ENDIF
C
C     Checkup for pathological situations in TS walk.
C
      IF (TS) THEN
C
C       Consider the case of zero force along
C       the reaction coordinate.
C
        IF (DABS(GISR(NMODE)).LE.SYMCHK) THEN
          IF (HIEIGS(NMODE).GT.0.D0) THEN
            WRITE(6,*)'  You need a nonzero force along'
            WRITE(6,*)'  mode # ',NMODE,' to go uphill. '
            call errquit('stepper:partit fatal error',911, UNKNOWN_ERR)
          ELSE
            TS = .FALSE.
            NMODE = 0
            OMIN = .TRUE.
            WRITE(6,*)'  There is no force along the reaction'
            WRITE(6,*)'  coordinate mode. '
            WRITE(6,*)'  Minimization along transverse modes'
            WRITE(6,*)'  will be performed. '
          ENDIF
        ENDIF
      ENDIF
C
C     Checkup for pathological situations in STD walk.
C
      IF (STD) THEN
C
C       Consider the case of zero force along
C       the reaction coordinate.
C
        IF (DABS(GISR(NMODE)).LE.SYMCHK) THEN
          IF (HIEIGS(NMODE).LT.0.D0) THEN
            WRITE(6,*)'  You need a nonzero force along'
            WRITE(6,*)'  mode # ',NMODE,' to go downhill. '
            call errquit('stepper:partit fatal error',911, UNKNOWN_ERR)
          ELSE
            STD = .FALSE.
            NMODE = 0
            OMIN = .TRUE.
            WRITE(6,*)'  There is no force along the reaction'
            WRITE(6,*)'  coordinate mode. '
            WRITE(6,*)'  Minimization along transverse modes'
            WRITE(6,*)'  will be performed. '
          ENDIF
        ENDIF
      ENDIF
C
C     Determine which modes are real and throw the rest out.
C     *** NOTE *** An eigenvalue with a corresponding force
C     less than SYMCHK is NOT counted as a real eigenvalue.
C
C     Determine number of modes with no forces
C     and softer than NMODE.
C
      NOFORC = 0
c      WRITE(6,*)'  Removing modes with zero forces. '
      call dfill(NAT3,0.0d00,GRAD,1)
      call dfill(NAT3TR,0.0d00,HESS,1)
      IREAL = 0
      IUNRL = 0
      DO 370 N = 1, NINTER
        IF(DABS(GISR(N)).GT.SYMCHK)THEN
          IREAL = IREAL + 1
          HESS(ISYM2(IREAL,IREAL)) = HIEIGS(N)
          GRAD(IREAL) = GISR(N)
          DO 380 ICART = 1,NAT3
            VC(ICART,IREAL) = VECSQ(ICART,N)
  380     CONTINUE
        ELSE
          IUNRL = IUNRL + 1
          IF (N.LT.NMODE) NOFORC = NOFORC + 1
        ENDIF
  370 CONTINUE
      if (oprint) then
         WRITE(6,*)'  Number of modes removed with force (in the',
     &        ' spectral representation)'
         WRITE(6,'(a,f10.6,a,i4)')'   less than: ',SYMCHK,' is: ',
     &        IUNRL
      endif
      NINTER = IREAL
      DO 400 N = 1,NINTER
        HIEIGS(N) = HESS(ISYM2(N,N))
        GISR(N) = GRAD(N)
  400 CONTINUE
      IF(NINTER.EQ.0)THEN
         if (oprint) then
            WRITE(6,*)'  There are no internal hessian eigenmodes '
            WRITE(6,*)'  remaining with corresponding forces '
            WRITE(6,'(a,f10.6)') '  greater than: ',SYMCHK
            WRITE(6,*)'  Convergence is set .TRUE. '
         endif
         RETURN
      ELSE
         if (oprint) then
            WRITE(6,*)'  PARTIT has partitioned the input hessian '
            WRITE(6,*)'  and gradient into two subspaces;  An '
            WRITE(6,*)'  "external" hessian and gradient composed '
            WRITE(6,*)'  of ',NEXTER,' independent translations '
            WRITE(6,*)'  and rotations, and an "internal" hessian '
            WRITE(6,*)'  and gradient with a dimension of ',NINTER
         endif
      ENDIF
      IF (NMODE.GT.0) THEN
C
C     Determine current NMODE.
C
        NMODE = NMODE - NOFORC
C
C       Calculate the number of negative hessian eigenvalues
C       and check whether transverse modes have positive hessian
C       eigenvalues.
C
        DO 410 N = 1, NINTER
          IF(HIEIGS(N).LT.0.D0)THEN
            IF (N.NE.NMODE) TROUB = .TRUE.
          ENDIF
  410   CONTINUE
C
C       Make RPATH active if a transverse mode has
C       a negative hessian eigenvalue.
C
        IF ((TROUB).AND.(ICALC.EQ.2)) THEN
          RPATH = .TRUE.
        ENDIF
      ENDIF
C
C     If RPATH is true and the current geometrical point is
C     not in the bottom of the valley then
C     minimize energy with respect to all modes which
C     are orthogonal to the reaction coordinate.
C
      IF ((RPATH).AND.(NMODE.GT.0)) THEN
        NMIN = NMODE - 1
        NPLUS = NMODE + 1
        IREAL = 0
        DO 420 N = 1, NMIN
          IF(DABS(GISR(N)).GT.RCONVG) IREAL = IREAL + 1
  420   CONTINUE
        DO 430 N = NPLUS, NINTER
          IF(DABS(GISR(N)).GT.RCONVG) IREAL = IREAL + 1
  430   CONTINUE
        IF ((IREAL.GT.0).OR.(TROUB)) THEN
          WRITE(6,*)'  RPATH is .TRUE.'
          IF (IREAL.GT.0) THEN
            WRITE(6,*)'  There are ',IREAL,' modes'
            WRITE(6,*)'  different than NMODE = ',NMODE
            WRITE(6,*)'  with forces bigger than ',RCONVG,'.'
          ENDIF
          IF (TROUB) THEN
            WRITE(6,*)'  There are transverse modes with '
            WRITE(6,*)'  negative hessian eigenvalues and '
            WRITE(6,*)'  nonzero forces. '
          ENDIF
          WRITE(6,*)'  A downhill step to the streambed'
          WRITE(6,*)'  will be performed. '
          call dfill(NAT3,0.0d00,GRAD,1)
          call dfill(NAT3TR,0.0d00,HESS,1)
          IREAL = 0
          DO 440 N = 1, NMIN
            IF(DABS(GISR(N)).GT.SYMCHK)THEN
              IREAL = IREAL + 1
              HIEIGS(IREAL) = HIEIGS(N)
              HESS(ISYM2(IREAL,IREAL)) = HIEIGS(IREAL)
              GISR(IREAL) = GISR(N)
              GRAD(IREAL) = GISR(IREAL)
              DO 450 M=1, NAT3
                VC(M,IREAL) = VC(M,N)
  450         CONTINUE
            ENDIF
  440     CONTINUE
          DO 460 N = NPLUS, NINTER
            IF(DABS(GISR(N)).GT.SYMCHK)THEN
              IREAL = IREAL + 1
              HIEIGS(IREAL) = HIEIGS(N)
              HESS(ISYM2(IREAL,IREAL)) = HIEIGS(IREAL)
              GISR(IREAL) = GISR(N)
              GRAD(IREAL) = GISR(IREAL)
              DO 470 M=1, NAT3
                VC(M,IREAL) = VC(M,N)
  470         CONTINUE
            ENDIF
  460     CONTINUE
          NINTER = IREAL
          TS = .FALSE.
          STD = .FALSE.
          OMIN = .TRUE.
          NMODE = 0
        ELSE
          WRITE(6,*)'  You are in the streambed. '
        ENDIF
      ENDIF
*debug:      call stpr_datestp(6,'partit: end')
      RETURN
 1000 FORMAT('  Eigenvector : ',I2,/,
     &       1X,(3(2X,F12.7)))
 1001 FORMAT(/)
 1002 FORMAT('  ',9('-'),' Hessian eigenvectors (3N-(6 ',
     &          ' or 5)) ',9('-'))
      END
C> @}
