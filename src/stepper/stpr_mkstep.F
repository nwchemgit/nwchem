C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_mkstep(
     &       HIEIGS,GISR,VC,STEP,SELECT,SCRVEC,ATMASS,
     &       ENERGY, CONVGE, CONVGG, CONVGGM, GRAD_MAX,
     &       E2NEW,NINTER,NAT3,NMODE,TRUSTI,BCKSTP,
     &       CNVGRD,rtdb)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      LOGICAL BCKSTP,CNVGRD,LPRINT, oprint, oprint_step
#include "chow.fh"
#include "util.fh"
      integer rtdb
      COMMON / THRESH / RCONVG,RSTD
      COMMON / CFACE / IWCTR,NATOM,ICALC
      DIMENSION HIEIGS(NINTER),GISR(NINTER),VC(NAT3,NINTER),
     &          STEP(NAT3),SELECT(NAT3),SCRVEC(NAT3),
     &          ATMASS(NATOM)
C
C     Convergence/symmetry tolerances.
C
#include "stpr_clean.fh"
      IATOM(I)  = (I+2)/3
c      LPRINT = .TRUE.
      LPRINT = util_print('makestep',print_debug)
      oprint = util_print('information', print_default)
      oprint_step = util_print('step', print_low)
      IF (MWC) THEN
        CONVGG=1.0D-08
      ENDIF
C
C     Write Hessian eigenvalues.
C
      IF(LPRINT)THEN
        WRITE(6,*)'  Hessian eigenvalues: '
        CALL stpr_matout(HIEIGS,NINTER,1,NINTER)
c        WRITE(6,*)(HIEIGS(I),I=1,NINTER)
      ENDIF
C
C     Write gradient components.
C
      IF(LPRINT)THEN
        WRITE(6,*)'  Gradient components in the spectral ',
     &            'representation: '
        CALL stpr_matout(GISR,NINTER,1,NINTER)
c        WRITE(6,*)(GISR(I),I=1,NINTER)
      ENDIF
C
C     Calculate gradient norm.
C
      GLNGTH = 0.0D0
      DO 10 I = 1,NINTER
        GLNGTH = GLNGTH + GISR(I)**2
   10 CONTINUE
      GLNGTH = SQRT(GLNGTH)
C
C     gradient maximum
C
      GMXMUM = GRAD_MAX
C
C     Count the number of negative eigenvalues.
C
      NNEG = 0
      DO 20 N = 1, NINTER
        IF(HIEIGS(N).LT.0.D0)THEN
          NNEG = NNEG + 1
        ENDIF
   20 CONTINUE
C
C     Check whether convergence has been reached.
C     Also set any TRUST info.
C
      CALL stpr_cheker(ENERGY,EOLD,E2OLD,EDIFF,
     &                 GMXMUM,GLNGTH,ERATIO,
     &                 CONVGE, CONVGG, CONVGGM,
     &                 TRUSTI,TRUST,TRUSTO,SOLD,NNEG,NWALK,
     &                 BCKSTP,CNVGRD)
C
C     Make the steepest descent step if appropriate
C
   25 CONTINUE
      IF (STD) THEN
         if (oprint)
     $        WRITE(6,*)'  Calculating the steepest descent step. '
         DO 30 I = 1,NINTER
           STEP(I) = -GISR(I)
   30    CONTINUE
C
C        Find the length of the steepest descent STEP.
C
         SLNGTH = 0.0D0
         DO 35 I = 1,NINTER
           SLNGTH = SLNGTH + STEP(I)**2
   35    CONTINUE
         SLNGTH = SQRT(SLNGTH)
C
C        Check and scale step if needed
C
C
C        Step too long ...
C
         IF (SLNGTH.GT.RSTD) THEN
           DO 7117 I=1,NINTER
             STEP(I) = STEP(I)*RSTD/SLNGTH
 7117      CONTINUE
         ELSE
C
C        Step shorter than RSTD and TRTR is true
C
           IF (TRTR) THEN
             DO 7127 I=1,NINTER
               STEP(I) = STEP(I)*RSTD/SLNGTH
 7127        CONTINUE
             GOTO 75
           ENDIF
C
C        Step shorter than RSTD, TROT is true, NMODE's
C        hessian eigenvalue still negative
C
           IF(TROT.AND.(HIEIGS(NMODE).LT.0.D0)) THEN
             DO 7137 I=1,NINTER
               STEP(I) = STEP(I)*RSTD/SLNGTH
 7137        CONTINUE
             GOTO 75
           ENDIF
         ENDIF
         GOTO 75
      ENDIF
C
C     TS or OMIN walks
C
C     Make the N-R step if:
C     - all hessian eigenvalues positive and minimization
C     - one hessian eigenvalue negative and walk to a transition state
C
      IF((NNEG.EQ.0.AND.OMIN).OR.
     &(NNEG.EQ.1.AND.TS.AND.(HIEIGS(1).LT.0.D0))) THEN
C
C       Determine the N-R step.
C
         if (oprint)
     $        WRITE(6,*)'  Calculating the Newton-Raphson step. '
        DO 40 I = 1,NINTER
          IF(ABS(HIEIGS(I)).GT.1.D-14)THEN
            STEP(I) = -GISR(I)/HIEIGS(I)
          ELSE
            IF (I.EQ.NMODE) THEN
              STEP(I) = DSIGN(TRUST,GISR(I))
            ELSE
              STEP(I) = DSIGN(TRUST,-GISR(I))
            ENDIF
          ENDIF
          IF ((ABS(STEP(I)).GT.TRUST).AND.(ICALC).EQ.2) THEN
            STEP(I) = DSIGN(TRUST,STEP(I))
          ENDIF
   40   CONTINUE
C
C     Minimization with some negative hessian eigenvalues
C
      ELSEIF (OMIN) THEN
        DO 50 I = 1,NINTER
          IF(ABS(HIEIGS(I)).GT.1.D-14)THEN
            IF (HIEIGS(I).GT.0.D0) THEN
              STEP(I) = -GISR(I)/HIEIGS(I)
            ELSE
              STEP(I) = DSIGN(TRUST,-GISR(I))
            ENDIF
          ELSE
            STEP(I) = DSIGN(TRUST,-GISR(I))
          ENDIF
          IF ((ABS(STEP(I)).GT.TRUST).AND.(ICALC).EQ.2) THEN
            STEP(I) = DSIGN(TRUST,STEP(I))
          ENDIF
   50   CONTINUE
C
C     Uphill along NMODE with all hessian eigenvalues positive
C
      ELSEIF (TS.AND.(NNEG.EQ.0)) THEN
        DO 60 I = 1,NINTER
          IF(ABS(HIEIGS(I)).GT.1.D-14)THEN
            IF (I.EQ.NMODE) THEN
              STEP(I) = DSIGN(TRUST,GISR(I))
            ELSE
              STEP(I) = -GISR(I)/HIEIGS(I)
            ENDIF
          ELSE
            IF (I.EQ.NMODE) THEN
              STEP(I) = DSIGN(TRUST,GISR(I))
            ELSE
              STEP(I) = DSIGN(TRUST,-GISR(I))
            ENDIF
          ENDIF
          IF ((ABS(STEP(I)).GT.TRUST).AND.(ICALC).EQ.2) THEN
            STEP(I) = DSIGN(TRUST,STEP(I))
          ENDIF
   60   CONTINUE
C
C     Uphill with fishy hessian structure
C
      ELSEIF (TS.AND.(NNEG.GT.0)) THEN
        DO 70 I = 1,NINTER
          IF(ABS(HIEIGS(I)).GT.1.D-14)THEN
            IF (I.EQ.NMODE) THEN
              IF (HIEIGS(I).LT.0.D0) THEN
                STEP(I) = -GISR(I)/HIEIGS(I)
              ELSE
                STEP(I) = DSIGN(TRUST,GISR(I))
              ENDIF
            ELSE
              IF (HIEIGS(I).GT.0.D0) THEN
                STEP(I) = -GISR(I)/HIEIGS(I)
              ELSE
                STEP(I) = DSIGN(TRUST,-GISR(I))
              ENDIF
            ENDIF
          ELSE
            IF (I.EQ.NMODE) THEN
              STEP(I) = DSIGN(TRUST,GISR(I))
            ELSE
              STEP(I) = DSIGN(TRUST,-GISR(I))
            ENDIF
          ENDIF
          IF ((ABS(STEP(I)).GT.TRUST).AND.(ICALC).EQ.2) THEN
            STEP(I) = DSIGN(TRUST,STEP(I))
          ENDIF
   70   CONTINUE
      ELSE
        WRITE(6,*) 'Unexpected type of walk. Program aborts'
        call errquit('stepper:mkstep: fatal error',911, CALC_ERR)
      ENDIF
   75 CONTINUE
      IF(LPRINT)THEN
         WRITE(6,*)'  The STEP in the spectral representation'
         CALL stpr_matout(STEP,NINTER,1,NINTER)
c     WRITE(6,*)(STEP(I),I=1,NINTER)
      ENDIF
C
C     Temporarily store this STEP in SCRVEC for later printing.
C
      DO 76 I = 1,NINTER
        SCRVEC(I) = STEP(I)
   76 CONTINUE
C
C     Find the length of STEP.
C
      SLNGTH = 0.0D0
      DO 80 I = 1,NINTER
        SLNGTH = SLNGTH + STEP(I)**2
   80 CONTINUE
      SLNGTH = SQRT(SLNGTH)
C
C     Calculation with first derivatives only
C     Scale down the step if needed
C
      IF ((SLNGTH.GT.TRUST).AND.(ICALC.EQ.1)) THEN
        DO 85 I=1,NINTER
           STEP(I) = STEP(I)*TRUST/SLNGTH
 85     CONTINUE
        IF(LPRINT)THEN
          WRITE(6,*)'  After scaling down ....',
     $          ' the STEP in the spectral representation'
          CALL stpr_matout(STEP,NINTER,1,NINTER)
c          WRITE(6,*)(STEP(I),I=1,NINTER)
        ENDIF
C
C       Find the length of STEP.
C
        SLNGTH = 0.0D0
        DO 90 I = 1,NINTER
          SLNGTH = SLNGTH + STEP(I)**2
 90     CONTINUE
        SLNGTH = SQRT(SLNGTH)
*        THRCLN=1.0D-04   ! RJH ... use data set 1e-6 to keep symm better
      ENDIF
      if (oprint)
     $     WRITE(6,'(a,f10.6)')
     $     '  The length (norm) of the STEP is ',SLNGTH
C
C     Determine TRUSTN for STD and first derivative walks
C
      IF (STD) THEN
        TRUSTN = MIN(RSTD,SLNGTH)
      ENDIF
      IF ((TS.OR.OMIN).AND.(ICALC.EQ.1)) THEN
        TRUSTN = MIN(TRUST,SLNGTH)
      ENDIF
C
C     Print small amount of useful information.
C
      if (oprint) then
         WRITE(6,1000)
         WRITE(6,1001)
         IEND = MIN(10,NINTER)
         DO 100 I = 1, IEND
            WRITE(6,1002)I,HIEIGS(I),GISR(I),SCRVEC(I),STEP(I)
 100     CONTINUE
 1000    FORMAT(/,'  The lowest few components in the bases ',/,
     &        '  generated by the spectral decomposition ',/,
     &        '  of the internal hessian (non-zero force). ')
 1001    FORMAT(/,'  Component  Hessian     Grad      Step',
     &        '   Scaled-Step',
     &        /,'  ---------  --------  --------  --------',
     &        ' -----------')
 1002    FORMAT('  ',4X,I2,4X,4(1X,D8.2,1X))
         WRITE(6,*)
      endif
C
C     Back transform STEP from the spectral decomp. to cart.
C     Store temporarily in SCRVEC.
C
      DO 120 I = 1,NAT3
        SCRVEC(I) = 0.0D0
        DO 110 J = 1,NINTER
          SCRVEC(I) = SCRVEC(I) + VC(I,J)*STEP(J)
  110   CONTINUE
        IF (MWC) THEN
          SCRVEC(I)=SCRVEC(I)/SQRT(ATMASS(IATOM(I)))
        ENDIF
  120 CONTINUE
C
C     End of stepping algorithm.  Wrap it up and return.
C
C     Clean up the STEP ... truncate any noise in the step
C     less then THRCLN.  Clean up the step in both the
C     cartesian basis and in the spectral representation.
C
      if (oprint)
     $     WRITE(6,'(a,f10.6)')
     $     '  Clean-up the step; truncate at: ',THRCLN
      CALL stpr_clnstp(STEP,THRCLN)
      CALL stpr_clnstp(SCRVEC,THRCLN)
C
C     Was this step meaningful?
C
      SLNGTH = 0.0D0
      DO 130 I = 1,NAT3
        SLNGTH = SLNGTH + SCRVEC(I)**2
 130  CONTINUE
      SLNGTH = SQRT(SLNGTH)
      if (oprint_step) then
         WRITE(6,'(a,f10.6)')
     $        '  In cartesian coordinates; norm of this step is: ',
     $        SLNGTH
      endif
c
c     RJH .. this was causing premature convergence.
c
c$$$      IF (SLNGTH.LE.CONVGG) THEN
c$$$         if (oprint_step) then
c$$$            WRITE(6,'(a,f10.6)')
c$$$     $           '  Norm of this step is less than ',CONVGG
c$$$            WRITE(6,*)'  This geometry walk is HISTORY. '
c$$$         endif
c$$$         CNVGRD = .TRUE.
c$$$        RETURN
c$$$      ENDIF
C
C     Find the largest component of STEP.
C
      SLNGTH = 0.0D0
      DO 140 I = 1,NINTER
        SLNGTH = MAX(ABS(SLNGTH),ABS(STEP(I)))
  140 CONTINUE
      IF ((TS.OR.OMIN).AND.(ICALC.EQ.2)) THEN
        TRUSTN = MIN(TRUST,SLNGTH)
      ENDIF
C
C     Get the predicted 2nd order energy.
C
C      ... Change STEP in the scaled coordinate back
C          to Y to calculate the energy expression (2nd order).
C
      FTERM = 0.0D0
      QTERM = 0.0D0
      DO 150 I=1,NINTER
        FTERM = FTERM + STEP(I)*GISR(I)
        QTERM = QTERM + 0.5D0*HIEIGS(I)*STEP(I)**2
  150 CONTINUE
      E2NEW = FTERM + QTERM
C
C     In the case of the STD walk check whether the predicted
C     energy change is negative
C
      IF (STD.AND.(E2NEW.GT.0.D0)) THEN
         if (oprint_step) then
            WRITE(6,*)'  For RSTD = ',RSTD,' the predicted energy '
            WRITE(6,*)'  change was positive',E2NEW
            RSTD=RSTD/2.D0
            WRITE(6,*)'  RSTD has been decreased to ',RSTD
         endif
        GOTO 25
      ENDIF
C
C     Move the step vector in the cartesian basis currently
C     stored in SCRVEC) into vector STEP replacing the step vector
C     in the spectral representation.
C
      call dfill(NAT3,0.0d00,STEP,1)
      DO 160 I = 1,NAT3
        STEP(I) = SCRVEC(I)
  160 CONTINUE
c
      call stpr_active(rtdb,step,nat3)
C
C     Put step parameters to save common and write them to tape 41.
C
      CALL stpr_putp(GMXMUM,GLNGTH,SLNGTH,EOLD,ENERGY,
     &           EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &           E2NEW,NWALK,NNEG)
      RETURN
      END
      subroutine stpr_active(rtdb,step,nat3)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer rtdb    ! [input] rtdb handle
      integer nat3    ! [input] 3*number of atoms (dimen. of step)
      double precision step(nat3) ! step vector
c
      logical oprint  ! print info here
      integer nat       ! number of atoms
      integer nactive   ! number of active atoms
      integer h_actlist ! ma handle for actlist
      integer k_actlist ! ma index for actlist
      integer h_step    ! ma handle for copy of step
      integer k_step    ! ma index for copy of step
      integer i         ! loop index
      integer atom      ! atom index
      integer indx      ! offset in step and h_step
      integer ma_type
c
* check for active atom list
*..
      if (.not.rtdb_ma_get(rtdb,'geometry:actlist', ma_type, 
     $     nactive,h_actlist)) return
c
      oprint = util_print('step',print_debug)
      if (oprint) then
        write(luout,*)' step before removal of active atoms'
        call output(step,1,nat3,1,1,nat3,1,1)
      endif
c
      nat = nat3/3
      if (nactive.le.0 .or. nactive.gt.nat) call errquit
     &    ('stpr_active: invalid number of active atoms',nactive,
     &       INPUT_ERR)
      if (.not.MA_get_index(h_actlist,k_actlist)) call errquit
     &    ('stpr_active: ma_get_index failed for actlist',911, MA_ERR)
*
* .. get space for step
*
      if (.not.MA_push_get(mt_dbl,nat3,'stpr:alloc:active:step',
     &    h_step,k_step)) call errquit
     &    ('stpr_active: ma_push_get failed for step',911, MA_ERR)
*
* .. copy step to scratch space
*
      call dcopy(nat3,step,1,dbl_mb(k_step),1)
      call dfill(nat3,0.0d00,step,1)
*
*.. modify step to have active atoms with steps
*
      do i = 0,(nactive-1)
        atom = int_mb(k_actlist+i)
        indx = 3*(atom-1) + 1
        step(indx)   = dbl_mb(k_step+indx-1)
        step(indx+1) = dbl_mb(k_step+indx)
        step(indx+2) = dbl_mb(k_step+indx+1)
      enddo
      if (.not.MA_free_heap(h_actlist)) call errquit
     &    ('stpr_active: ma_free_heap failed for actlist',911, MA_ERR)
      if (.not.MA_pop_stack(h_step)) call errquit
     &    ('stpr_active: ma_pop_stack failed for h_step',911, MA_ERR)
c
      if (oprint) then
        write(luout,*)' '
        write(luout,*)' step after removal of active atoms'
        call output(step,1,nat3,1,1,nat3,1,1)
      endif
c
      end
C> @}
