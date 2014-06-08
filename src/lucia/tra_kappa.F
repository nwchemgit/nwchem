*----------------------------------------------------------------------*
      SUBROUTINE MKNEWCMO(LUKAP,LUUMAT,LUCMO,IOXCC,NOXC,NSPIN)
*
*     On input  :  LUKAP    kappa   OR
*                  LUMAT    U
*                  LUCMO    old CMO (if negative number is provided:
*                             read from environment, e.g. SIRIFC)
*
*     On output :  LUCMO    new CMO = CMO(old)*U
*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cintfo.inc'
      
      INTEGER, PARAMETER :: NTEST = 00

      DIMENSION IOXCC(*), NOXC(2)

      REAL(8) INPROD

      LBLK = -1
      IF (NTEST.GE.10) THEN
        WRITE(6,*) 'MKNEWCMO AT WORK'
        WRITE(6,*) '================'
        WRITE(6,*) ' NSPIN: ',NSPIN
        WRITE(6,*) ' NOXC:  ', NOXC(1:NSPIN)
        WRITE(6,*) ' LUKAP, LUUMAT, LUCMO: ',
     &               LUKAP, LUUMAT, LUCMO
      END IF

      NOXC_TOT = NOXC(1)
      IF (NSPIN.EQ.2) NOXC_TOT = NOXC_TOT+NOXC(2)

      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MKCMO ')

      LEN = 0
      LCMO = 0
      DO ISM = 1, NSMOB
        LEN = LEN + NTOOBS(ISM)*NTOOBS(ISM)
        LCMO = LCMO + NAOS_ENV(ISM)*NTOOBS(ISM)
      END DO

      LEN1 = MAX(LEN,NOXC_TOT)
      LEN2 = MAX(LEN,LCMO)
      CALL MEMMAN(KLBUF1,LEN1,'ADDL  ',2,'BUFF1 ')
      CALL MEMMAN(KLBUF2,LEN2,'ADDL  ',2,'BUFF2 ')
      CALL MEMMAN(KLCMO,LCMO,'ADDL  ',2,'CMO   ')

      IF (LUUMAT.LE.0) THEN
        CALL MEMMAN(KKAP,NOXC_TOT,'ADDL  ',2,'KAPPA ')
        CALL VEC_FROM_DISC(WORK(KKAP),NOXC_TOT,1,LBLK,LUKAP)
      END IF

      ! note: no spin-dependence for input CMO implemented
      IF (LUCMO.GT.0) THEN
        CALL VEC_FROM_DISC(WORK(KLCMO),LCMO,1,LBLK,LUCMO)
      ELSE
        CALL GET_CMOAO_ENV(WORK(KLCMO))
      END IF

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'Old CMO matrix'
        CALL APRBLM2(WORK(KLCMO),NAOS_ENV,NTOOBS,NSMOB,0)          
      END IF

      DO ISPIN = 1, NSPIN
        KOFF  = (ISPIN-1)*NOXC(1)
        IXOFF = (ISPIN-1)*NOXC(1)*2 + 1

        IF (LUUMAT.LE.0) THEN
          CALL EXPAND_KAP(WORK(KKAP+KOFF),WORK(KLBUF2),
     &                                 IOXCC(IXOFF),NOXC(ISPIN))
          CALL EXPT_KAPPA(WORK(KLBUF2),WORK(KLBUF1))
        ELSE
          IF (ISPIN.EQ.1) IREW = 1
          IF (ISPIN.EQ.2) IREW = 0
          CALL VEC_FROM_DISC(WORK(KLBUF1),LEN,IREW,LBLK,LUUMAT)
        END IF

        ! KLCMO x KLBUF1 --> KLBUF2
        IOFF1 = 0
        IOFF2 = 0
        DO ISM = 1, NSMOB
          NDIM1 = NTOOBS(ISM)
          NDIM2 = NAOS_ENV(ISM)
          CALL MATML7(WORK(KLBUF2+IOFF2),WORK(KLCMO+IOFF2),
     &         WORK(KLBUF1+IOFF1),
     &         NDIM2,NDIM1,NDIM2,NDIM1,NDIM1,NDIM1,0d0,1d0,0)
          IOFF1 = IOFF1 + NDIM1*NDIM1
          IOFF2 = IOFF2 + NDIM1*NDIM2
        END DO
 
        ! save new U
        ! alpha and beta part on consecutive records
        IF (ISPIN.EQ.1) IREW = 1
        IF (ISPIN.EQ.2) IREW = 0
        CALL VEC_TO_DISC(WORK(KLBUF2),LCMO,IREW,LBLK,ABS(LUCMO))

        IF (NTEST.GE.100) THEN
          WRITE(6,*) 'New CMO matrix (ISPIN=',ISPIN,')'
          CALL APRBLM2(WORK(KLBUF2),NAOS_ENV,NTOOBS,NSMOB,0)          
        END IF

      END DO ! ISPIN

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'MKCMO ')

      IF (NTEST.GE.10) WRITE(6,*) 'LEAVING MKCMO'

      RETURN
      END
*----------------------------------------------------------------------*
      SUBROUTINE KAP2U(IMODE,LUDKAP,LUKAP,LUUMAT,IOXCC,NOXC,NSPIN)
*
*     On input  :  LUDKAP   d_kappa
*                  LUKAP    kappa_old
*                  LUMAT    Uold
*
*     On output :  LUDKAP   step_kappa = kappa - kappa_old
*                  LUKAP    kappa
*                  LUMAT    Unew
*
*      where Unew = exp(d_kappa) Uold
*            kappa = ln(Unew)
*            step_kappa = kappa_new - kappa_old
*
*    imode = 1   full program
*    imode = 2   do not save Unew
*    imode = 3   just get Unew from kappa (!) and save
*    imode = -3  just get Unew from dkappa (!) and save
*    imode = 4   do not save kappa_new and Unew
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cintfo.inc'
      
      INTEGER, PARAMETER :: NTEST = 00

      DIMENSION IOXCC(*), NOXC(2)

      REAL(8) INPROD

      LBLK = -1
      IF (NTEST.GE.10) THEN
        WRITE(6,*) 'KAP2U AT WORK'
        WRITE(6,*) '============='
        WRITE(6,*) ' IMODE, NSPIN: ',IMODE, NSPIN
        WRITE(6,*) ' NOXC: ', NOXC(1:NSPIN)
        WRITE(6,*) ' LUDKAP, LUKAP, LUUMAT: ',
     &               LUDKAP, LUKAP, LUUMAT
      END IF

      NOXC_TOT = NOXC(1)
      IF (NSPIN.EQ.2) NOXC_TOT = NOXC_TOT+NOXC(2)

      IF (NSPIN.EQ.2) LUSCR = IOPEN_NUS('KAP2U_SCR')

      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'KAP2U ')

      LEN = 0
      DO ISM = 1, NSMOB
        LEN = LEN + NTOOBS(ISM)*NTOOBS(ISM) 
      END DO

      LEN1 = MAX(LEN,NOXC_TOT)
      CALL MEMMAN(KLBUF1,LEN1,'ADDL  ',2,'BUFF1 ')
      CALL MEMMAN(KLBUF2,LEN,'ADDL  ',2,'BUFF2 ')

      CALL MEMMAN(KKAP,NOXC_TOT,'ADDL  ',2,'KAPPA ')

      IF (IMODE.EQ.3) THEN
        CALL VEC_FROM_DISC(WORK(KKAP),NOXC_TOT,1,LBLK,LUKAP)
      ELSE
        CALL VEC_FROM_DISC(WORK(KKAP),NOXC_TOT,1,LBLK,LUDKAP)
      END IF

      DO ISPIN = 1, NSPIN
        KOFF  = (ISPIN-1)*NOXC(1)
        IXOFF = (ISPIN-1)*NOXC(1)*2 + 1
        CALL EXPAND_KAP(WORK(KKAP+KOFF),WORK(KLBUF2),
     &                               IOXCC(IXOFF),NOXC(ISPIN))
        CALL EXPT_KAPPA(WORK(KLBUF2),WORK(KLBUF1))

        IF (IMODE.NE.3) THEN
          CALL MEMMAN(KLUM,LEN,'ADDL  ',2,'UMAT  ')
          ! load U to KLUM
          IF (ISPIN.EQ.1)
     &       CALL VEC_FROM_DISC(WORK(KLUM),LEN,1,LBLK,LUUMAT)
          ! park beta part of Umatrix on LUSCR ...
          IF (ISPIN.EQ.1.AND.NSPIN.EQ.2) THEN
            CALL REWINO(LUSCR)
            CALL COPVCD(LUUMAT,LUSCR,WORK(KLBUF2),0,LBLK)
          END IF
          ! ... and fetch it from there in the second round
          IF (ISPIN.EQ.2)
     &       CALL VEC_FROM_DISC(WORK(KLUM),LEN,1,LBLK,LUSCR)


          ! KLUM x KLBUF1 --> KLBUF2
          IOFF = 0
          DO ISM = 1, NSMOB
            NDIM = NTOOBS(ISM)
            CALL MATML7(WORK(KLBUF2+IOFF),WORK(KLBUF1+IOFF),
     &         WORK(KLUM+IOFF),
     &         NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0d0,1d0,0)
            IOFF = IOFF + NDIM*NDIM
          END DO
 
        END IF
        
        ! save new U
        ! alpha and beta part on consecutive records
        IF (ISPIN.EQ.1) IREW = 1
        IF (ISPIN.EQ.2) IREW = 0
        IF (IMODE.EQ.1.OR.IMODE.EQ.-3) THEN
c          CALL VEC_TO_DISC(WORK(KLUM),LEN,IREW,LBLK,LUUMAT)
          CALL VEC_TO_DISC(WORK(KLBUF2),LEN,IREW,LBLK,LUUMAT)
        ELSE IF (IMODE.EQ.3) THEN
          CALL VEC_TO_DISC(WORK(KLBUF1),LEN,IREW,LBLK,LUUMAT)
        END IF

        IF (ABS(IMODE).NE.3) THEN
          ! kappa_new = ln(U)
          CALL LOG_UMAT(WORK(KLBUF2),WORK(KLBUF1))

          CALL COMPRS_KAP(WORK(KLBUF1),WORK(KKAP+KOFF),
     &                                 IOXCC(IXOFF),NOXC(ISPIN))
        END IF

      END DO ! ISPIN

      IF (ABS(IMODE).NE.3) THEN

        ! get old kappa
        CALL VEC_FROM_DISC(WORK(KLBUF1),NOXC_TOT,1,LBLK,LUKAP)

        CALL VECSUM(WORK(KLBUF1),WORK(KLBUF1),WORK(KKAP),
     &              -1d0,1d0,NOXC_TOT)

        ! save new kappa and difference to old kappa
        CALL VEC_TO_DISC(WORK(KLBUF1),NOXC_TOT,1,LBLK,LUDKAP)
        IF (IMODE.NE.4)
     &       CALL VEC_TO_DISC(WORK(KKAP),NOXC_TOT,1,LBLK,LUKAP)

      END IF

      IF (NSPIN.EQ.2) CALL RELUNIT(LUSCR,'delete')

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'KAP2U ')

      IF (NTEST.GE.10) WRITE(6,*) 'LEAVING KAP2U'

      RETURN
      END
      SUBROUTINE TRA_KAPPA(LUKAP,LUUMAT,IOXCC,NOXC,NSPIN,
     &                     IRELOAD,LU1INT,LU2INT)
*
*     from TRA_SIMH
*     AK, April 2004
*
*     primitive MO -> MO transformation for
*
*     Hnew = exp(-kappa)H exp(kappa)
*
*     if LUUMAT .gt. 0 :  Hnew = U^+HU
*
*     LUKAP -> XKAP(NOXC) is coming along in compressed form, IOXC contains the
*     in information to decompress.
*
* NSPIN
* = 1 => orbital transformation assuming that alpha and beta parts of kappa
*        are identical
* = 2 => spin-orbital transformation allowing differences in kappa(alpha)
*        and kappa(beta)
*
* IRELOAD == 1: reload original 1ints and 2ints from file 
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      DIMENSION IOXCC(*), NOXC(NSPIN)
      REAL(8) INPROD

*
      CALL ATIM(CPU0,WALL0)

      IDUM = 0
      NTEST = 00
      
      NOXC_TOT = NOXC(1)
      IF (NSPIN.EQ.2) NOXC_TOT = NOXC_TOT + NOXC(2)

      LBLK = -1

      ICC_EXC_M = ICC_EXC
      ICC_EXC = 0

      CALL LUCIAQENTER('KAPTR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRAKAP')
*
      IF(NSPIN.EQ.1) THEN
        WRITE(6,*) ' Orbital based transformation '
      ELSE 
        WRITE(6,*) ' Spinorbital based transformation '
      END IF

      IF(NTEST.GE.10) THEN
        IF (LUUMAT.LE.0) THEN
          WRITE(6,*) ' U is built from kappa (LUKAP = ',LUKAP,')'
        ELSE
          WRITE(6,*) ' U is read from disc (LUUMAT = ',LUUMAT,')'        
        END IF
      END IF

      I_UNRORB = 0 ! getint should only fetch KINT2 integrals
      IF (I_USE_SIMTRH.EQ.1) THEN
        WRITE(6,*) 'Boy, you are in the wrong trafo-routine!'
        STOP 'TRA_KAPPA'
      END IF

      LEN = 0
      DO ISM = 1, NSMOB
        LEN = LEN + NTOOBS(ISM)*NTOOBS(ISM) 
      END DO
      CALL MEMMAN(KLUM,LEN,'ADDL  ',2,'UTRA  ')
      IF (NSPIN.EQ.2)
     &     CALL MEMMAN(KLUMB,LEN,'ADDL  ',2,'UTRA B')

      IF (LUUMAT.LE.0) THEN
        ! work memory for exp(kappa)
        IDUM = 0
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPKAP')
        CALL MEMMAN(KLKAP,NOXC_TOT,'ADDL  ',2,'KAPF ')
        CALL MEMMAN(KLBUF1,LEN,'ADDL  ',2,'BUFF1 ')

        ! load KAPPA
        CALL VEC_FROM_DISC(WORK(KLKAP),NOXC_TOT,1,-1,LUKAP)
      END IF

      DO ISPIN = 1, NSPIN
        
        IF (ISPIN.EQ.1) KU_ADR = KLUM
        IF (ISPIN.EQ.2) KU_ADR = KLUMB

        KOFF = (ISPIN-1)*NOXC(1)
        IXOFF = (ISPIN-1)*NOXC(1)*2+1

        IF (LUUMAT.LE.0) THEN
          ! Obtain KAPPA in matrix-form
          CALL EXPAND_KAP(WORK(KLKAP+KOFF),WORK(KLBUF1),
     &                    IOXCC(IXOFF),NOXC(ISPIN))
          ! Obtain transformation matrix U = exp(kappa)
          CALL EXPT_KAPPA(WORK(KLBUF1),WORK(KU_ADR))
        ELSE
          LBLK = -1
          IF (ISPIN.EQ.1) IREW = 1
          IF (ISPIN.EQ.2) IREW = 0
          CALL VEC_FROM_DISC(WORK(KU_ADR),LEN,IREW,LBLK,LUUMAT)
        END IF

        IF (NTEST.GE.100) THEN
          WRITE(6,*) 'U-matrix:'
          CALL APRBLM2(WORK(KU_ADR),NTOOBS,NTOOBS,NSMOB,0)
        END IF
      END DO

      IF (LUUMAT.LE.0) THEN
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXPKAP')
      END IF

*----------------------------------------------------------------------*
* Transform one-electron integrals
*----------------------------------------------------------------------*
      IF(IRELOAD.EQ.1) THEN
        CALL VEC_FROM_DISC(WORK(KINT1),NINT1,1,LBLK,LU1INT)
        IF (NSPIN.EQ.2) CALL COPVEC(WORK(KINT1),WORK(KINT1B),NINT1)
      END IF

      DO ISPIN = 1, NSPIN
        IF (ISPIN.EQ.1) THEN
          KU_ADR = KLUM
          KI1_ADR = KINT1
        ELSE
          KU_ADR = KLUMB
          KI1_ADR = KINT1B
        END IF

        CALL TRA1_SIMPLE2(WORK(KU_ADR),WORK(KI1_ADR))
      END DO

*----------------------------------------------------------------------*
* Transform two-electron integrals
*----------------------------------------------------------------------*
      IF (NSPIN.EQ.1) NSPC = 1
      IF (NSPIN.EQ.2) NSPC = 3

      IF (IRELOAD.EQ.1)
     &     CALL VEC_FROM_DISC(WORK(KINT2),NINT2,1,LBLK,LU2INT)

      ! do the loop reversly, as KINT2 is implicitely assumed as
      ! source of the untransformed integrals; as currently the
      ! transformation proceeds symmetry-block by symmetry-block
      ! (well....) the last transformation may write the transformed
      ! integrals onto its source array
      DO ISPC = NSPC, 1, -1
        IF (ISPC.EQ.1) THEN
          KU_ADR1 = KLUM
          KU_ADR2 = KLUM
          KI2_ADR = KINT2
        ELSE IF (ISPC.EQ.2) THEN
          KU_ADR1 = KLUMB
          KU_ADR2 = KLUMB
          KI2_ADR = KINT2BB
        ELSE IF (ISPC.EQ.3) THEN
          KU_ADR1 = KLUM
          KU_ADR2 = KLUMB
          KI2_ADR = KINT2AB
          KPI2_ADR = KPINT2AB
        END IF

        IF (ISPC.EQ.1.OR.ISPC.EQ.2)
     &       CALL TRA2_SIMPLE(WORK(KU_ADR1),WORK(KI2_ADR))

        IF (ISPC.EQ.3)
     &       CALL TRA2Q_SIMPLE(WORK(KU_ADR1),WORK(KU_ADR1),
     &                         WORK(KU_ADR2),WORK(KU_ADR2),
     &                         1,0,WORK(KI2_ADR),WORK(KPI2_ADR))

      END DO

*----------------------------------------------------------------------*
* hole contributions to core energy and one-electron integrals
*----------------------------------------------------------------------*
      IF(IUSE_PH.EQ.1) THEN
        IF(NSPIN.EQ.1) THEN
          CALL FI(WORK(KINT1),ECORE_HEX,1)
        ELSE 
          I_UNRORB = 1
          CALL FI_HS_AB(WORK(KINT1),WORK(KINT1B),
     &                  WORK(KFI_AL),WORK(KFI_BE),ECORE_HEX,1) 
        END IF
      ELSE
         ECORE_HEX = 0.0D0
      END IF

      ECORE = ECORE_ORIG + ECORE_HEX
      WRITE(6,*) ' TRA_KAPPA : Updated core energy ',ECORE
     
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRAKAP')
      CALL LUCIAQEXIT('KAPTR')

      CALL ATIM(CPU,WALL)
      CALL PRTIM(6,'time in transformation',CPU-CPU0,WALL-WALL0)

      IF (NSPIN.EQ.2) I_UNRORB = 1
      ICC_EXC = ICC_EXC_M

      RETURN
      END
      SUBROUTINE EXPAND_KAP(XKAP_CMP,XKAP_FULL,IOXCC,NOXC)!,NTOOB)

      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'

      INTEGER, PARAMETER :: NTEST = 00

      DIMENSION XKAP_CMP(*), XKAP_FULL(*)
      DIMENSION IOXCC(2,*)

      INTEGER IOFF(NSMOB)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPKAP')
      CALL MEMMAN(KLBUF,NTOOB*NTOOB,'ADDL  ',2,'BUFFA ')

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'KAPPA in compressed form:'
        DO IEXC = 1, NOXC
          WRITE(6,'(X,2I6,E20.10)')
     &         IOXCC(1,IEXC), IOXCC(2,IEXC), XKAP_CMP(IEXC)
        END DO
      END IF

      IDX = 0
      DO ISM = 1, NSMOB
        IOFF(ISM) = IDX
        XKAP_FULL(IDX+1:IDX+NTOOBS(ISM)*NTOOBS(ISM)) = 0d0
        IDX = IDX + NTOOBS(ISM)*NTOOBS(ISM)        
      END DO

      DO IEXC = 1, NOXC
        II = IOXCC(1,IEXC)
        JJ = IOXCC(2,IEXC)

        ISM = ISMFTO(II)
        II = IREOTS(II) - IBSO(ISM) + 1
        JJ = IREOTS(JJ) - IBSO(ISM) + 1
        IJDX = IOFF(ISM) + (JJ-1)*NTOOBS(ISM)+II
        JIDX = IOFF(ISM) + (II-1)*NTOOBS(ISM)+JJ
        XKAP_FULL(IJDX) =  XKAP_CMP(IEXC)
        XKAP_FULL(JIDX) = -XKAP_CMP(IEXC)
      END DO

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'kappa in blocked matrix form:'
        CALL APRBLM2(XKAP_FULL,NTOOBS,NTOOBS,NSMOB,0)
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'EXPKAP')
      
      RETURN
      END
      SUBROUTINE COMPRS_KAP(XKAP_FULL,XKAP_CMP,IOXCC,NOXC)

      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'

      INTEGER, PARAMETER :: NTEST = 00

      DIMENSION XKAP_CMP(*), XKAP_FULL(*)
      DIMENSION IOXCC(2,*)

      INTEGER IOFF(NSMOB)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPKAP')
      CALL MEMMAN(KLBUF,NTOOB*NTOOB,'ADDL  ',2,'BUFFA ')

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'kappa in blocked matrix form:'
        CALL APRBLM2(XKAP_FULL,NTOOBS,NTOOBS,NSMOB,0)
      END IF

      IDX = 0
      DO ISM = 1, NSMOB
        IOFF(ISM) = IDX
        IDX = IDX + NTOOBS(ISM)*NTOOBS(ISM)        
      END DO

      DO IEXC = 1, NOXC
        II = IOXCC(1,IEXC)
        JJ = IOXCC(2,IEXC)

        ISM = ISMFTO(II)
        II = IREOTS(II) - IBSO(ISM) + 1
        JJ = IREOTS(JJ) - IBSO(ISM) + 1
        IJDX = IOFF(ISM) + (JJ-1)*NTOOBS(ISM)+II
        XKAP_CMP(IEXC) = XKAP_FULL(IJDX) 
      END DO

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'KAPPA in compressed form:'
        DO IEXC = 1, NOXC
          WRITE(6,'(X,2I6,E20.10)')
     &         IOXCC(1,IEXC), IOXCC(2,IEXC), XKAP_CMP(IEXC)
        END DO
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'EXPKAP')
      
      RETURN
      END
      SUBROUTINE EXPT_KAPPA(XKAP,XEXPKAP)!,NTOOB)

      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'constants.inc'

      INTEGER, PARAMETER :: NTEST = 00
      DIMENSION XKAP(*), XEXPKAP(*)
      LOGICAL CONVERGED
      REAL(8), EXTERNAL :: INPROD

c      THRSH = 1D-12
      THRSH = 1D-20

      IF (NTEST.GT.100) THEN
        WRITE(6,*) 'On input: kappa = '
        CALL APRBLM2(XKAP,NTOOBS,NTOOBS,NSMOB,0)
      END IF

      ! Init with 1 + kappa
      NTOT = 0
      NMAX = 0
      IOFF = 0
      DO ISM = 1, NSMOB
        DO II = 1, NTOOBS(ISM)*NTOOBS(ISM)
          XEXPKAP(IOFF+II) = XKAP(IOFF+II)
        END DO
        IDX = 1
        DO II = 1, NTOOBS(ISM)
          XEXPKAP(IOFF+IDX) = XEXPKAP(IOFF+IDX) + 1D0
          IDX = IDX+NTOOBS(ISM)+1
        END DO
        NTOT = NTOT + NTOOBS(ISM)*NTOOBS(ISM)
        NMAX = MAX(NMAX,NTOOBS(ISM)*NTOOBS(ISM))
        IOFF = IOFF + NTOOBS(ISM)*NTOOBS(ISM)
      END DO

      XNRM = SQRT(INPROD(XKAP,XKAP,NTOT))


      IF (XNRM.LT.THRSH) RETURN

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPTKP')
      CALL MEMMAN(KLBUF,NMAX,'ADDL  ',2,'BUFFA ')
      CALL MEMMAN(KLUUU,NMAX,'ADDL  ',2,'BUFFB ')

      IOFF = 0
      DO ISM = 1, NSMOB

        IF (NTEST.GE.10)
     &         WRITE(6,*) 'Symmetry block ',ISM

        ! Init U^n array
        NDIM = NTOOBS(ISM)
        IF (NDIM.EQ.0) EXIT

        NDIM2 = NDIM*NDIM
        DO II = 1, NDIM2
          WORK(KLUUU-1+II) = XKAP(IOFF+II)
        END DO
        FAC = 1D0
        N = 1
        MAXN = 100

        XNRM = 100d0 ! just init
        
        DO WHILE (XNRM.GE.THRSH)
          N = N+1
          FAC = 1d0/DBLE(N)
          IF (N.GT.MAXN) THEN
            WRITE(6,*) 'EXPT_KAPPA did not converge in ',MAXN,
     &           ' iterations for symmetry block ',ISM
            stop 'problem in expt_kappa'
          END IF

          CALL MATML7(WORK(KLBUF),WORK(KLUUU),XKAP(IOFF+1),
     &              NDIM,NDIM,
     &              NDIM,NDIM,
     &              NDIM,NDIM,0d0,FAC,0)

          XNRM = SQRT(INPROD(WORK(KLBUF),WORK(KLBUF),NDIM2))
        
          IF (NTEST.GE.10)
     &         WRITE(6,*) ' N = ',N,'  |1/N! U^N| = ',XNRM

          CALL VECSUM(XEXPKAP(IOFF+1),XEXPKAP(IOFF+1),
     &                WORK(KLBUF),1d0,1d0,NDIM*NDIM)

          CALL COPVEC(WORK(KLBUF),WORK(KLUUU),NDIM2)

        END DO

        IOFF = IOFF + NDIM2

      END DO ! ISM

      IF (NTEST.EQ.100) THEN
        WRITE(6,*) 'exp(k) in blocked matrix form:'
        CALL APRBLM2(XEXPKAP,NTOOBS,NTOOBS,NSMOB,0)
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'EXPTKP')
      
      RETURN
      END
      SUBROUTINE LOG_UMAT(XUMAT,XLOGUMAT)

      ! note XUMAT will be XUMAT - 1 on exit
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'

      INTEGER, PARAMETER :: NTEST = 00
      DIMENSION XUMAT(*), XLOGUMAT(*)
      LOGICAL CONVERGED
      REAL(8), EXTERNAL :: INPROD

      THRSH = 1D-12

      IF (NTEST.EQ.100) THEN
        WRITE(6,*) 'on entry in log_umat:'
        WRITE(6,*) 'U in blocked matrix form:'
        CALL APRBLM2(XUMAT,NTOOBS,NTOOBS,NSMOB,0)
      END IF

c      ! Init with U - 1
c      NTOT = 0
      NMAX = 0
c      IOFF = 0
c      XMAX = 0D0
      DO ISM = 1, NSMOB
c        IDX = 1
c        DO II = 1, NTOOBS(ISM)
c          XUMAT(IOFF+IDX) = XUMAT(IOFF+IDX) - 1D0
c          IDX = IDX+NTOOBS(ISM)+1
c        END DO
c        XMAX =
c     &       MAX(XMAX,FNDMNX(XUMAT(IOFF+1),NTOOBS(ISM)*NTOOBS(ISM),2))
c        DO II = 1, NTOOBS(ISM)*NTOOBS(ISM)
c          XLOGUMAT(IOFF+II) = XUMAT(IOFF+II)
c        END DO
c        NTOT = NTOT + NTOOBS(ISM)*NTOOBS(ISM)
        NMAX = MAX(NMAX,NTOOBS(ISM)*NTOOBS(ISM))
c        IOFF = IOFF + NTOOBS(ISM)*NTOOBS(ISM)
      END DO
c
c      XNRM = SQRT(INPROD(XLOGUMAT,XLOGUMAT,NTOT))
c
c      IF (XNRM.LT.THRSH) RETURN


      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LOG U ')
      CALL MEMMAN(KLBUF,NMAX,'ADDL  ',2,'BUFFA ')
      CALL MEMMAN(KLBUG,NMAX,'ADDL  ',2,'BUFFB ') ! no bug, hopefully
      CALL MEMMAN(KLUUU,NMAX,'ADDL  ',2,'BUFFC ')


      IOFF = 0
      DO ISM = 1, NSMOB

        IF (NTEST.GE.10)
     &         WRITE(6,*) 'Symmetry block ',ISM

        ! Init U^n array
        NDIM = NTOOBS(ISM)
        IF (NDIM.EQ.0) EXIT

        NDIM2 = NDIM*NDIM

c new: via eigenvalue spectrum of U
        CALL LOGUMAT(NDIM,XLOGUMAT(IOFF+1),XUMAT(IOFF+1),
     &       WORK(KLBUF),WORK(KLUUU),WORK(KLBUG))

c old: Taylor expansion (very limited convergence radius!)
c        DO II = 1, NDIM2
c          WORK(KLUUU-1+II) = XUMAT(IOFF+II)
c        END DO
c        FAC1 = 1D0
c        N = 1
c        MAXN = 100
c
c        XNRM = 100d0 ! just init
c        
c        DO WHILE (XNRM.GE.THRSH)
c          N = N+1
c          FAC1 = FAC1*(-1D0)
c          FAC2 = 1d0/DBLE(N)
c          IF (N.GT.MAXN) THEN
c            WRITE(6,*) 'LOG_UMAT did not converge in ',MAXN,
c     &           ' iterations for symmetry block ',ISM
c            WRITE(6,*) 'largest element of U: ', XMAX
c            stop 'problem in log_umat'
c          END IF
c
c          CALL MATML7(WORK(KLBUF),WORK(KLUUU),XUMAT(IOFF+1),
c     &              NDIM,NDIM,
c     &              NDIM,NDIM,
c     &              NDIM,NDIM,0D0,1D0,0)
c
c          XNRM = SQRT(INPROD(WORK(KLBUF),WORK(KLBUF),NDIM2))
c        
c          IF (NTEST.GE.10)
c     &         WRITE(6,*) ' N = ',N,'  |1/N U^N| = ',XNRM
c
c          CALL VECSUM(XLOGUMAT(IOFF+1),XLOGUMAT(IOFF+1),
c     &                WORK(KLBUF),1d0,FAC1*FAC2,NDIM2)
c
c          CALL COPVEC(WORK(KLBUF),WORK(KLUUU),NDIM2)
c
c        END DO

        IOFF = IOFF + NDIM2

      END DO ! ISM

      IF (NTEST.EQ.100) THEN
        WRITE(6,*) 'log(U) in blocked matrix form:'
        CALL APRBLM2(XLOGUMAT,NTOOBS,NTOOBS,NSMOB,0)
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'LOG U ')
      
      RETURN
      END
      SUBROUTINE COMPRESS_T1S(T,T1,ITSM,IC_OCC,IA_OCC,
     &                ISPOBEX_TP,ISOX_TO_OX,IEXC_FOR_OX,
     &                IBSPOBEX_TP,NSPOBEX_TP,IAB,ISYMMET)
*
* Slave routine for compressing the p/h part of complete matrix T1
* to singles amplitudes in the string format stored T
*
* IAB = 1 => alpha excitations are collected 
* IAB = 2 => beta excitations are collected 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cc_exc.inc'
*. General input 
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISOX_TO_OX(NSPOBEX_TP)
      INTEGER IEXC_FOR_OX(*)
      INTEGER IBSPOBEX_TP(NSPOBEX_TP)
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION T1(*)
*. Scratch through input
      INTEGER IC_OCC(NTOOB),IA_OCC(NTOOB)
*. Local scratch
      INTEGER IC_GRP(MXPNGAS),IA_GRP(MXPNGAS),IOFF(8)
*

      IF (ITSM.NE.1) THEN
        WRITE(6,*) 'ITSM.NE.1 in COMPRESS_T1S'
        STOP 'EXPAND_T1S'
      END IF

      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'input T1 in blocked matrix form '
        CALL APRBLM2(T1,NTOOBS,NTOOBS,NSMOB,0)
      END IF

      IDX = 0
      DO ISM = 1, NSMOB
        IOFF(ISM) = IDX
        IDX = IDX + NTOOBS(ISM)*NTOOBS(ISM)        
      END DO
      LEN = IDX

      IF(MSCOMB_CC.EQ.1) THEN
*. Combinations, 
        FACTOR = 1.0D0/SQRT(2.0D0)
      ELSE
*. Individual terms
        FACTOR = 1.0D0
      END IF
*
      DO ISOXTP = 1, NSPOBEX_TP
*
      IF(IEXC_FOR_OX(ISOX_TO_OX(ISOXTP)).EQ.1) THEN
*. Is single excitation alpha or beta excitation ?
        NCA_OP = IELSUM(ISPOBEX_TP(1+0*NGAS,ISOXTP),NGAS)
        NCB_OP = IELSUM(ISPOBEX_TP(1+1*NGAS,ISOXTP),NGAS)
        NAA_OP = IELSUM(ISPOBEX_TP(1+2*NGAS,ISOXTP),NGAS)
        NAB_OP = IELSUM(ISPOBEX_TP(1+3*NGAS,ISOXTP),NGAS)
        IELMNT = IBSPOBEX_TP(ISOXTP)
        IF(NCA_OP.EQ.1.AND.IAB.EQ.1.OR.
     &     NCB_OP.EQ.1.AND.IAB.EQ.2     )THEN
          IF(NCA_OP.EQ.1.AND.IAB.EQ.1) THEN
*. Transform from occupations to groups
            CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ISOXTP),IC_GRP,1)
            CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ISOXTP),IA_GRP,1)
            NC_OP = NCA_OP
            NA_OP = NAA_OP
          ELSE IF(NCB_OP.EQ.1.AND.IAB.EQ.2) THEN
*. Transform from occupations to groups
            CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ISOXTP),IC_GRP,1)
            CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ISOXTP),IA_GRP,1)
            NC_OP = NCB_OP
            NA_OP = NAB_OP
          END IF
*
          DO ISM_C = 1, NSMST
            ISM_A = MULTD2H(ISM_C,ITSM)
*. Obtain occupations
            IDUM = 0
            CALL GETSTR2_TOTSM_SPGP(IC_GRP,NGAS,ISM_C,NC_OP,NSTR_C,
     &           IC_OCC, NTOOB,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IA_GRP,NGAS,ISM_A,NA_OP,NSTR_A,
     &           IA_OCC, NTOOB,0,IDUM,IDUM)
*. And read from complete matrix 
C?           WRITE(6,*) ' ISM_C, ISM_A, NSTR_A, NSTR_C =', 
C?   &                    ISM_C, ISM_A, NSTR_A, NSTR_C 
            IF (ISYMMET.EQ.0) THEN
              DO ISA = 1, NSTR_A
                DO ISC = 1, NSTR_C
                 IORB_A = IA_OCC(ISA)
                 IORB_C = IC_OCC(ISC)
                 II = IREOTS(IORB_C) - IBSO(ISM_C) + 1
                 JJ = IREOTS(IORB_A) - IBSO(ISM_A) + 1
                 IJDX = IOFF(ISM_A) + (JJ-1)*NTOOBS(ISM_A)+II
                 T(IELMNT) = T1(IJDX)/FACTOR
                 IELMNT = IELMNT + 1
                END DO
              END DO
            ELSE IF (ISYMMET.EQ.1) THEN
              DO ISA = 1, NSTR_A
                DO ISC = 1, NSTR_C
                 IORB_A = IA_OCC(ISA)
                 IORB_C = IC_OCC(ISC)
                 II = IREOTS(IORB_C) - IBSO(ISM_C) + 1
                 JJ = IREOTS(IORB_A) - IBSO(ISM_A) + 1
                 IJDX = IOFF(ISM_A) + (JJ-1)*NTOOBS(ISM_A)+II
                 JIDX = IOFF(ISM_C) + (II-1)*NTOOBS(ISM_C)+JJ
                 T(IELMNT) = (T1(IJDX)+T1(JIDX))*(0.5d0/FACTOR)
                 IELMNT = IELMNT + 1
                END DO
              END DO
            ELSE IF (ISYMMET.EQ.-1) THEN
              DO ISA = 1, NSTR_A
                DO ISC = 1, NSTR_C
                 IORB_A = IA_OCC(ISA)
                 IORB_C = IC_OCC(ISC)
                 II = IREOTS(IORB_C) - IBSO(ISM_C) + 1
                 JJ = IREOTS(IORB_A) - IBSO(ISM_A) + 1
                 IJDX = IOFF(ISM_A) + (JJ-1)*NTOOBS(ISM_A)+II
                 JIDX = IOFF(ISM_C) + (II-1)*NTOOBS(ISM_C)+JJ
                 T(IELMNT) = (T1(IJDX)-T1(JIDX))*(0.5d0/FACTOR)
                 IELMNT = IELMNT + 1
                END DO
              END DO
            ELSE IF (ISYMMET.EQ.-2) THEN
              DO ISA = 1, NSTR_A
                DO ISC = 1, NSTR_C
                 IORB_A = IA_OCC(ISA)
                 IORB_C = IC_OCC(ISC)
                 II = IREOTS(IORB_C) - IBSO(ISM_C) + 1
                 JJ = IREOTS(IORB_A) - IBSO(ISM_A) + 1
                 IJDX = IOFF(ISM_A) + (JJ-1)*NTOOBS(ISM_A)+II
                 JIDX = IOFF(ISM_C) + (II-1)*NTOOBS(ISM_C)+JJ
                 T(IELMNT) = (T1(IJDX)-T1(JIDX))*(2.0d0)
c                 T(IELMNT) = (T1(IJDX)-T1(JIDX))*(2.0d0/FACTOR)
                 IELMNT = IELMNT + 1
                END DO
              END DO
            END IF
          END DO
        END IF
*       ^ End of correct single excitation
      END IF
*     ^ End if single excitation
      END DO
*     ^ End of loop over excitation types
*
      RETURN
      END
      SUBROUTINE ZERO_T1(T)
*
* Convenience-wrapper for
* routine zeroing t1 coefficients
*
      INCLUDE 'wrkspc.inc'
c      INCLUDE 'orbinp.inc'
      INCLUDE 'ctcc.inc'
c      INCLUDE 'csm.inc'

      CALL ZERO_T1S(T,
     &     WORK(KLSOBEX),WORK(KLSOX_TO_OX),WORK(KLCOBEX_TP),
     &     WORK(KLIBSOBEX),WORK(KLLSOBEX),NSPOBEX_TP)      

      END
      SUBROUTINE ZERO_T1S(T,
     &                ISPOBEX_TP,ISOX_TO_OX,IEXC_FOR_OX,
     &                IBSPOBEX_TP,ILSPOBEX_TP,NSPOBEX_TP)
*
* Routine for zeroing t1 coefficients
*
* ITASK = 0 : SET TO ZERO 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cc_exc.inc'
*. General input 
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISOX_TO_OX(NSPOBEX_TP)
      INTEGER IEXC_FOR_OX(*)
      INTEGER IBSPOBEX_TP(NSPOBEX_TP)
      INTEGER ILSPOBEX_TP(NSPOBEX_TP)
*. Specific input
      DIMENSION T(*)

      DO ISOXTP = 1, NSPOBEX_TP
*
*. Is single excitation alpha or beta excitation ?
        IF(IEXC_FOR_OX(ISOX_TO_OX(ISOXTP)).EQ.1) THEN
*
          IOFF = IBSPOBEX_TP(ISOXTP)
          NDIM = ILSPOBEX_TP(ISOXTP)

          T(IOFF:IOFF-1+NDIM) = 0D0
*
        END IF
      END DO
*
      NTEST = 00
      IF (NTEST.GE.1000) THEN
        WRITE(6,*) 'T after 0ing T1 parts:'
        CALL WRT_CC_VEC2(T,6,'GEN_CC')
      END IF
*
      RETURN
      END
      SUBROUTINE TRA1_SIMPLE2(C,XINT)
*
* Transform one-electron integrals
*
* Input integrals in KINT1
* Output integrals also in KINT1
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
*. Specific Input
      DIMENSION C(*)
*
      NTEST = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA1_S')
*. Largest symmetry block of orbitals 
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Largest number of orbitals in sym block ', MXSOB
      END IF
*. Space for 1e- integrals and two sym blocks
      LENH1 = NTOOB ** 2
      CALL MEMMAN(KLH1,LENH1,'ADDL  ',2,'LH1   ')
      LENSCR = 2 * MXSOB ** 2
      CALL MEMMAN(KLSCR,LENSCR,'ADDL   ',2,'H1SCR ')
      WRITE(6,*) ' LENH1 LENSCR ', LENH1, LENSCR
*. and do it
C     TRAN_SYM_BLOC_MAT3 (AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      WRITE(6,*) ' Before TRAN_SYM '
      CALL TRAN_SYM_BLOC_MAT3(XINT,C,NSMOB,NTOOBS,NTOOBS,
     &                        WORK(KLH1),WORK(KLSCR),1) 
      WRITE(6,*) ' After TRAN_SYM '
*. Copy to original array
      CALL COPVEC(WORK(KLH1),XINT,NINT1)
*. Flush memory
      WRITE(6,*) ' Returning from TRA1 '
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'TRA1_S')
*
      RETURN
      END
c a clone of expand_t1s in tra_simple.f
      SUBROUTINE EXPAND_T1S_NEW(T,T1,ITSM,IC_OCC,IA_OCC,
     &                ISPOBEX_TP,ISOX_TO_OX,IEXC_FOR_OX,
     &                IBSPOBEX_TP,NSPOBEX_TP,IAB)
*
* Slave routine for expanding t1 coefficients to complete matrix
*
* Jeppe Olsen, Korshojen, august 2000
* IAB transferred to argument list, July 2002 to allow for OS
*
* IAB = 1 => alpha excitations are collected 
* IAB = 2 => beta excitations are collected 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cc_exc.inc'
*. General input 
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISOX_TO_OX(NSPOBEX_TP)
      INTEGER IEXC_FOR_OX(*)
      INTEGER IBSPOBEX_TP(NSPOBEX_TP)
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION T1(*)
*. Scratch through input
      INTEGER IC_OCC(NTOOB),IA_OCC(NTOOB)
*. Local scratch
      INTEGER IC_GRP(MXPNGAS),IA_GRP(MXPNGAS),IOFF(8)
*

      IF (ITSM.NE.1) THEN
        WRITE(6,*) 'ITSM.NE.1 in EXPAND_T1S'
        STOP 'EXPAND_T1S'
      END IF

      IDX = 0
      DO ISM = 1, NSMOB
        IOFF(ISM) = IDX
        IDX = IDX + NTOOBS(ISM)*NTOOBS(ISM)        
      END DO
      LEN = IDX

      IF(MSCOMB_CC.EQ.1) THEN
*. Combinations, 
        FACTOR = 1.0D0/SQRT(2.0D0)
      ELSE
*. Individual terms
        FACTOR = 1.0D0
      END IF
*
      ZERO = 0.0D0
      CALL SETVEC(T1,ZERO,LEN)
*
      DO ISOXTP = 1, NSPOBEX_TP
*
      IF(IEXC_FOR_OX(ISOX_TO_OX(ISOXTP)).EQ.1) THEN
*. Is single excitation alpha or beta excitation ?
        NCA_OP = IELSUM(ISPOBEX_TP(1+0*NGAS,ISOXTP),NGAS)
        NCB_OP = IELSUM(ISPOBEX_TP(1+1*NGAS,ISOXTP),NGAS)
        NAA_OP = IELSUM(ISPOBEX_TP(1+2*NGAS,ISOXTP),NGAS)
        NAB_OP = IELSUM(ISPOBEX_TP(1+3*NGAS,ISOXTP),NGAS)
        IELMNT = IBSPOBEX_TP(ISOXTP)
        IF(NCA_OP.EQ.1.AND.IAB.EQ.1.OR.
     &     NCB_OP.EQ.1.AND.IAB.EQ.2     )THEN
          IF(NCA_OP.EQ.1.AND.IAB.EQ.1) THEN
*. Transform from occupations to groups
            CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ISOXTP),IC_GRP,1)
            CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ISOXTP),IA_GRP,1)
            NC_OP = NCA_OP
            NA_OP = NAA_OP
          ELSE IF(NCB_OP.EQ.1.AND.IAB.EQ.2) THEN
*. Transform from occupations to groups
            CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ISOXTP),IC_GRP,1)
            CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ISOXTP),IA_GRP,1)
            NC_OP = NCB_OP
            NA_OP = NAB_OP
          END IF
*
          DO ISM_C = 1, NSMST
            ISM_A = MULTD2H(ISM_C,ITSM)
*. Obtain occupations
            IDUM = 0
            CALL GETSTR2_TOTSM_SPGP(IC_GRP,NGAS,ISM_C,NC_OP,NSTR_C,
     &           IC_OCC, NTOOB,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IA_GRP,NGAS,ISM_A,NA_OP,NSTR_A,
     &           IA_OCC, NTOOB,0,IDUM,IDUM)
*. And scatter to complete matrix 
C?           WRITE(6,*) ' ISM_C, ISM_A, NSTR_A, NSTR_C =', 
C?   &                    ISM_C, ISM_A, NSTR_A, NSTR_C 
             DO ISA = 1, NSTR_A
               DO ISC = 1, NSTR_C
                 IORB_A = IA_OCC(ISA)
                 IORB_C = IC_OCC(ISC)
                 II = IREOTS(IORB_C) - IBSO(ISM_C) + 1
                 JJ = IREOTS(IORB_A) - IBSO(ISM_A) + 1
                 IJDX = IOFF(ISM_A) + (JJ-1)*NTOOBS(ISM_A)+II
                 T1(IJDX) = T(IELMNT)*FACTOR
c and delete original
c                 T(IELMNT) = 0D0
                 IELMNT = IELMNT + 1
               END DO
             END DO
          END DO
        END IF
*       ^ End of correct single excitation
      END IF
*     ^ End if single excitation
      END DO
*     ^ End of loop over excitation types
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T1 in blocked matrix form '
        CALL APRBLM2(T1,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
c $Id$
