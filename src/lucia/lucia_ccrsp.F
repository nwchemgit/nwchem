      SUBROUTINE MV10(CCIN,CCOUT)
*
* Initial outer routine for cc Jacobian times vector
*
* Jeppe Olsen, June 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'lorr.inc'
      INCLUDE 'crun.inc'
C     COMMON/LORR/L_OR_R
*
C      JAC_T_VEC(L_OR_R,CC_AMP,JAC_VEC,TVEC,VEC1,VEC2,
C    &                     CCVEC)
      CALL JAC_T_VEC(L_OR_R,WORK(KCC1),CCOUT,CCIN,
     &     WORK(KVEC1P),WORK(KVEC2P),WORK(KCC5))
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =============================== '
        WRITE(6,*) ' Input and output from JAC_T_VEC ' 
        WRITE(6,*) ' =============================== '
        WRITE(6,*) 
        WRITE(6,*) ' L_OR_R: ', L_OR_R
        CALL WRTMAT(CCIN,1,N_CC_AMP,1,LEN_T_VEC_MX)
        WRITE(6,*)
        CALL WRTMAT(CCOUT,1,N_CC_AMP,1,LEN_T_VEC_MX)
      END IF
*
      RETURN
      END 
      SUBROUTINE CC_EXC_E(CCVEC1,CCVEC2,VEC1,VEC2,
     &                    LUCCAMP,IEXC_RESTRT)
*
* Master routine for coupled cluster linear response calculation 
* of excitation energies
*
* Output : CC-excitation operators are stored on LU_CCEXC_OP
*
* Jeppe Olsen, May 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'lorr.inc'
      INCLUDE 'opti.inc'
*. Scratch vectors ( To show Jesper that I am on the right track )
      DIMENSION VEC1(*),VEC2(*),CCVEC1(*),CCVEC2(*)
*. Local scratch, assuming atmost 1000 roots per sym 
      INTEGER IREO(1000)
      LOGICAL DID_RHS
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CC_EXC ')
*
      WRITE(6,*)
      WRITE(6,*) ' *******************************************'
      WRITE(6,*) ' *                                         *'
      WRITE(6,*) ' * CCLR calculation of excitation energies *'
      WRITE(6,*) ' *                                         *'
      WRITE(6,*) ' *                   Version of May 2000   *'
      WRITE(6,*) ' *                   Revised   June 2004   *'
      WRITE(6,*) ' *                                         *'
      WRITE(6,*) ' *******************************************'
      WRITE(6,*)
      IF(IEXC_RESTRT.NE.0) THEN
        WRITE(6,*) ' Restarted calculation '
      END IF
      IEXC_RESTRT_L = IEXC_RESTRT
*
* Save - if required - excitation vectors from previous run
*
      LU_INIAMP = IOPEN_NUS('CC_INIAMP')
      LU_SCR0   = IOPEN_NUS('CC_EXCSCR0')
      LU_SCR1   = IOPEN_NUS('CC_EXCSCR1')
      LU_SCR2   = IOPEN_NUS('CC_EXCSCR2')

      DID_RHS = .FALSE.
* Loop over right and left equations
      DO ISIDE = 2, 1, -1

*. should this side be done?
        IF (IAND(ICCEX_SLEQ,ISIDE).NE.ISIDE) CYCLE

        WRITE(6,'(X,"+",44("-"),"+")')
        IF (ISIDE.EQ.2)
     &       WRITE(6,*)'| Solving right-hand-side eigenvalue problem |'
        IF (ISIDE.EQ.1)
     &       WRITE(6,*)'| Solving left-hand-side eigenvalue problem  |'
        WRITE(6,'(X,"+",44("-"),"+")')

        IF (ISIDE.EQ.2) THEN
          LU_CCEXC_OP_R = IOPEN_NUS('CC_EXCAMP_R')
          LU_CCEXC_OP = LU_CCEXC_OP_R
          LU_CCEXC_OP_REST = LU_CCEXC_OP_R
        END IF
        IF (ISIDE.EQ.1) THEN
          LU_CCEXC_OP_L = IOPEN_NUS('CC_EXCAMP_L')
          LU_CCEXC_OP = LU_CCEXC_OP_L
          LU_CCEXC_OP_REST = LU_CCEXC_OP_L
          ! restart LHS EVP from RHS solutions
          IF (IEXC_RESTRT.EQ.0.AND.DID_RHS)
     &         LU_CCEXC_OP_REST = LU_CCEXC_OP_R
        END IF

      IF(IEXC_RESTRT.NE.0.OR.DID_RHS) THEN
*. Copy initial sets of excitation operators from file 
*. LU_CCEXC_OP to LU_INIAMP
        IEXC_RESTRT_L = 1

        IF (IEXC_RESTRT.EQ.0.AND.DID_RHS)
     &    WRITE(6,*) ' Initializing with right-hand-side amplitudes ...'

        CALL REWINO(LU_CCEXC_OP_REST)
        CALL REWINO(LU_INIAMP)
*
        DO ISM = 1, NSMST
          IF(NEXC_PER_SYM(ISM).NE.0) THEN 
C          IMSCOMB_CC = 0
           CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ISM,
     &          MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &          WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &          MSCOMB_CC,MX_SBSTR,
     &          WORK(KISPX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &          NTCONF,IPRCC)
C               FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
           NROOT_CC = NEXC_PER_SYM(ISM)
           CALL IFRMDS(LEN_T_VEC_PREV,1,-1,LU_CCEXC_OP_REST)
           DO IROOT_CC = 1, NROOT_CC
             ZERO = 0.0D0
             CALL SETVEC(CCVEC1,ZERO,LEN_T_VEC)
             CALL FRMDSC(CCVEC1,LEN_T_VEC_PREV,-1,LU_CCEXC_OP_REST,
     &                   IMZERO,I_AM_PACKED)
             CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_INIAMP)
           END DO
          END IF 
        END DO
      END IF
*     ^ End of this run is a restart 

      CALL REWINO(LU_CCEXC_OP)
      IF(IEXC_RESTRT_L.NE.0) CALL REWINO(LU_INIAMP)
      DO ISM = 1, NSMST
        IF(NEXC_PER_SYM(ISM).NE.0) THEN 
          WRITE(6,*)
          WRITE(6,*) ' ============================================='
          WRITE(6,*) ' Information for excitations of symmetry',ISM
          WRITE(6,*) ' ============================================='
          WRITE(6,*)
          WRITE(6,*) ' Number of roots required : ', NEXC_PER_SYM(ISM)
          NROOT_CC = NEXC_PER_SYM(ISM)
          ITEX_SM = ISM
*
*. Number of Amplitudes for this symmetry
*
C       IMSCOMB_CC = 0
           CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ISM,
     &          MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &          WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &          MSCOMB_CC,MX_SBSTR,
     &          WORK(KISPX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &          NTCONF,IPRCC)
        WRITE(6,*) ' Number of CC amplitudes ', LEN_T_VEC
        IF (LEN_T_VEC.GT.N_CC_AMP) STOP 'DIMENSIONS'
*
*. Set up Diagonal for this symmetry
*
        IMOD = 2 ! use alpha/beta Fock-matrix
        CALL GENCC_F_DIAG_M(IMOD,WORK(KLSOBEX),NSPOBEX_TP,CCVEC1,ISM,
     &                      XDUM,IDUM,IDUM,0,
     &                      VEC1,VEC2,MX_ST_TSOSO,MX_ST_TSOSO_BLK)
*. Save on LUDIA for future use
        CALL REWINO(LUDIA)
        CALL TODSC(CCVEC1,LEN_T_VEC,-1,LUDIA)
*
*. Initialization : Restart or lowest diagonal elements 
*
*
        IF(IEXC_RESTRT_L.EQ.0) THEN
*. Lowest elements of diagonal
          CALL DUMSORT(CCVEC1,LEN_T_VEC,NROOT_CC,IREO,CCVEC2)
*. Save corresponding initial guesses
          CALL REWINO(LU_SCR0)
          DO IROOT = 1, NROOT_CC
            ZERO = 0.0D0
            CALL SETVEC(CCVEC1,ZERO,LEN_T_VEC)
            CCVEC1(IREO(IROOT)) = 1.0D0
            CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_SCR0)
C?          WRITE(6,*) ' Initial vector for IROOT = ', IROOT
C?          CALL WRTMAT(CCVEC1,1,LEN_T_VEC,1,LEN_T_VEC) 
          END DO
        ELSE IF(IEXC_RESTRT_L.NE.0) THEN
          CALL REWINO(LU_SCR0)
*. Read initial vectors from LU_INIAMP
          DO IROOT =1, NROOT_CC
           CALL FRMDSC(CCVEC1,LEN_T_VEC,-1,LU_INIAMP,
     &                 IMZERO,I_AM_PACKED)
           CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_SCR0)
          END DO
        END IF
*       ^ End of shift between different forms of initialization
*. Allocate scratch space for the diagonalization 
*. Length of APROJ : MAXVEC**2
*. Length of AVEC  : MAXVEC**2
*  Length of WORK  : 4*MAXVEC**2 + 7*MAXVEC
*  Length of H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
        MAXVEC = MXCIV*NROOT_CC
        CALL MEMMAN(KLAPROJ,MAXVEC**2,'ADDL  ',2,'APROJ ')
        CALL MEMMAN(KLAVEC ,MAXVEC**2,'ADDL  ',2,'AVEC  ')
        LWORK = 4*MAXVEC**2 + 7*MAXVEC
        CALL MEMMAN(KLWORK,LWORK,'ADDL  ',2,'LWORK ')
        CALL MEMMAN(KLSCR ,0,'ADDL  ',2,'LSCR  ')
        LEN = NROOT_CC*MAXIT
        CALL MEMMAN(KLRNRM,LEN,'ADDL  ',2,'RNRM  ')
        CALL MEMMAN(KLEIG ,LEN,'ADDL  ',2,'EIG   ')
        CALL MEMMAN(KLFINEIG,NROOT_CC,'ADDL  ',2,'FINEIG')
*. And then to the work
C     GENMAT_DIAG(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
C    &            LU3,LUDIA,NROOT,MAXVEC,NINVEC,
C    &            APROJ,AVEC,WORK,IPRT,
C    &            NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
C    &            IOLSEN,IPICO)
*. We will obtain right eigenvectors. Tell program 
*. to perform Jacobian times right vectors
        L_OR_R = ISIDE
        IPRDIA_L = 10
        CONV = CCEX_CONV

        CALL GENMAT_DIAG(L_OR_R,
     &       CCVEC1,CCVEC2,VEC1,VEC2,
     &       LU_SCR0,LU_SCR1,LUCCAMP,
     &       WORK(KLRNRM),WORK(KLEIG),WORK(KLFINEIG),MAXIT,LEN_T_VEC,
     &                                                     N_CC_AMP,
     &       LU_SCR2,LUDIA,NROOT_CC,MAXVEC,NROOT_CC,
     &       WORK(KLAPROJ),WORK(KLAVEC),WORK(KLWORK),IPRDIA_L,
     &       0,0.0D0,0,0,0,0,0.0D0,0.0D0,
     &       0,0,CONV)
*
* Analyze excitation vectors - and copy to LU_CCEXC
*
        CALL REWINO(LU_SCR0)
        CALL ITODS(LEN_T_VEC,1,-1,LU_CCEXC_OP)
        DO IROOT = 1, NROOT_CC
          WRITE(6,'(/,X,"+",75("="),"+")')
          WRITE(6,'(X,"|",2X,A,I2,A,I4,39X,"|")')
     &         'Analysis for symmetry ',ISM,' root ',IROOT
          WRITE(6,'(X,"|",2X,A,F20.12,34X,"|")')
     &         'Excitation energy: ',WORK(KLFINEIG+IROOT-1)
          WRITE(6,'(X,"|",2X,A,F20.12,34X,"|")')
     &         'Residual norm:     ',WORK(KLRNRM+IROOT-1)
          WRITE(6,'(X,"+",75("="),"+")')
          CALL FRMDSC(CCVEC1,LEN_T_VEC,-1,LU_SCR0,IMZERO,IAMPACK)
          CALL ANA_GENCC(CCVEC1,ISM)
          CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_CCEXC_OP)
        END DO
*
       END IF
      END DO
*
        IF (ISIDE.EQ.2) DID_RHS = .TRUE.

      END DO ! ISIDE
*
      IF (IAND(1,ICCEX_SLEQ).EQ.1) CALL RELUNIT(LU_CCEXC_OP_L,'KEEP')
      IF (IAND(2,ICCEX_SLEQ).EQ.2) CALL RELUNIT(LU_CCEXC_OP_R,'KEEP')
      CALL RELUNIT(LU_INIAMP,'DELETE')
      CALL RELUNIT(LU_SCR0,'DELETE')
      CALL RELUNIT(LU_SCR1,'DELETE')
      CALL RELUNIT(LU_SCR2,'DELETE')
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CC_EXC ')
      RETURN
      END

      SUBROUTINE GET_NONSYM_SUBMAT(LUC,LUHC,NVEC,NVECP,SUBMAT,NDIM,
     &                             VEC1,VEC2,SUBMATP)
*
* Obtain subspace matrix from trialvectors and matrix times 
* trial vectors. No assumption about symmetry of matrix
*
* Jeppe Olsen, May 2000
*
* Version assuming the ability to hold two vectors in core
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input : Previous subspace matrix
      DIMENSION SUBMATP(*)
*. Output : New subspace matrix 
      DIMENSION SUBMAT(*)
*. Scratch 
      DIMENSION VEC1(*),VEC2(*)
*
      LBLK = -1
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Input vectors to GET_NONSYM_SUBMAT '
        WRITE(6,*) ' ==================================='
*
        WRITE(6,*) ' C vectors '
        CALL REWINO(LUC)
        DO IVEC = 1, NVEC
          CALL FRMDSC(VEC1,NDIM,LBLK,LUC,IMZERO,IAMPACK)
          WRITE(6,*) ' Input C vector number = ', IVEC
          CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END DO
*
        WRITE(6,*) ' HC vectors '
        CALL REWINO(LUHC)
        DO IVEC = 1, NVEC
          CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
          WRITE(6,*) ' Input HC vector number = ', IVEC
          CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END DO
*
      END IF
*
      IMZERO = 0
      IAMPACK = 0
*.Reform previous matrix from (NVECP,NVECP) to (NVEC,NVEC) format
      IF(NVECP.GT.0) THEN
C            COPMAT(AIN,AOUT,NIN,NOUT)
        CALL COPMAT(SUBMATP,SUBMAT,NVECP,NVEC)
      END IF
*
*. Add new columns of subspace matrix
*
C?    WRITE(6,*) ' NVEC, NVECP = ', NVEC,NVECP
      CALL REWINO(LUHC)
      DO IVEC = 1, NVECP
        CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
      END DO
      DO J = 1, NVEC-NVECP
*. Read in Hc(nvecp+j)
C?      WRITE(6,*) ' J = ', J
        CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
        CALL REWINO(LUC)
        DO I = 1, NVEC
C?        WRITE(6,*) ' I = ', I
*. Read in c(i)
          CALL FRMDSC(VEC2,NDIM,LBLK,LUC,IMZERO,IAMPACK)
* Submat(i,nvecp+j) = c(i)t Hc(nvecp+j)
          SUBMAT((NVECP+J-1)*NVEC+I) = INPROD(VEC2,VEC1,NDIM)
        END DO
      END DO
*
*. Add new rows of subspace matrix 
*
C?    WRITE(6,*) ' Part 2 ' 
      CALL REWINO(LUC)
      DO I = 1, NVECP
        CALL FRMDSC(VEC2,NDIM,LBLK,LUC,IMZERO,IAMPACK)
      END DO
C     CALL SKPRC3(NVECP,LUC)
      DO I = 1, NVEC-NVECP
C?      WRITE(6,*) ' I = ', I
*. Read in c(nvecp+i)
        CALL FRMDSC(VEC2,NDIM,LBLK,LUC,IMZERO,IAMPACK)
        CALL REWINO(LUHC)
        DO J = 1, NVECP
C?        WRITE(6,*) ' J = ', J
*. Read in Hc(j)
          CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
* Submat(i,nvecp+j) = c(i)t Hc(nvecp+j)
          SUBMAT((J-1)*NVEC+I+NVECP) = INPROD(VEC2,VEC1,NDIM)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated subspace matrix '
        CALL WRTMAT(SUBMAT,NVEC,NVEC,NVEC,NVEC)
      END IF
*
      RETURN
      END
      SUBROUTINE GENMAT_DIAG(L_OR_R,
     &                  VEC1,VEC2,SCRVEC1,SCRVEC2,
     &                  LU1,LU2,LUCCAMP,!LUINTM1,LUINTM2,LUINTM3,
     &                  RNRM,EIG,FINEIG,MAXIT,NVAR,NVAR0,
     &                  LU3,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &                  IOLSEN,IPICO,TEST)
*
* p.t. matrix vector routine is hardwired 
*
* Iterative solution of eigenvalue problem for general matrix
*
* Final and intermediate eigenvalues and eigenvectors are assumed real
*
* MIN version requiring two vectors in core
*
*
* Jeppe Olsen May 2000, from MINDV4
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors, each must be dimensioned to hold
*                    complete vector
*        LU2,LU3   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
* On input LU1 is supposed to hold initial guess to eigenvectors
*
* IOLSEN : Use inverse iteration modified Davidson
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       LOGICAL CONVER
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION FINEIG(1)
       DIMENSION H0(*),IPNTR(*)
*
*. Scratch through argument list
*
*. Length of APROJ : MAXVEC**2
*. Length of AVEC  : MAXVEC**2
*  Length of WORK  : 4*MAXVEC**2 + 7*MAXVEC
*  Length of H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0SCR(*)
*. Local scratch - atmost 1000 roots
       LOGICAL RTCNV(1000)
*. Ad hoc arrays for root hooming
       DIMENSION XJEP(3000),IXJEP(3000)
*
*
       CALL ATIM(CPU0,WALL0)
*
       IOLSTM = 0      
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Variational'
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Perturbational'
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Normal Davidson, Variational '
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Normal Davidson, Perturbational'
C      IF( MAXVEC .LT. 2 * NROOT ) THEN
C        WRITE(6,*) ' SORRY MINDV2 WOUNDED , MAXVEC .LT. 2*NROOT '
C        STOP ' ENFORCED STOP IN MINDV2'
C      END IF
       WRITE(6,*) ' Convergence threshold for residual : ', TEST
*
       I_DO_ROOTHOOMING = 0
       IF(I_DO_ROOTHOOMING.EQ.1) THEN
         WRITE(6,*) ' Root hooming active '
       END IF
*
*. Scratch files
       LUINTM1 = IOPEN_NUS('CCGMATSCR1')
       LUINTM2 = IOPEN_NUS('CCGMATSCR2')
       LUINTM3 = IOPEN_NUS('CCGMATSCR3')
*
       IF(IPICO.NE.0) THEN
         MAXVEC = 2*NROOT
       END IF
*. Division of scratch memory 
      KAPROJ = 1
      KFREE = KAPROJ + MAXVEC**2
*
      KARVAL = KFREE
      KFREE = KARVAL + MAXVEC
*
      KAIVAL = KFREE
      KFREE = KAIVAL + MAXVEC
*
      KARVEC = KFREE
      KFREE = KARVEC + MAXVEC**2
*
      KAIVEC = KFREE
      KFREE = KAIVEC + MAXVEC**2
*
      KZ = KFREE
      KFREE = KZ    + MAXVEC**2
*
      KW = KFREE
      KFREE = KW    + MAXVEC
*
      KSCR1 = KFREE
      KFREE = KSCR1 + 4*MAXVEC
*
C     TEST = 1.0D-10
      CONVER = .FALSE.
      DO 1234 MACRO = 1,1
*
*.   INITAL ITERATION
*
* ===========================================================
*. The initial vectors does not neccessarily constitute an 
*. orthonormal basis. Start by orthogonalizing
* ===========================================================
*
       DO IROOT = 1, NINVEC
*. Read vector IROOT in
         CALL REWINO(LU1)
         DO ISKP = 1, IROOT
           CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
         END DO
*. Diagonal element
         WORK(KAPROJ-1+(IROOT-1)*NINVEC+IROOT) = 
     &   INPROD(VEC1,VEC1,NVAR)
         DO JROOT = IROOT+1, NINVEC
           CALL FRMDSC(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
           WORK(KAPROJ-1+(IROOT-1)*NINVEC+JROOT) = 
     &     INPROD(VEC1,VEC2,NVAR)
           WORK(KAPROJ-1+(JROOT-1)*NINVEC+IROOT) = 
     &     WORK(KAPROJ-1+(IROOT-1)*NINVEC+JROOT)   
         END DO
       END DO
*
       IF(IPRT.GE.1) THEN
        WRITE(6,*) ' Overlap matrix of initial basis '
        CALL WRTMAT(WORK(KAPROJ),NROOT,NROOT,NROOT,NROOT)
       END IF
*. Orthonormal basis for the NROOT lowest roots
       CALL MGS3(WORK(KARVEC),WORK(KAPROJ),NINVEC,WORK(KAIVEC))
*. Orthogonalize, save orthogonlized vectors on LU3
       CALL REWINO( LU3)
       DO IROOT = 1, NINVEC
         CALL REWINO( LU1)
         CALL SETVEC(VEC1,0.0D0,NVAR)
         DO IVEC = 1, NINVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR =  WORK(KARVEC-1+(IROOT-1)*NINVEC+IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
         END DO
         CALL TODSC(VEC1,NVAR,-1  ,LU3)
       END DO
*
       CALL REWINO( LU1)
       CALL REWINO( LU3)
       DO IVEC = 1,NINVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
         CALL TODSC (VEC1,NVAR,-1,  LU1)
       END DO
*
       ITER = 1
       CALL REWINO( LU1 )
       CALL REWINO( LU2 )
* Mat * initial  vectors 
       DO JVEC = 1,NINVEC
         IADD_RHS = 0
         IF (JVEC.EQ.1) THEN
           LUINTM1_L = -LUINTM1
           LUINTM2_L = -LUINTM2
         ELSE
           LUINTM1_L = LUINTM1
           LUINTM2_L = LUINTM2
         END IF
         IWRMOD=1
         XDUM = -1
         LUDUM = 0
         CALL JAC_T_VEC2(L_OR_R,0,IWRMOD,0,0,
     &        VEC1,VEC2,SCRVEC1,SCRVEC2,
     &        NVAR,NVAR0,
     &        XDUM,XNORM,RNORM,
     &        LUCCAMP,LUDUM,LU1,LU2,LUDUM,
     &        LUINTM1_L,LUINTM2_L,LUINTM3)
c         CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
c         CALL MV10(VEC1,VEC2)
c         CALL TODSC(VEC2,NVAR,-1  ,LU2)
       END DO 
*. Projected matrix
       CALL GET_NONSYM_SUBMAT(LU1,LU2,NINVEC,0,APROJ,NVAR,
     &                        VEC1,VEC2,0.0D0)
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL WRTMAT(APROJ,NINVEC,NINVEC,NINVEC,NINVEC)
       END IF
*  Diagonalize initial subspace matrix  
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*NINVEC)
       CALL EIGGMT3(WORK(KAPROJ),NINVEC,WORK(KARVAL),WORK(KAIVAL),
     &              WORK(KARVEC),WORK(KAIVEC),WORK(KZ),WORK(KW),
     &              WORK(KSCR1),1,1)
C          EIGGMT3(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
C    &                   Z,W,SCR,IORD,IEIGVC)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KARVAL-1+IROOT)
   20  CONTINUE
       CALL COPVEC(WORK(KARVEC),AVEC,NINVEC**2)
*
       IF( IPRT  .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Initial set of eigenvalues '
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF),IROOT=1,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
         WRITE(6,'(">>>",A)')
     &        '  Iter. Root     exc. energy     res. norm  conv.'
         WRITE(6,'(">>>",A)')
     &        '-------------------------------------------------'
*
** LOOP OVER ITERATIONS
*
 1000 CONTINUE
      IF(IPRT  .GE. 10 ) THEN
       WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
      END IF
 
 
        ITER = ITER + 1
*
** 1          NEW DIRECTION TO BE INCLUDED
*
*   1.1 : R = H*X - EIGAPR*X
       IADD = 0
       CONVER = .TRUE.
       RNRMMX = 0d0
       DO 100 IROOT = 1, NROOT
         CALL SETVEC(VEC1,0.0D0,NVAR)
*
         CALL REWINO( LU2)
         DO 60 IVEC = 1, NVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   60    CONTINUE
*
         EIGAPR = EIG(ITER-1,IROOT)
         CALL REWINO( LU1)
         DO 50 IVEC = 1, NVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR = (-EIGAPR)*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   50    CONTINUE
           IF ( IPRT  .GE.600 ) THEN
             WRITE(6,*) '  ( HX - EX ) '
             CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
         RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
         RNRM(ITER-1,IROOT) = RNORM
         RNRMMX = MAX(RNRMMX,RNORM)
         IF(RNORM.LT. TEST ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
*
         WRITE(6,'(">>>",5X,I5,F20.12,2X,E10.4,4X,L)')
     &        IROOT,EIGAPR,RNORM,RTCNV(IROOT)
*
         IF( ITER .GT. MAXIT) GOTO 100
*.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
           CALL REWINO( LUDIA)
           CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL H0M1TV(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
              CALL REWINO(LU3)
              CALL TODSC(VEC1,NVAR,-1,LU3)
* Current eigen vector
              CALL REWINO( LU1)
              CALL SETVEC(VEC1,0.0D0,NVAR)
              DO 59 IVEC = 1, NVEC
                CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
                FACTOR = AVEC((IROOT-1)*NVEC+ IVEC)
                CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   59         CONTINUE
              IF ( IPRT  .GE. 600 ) THEN
                WRITE(6,*) ' And X  '
                CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
              END IF
              CALL TODSC(VEC1,NVAR,-1,LU3)
* (H0 - E )-1  * X
              CALL REWINO( LUDIA)
              CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
              CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
     &                   H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
              CALL TODSC(VEC2,NVAR,-1,LU3)
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPROD(VEC2,VEC1,NVAR)
* is X an eigen vector for (H0 - 1 ) - 1
              CALL VECSUM(VEC2,VEC1,VEC2,GAMMA,-1.0D0,NVAR)
              VNORM = SQRT(MAX(0.0D0,INPROD(VEC2,VEC2,NVAR)))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                CALL REWINO(LU3)
                CALL FRMDSC(VEC2,NVAR,-1,LU3,IMZERO,IAMPACK)
                DELTA = INPROD(VEC1,VEC2,NVAR)
                CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                FACTOR = (-DELTA)/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSUM(VEC1,VEC1,VEC2,FACTOR,1.0D0,NVAR)
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) '  Modified new trial vector '
                  CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
                END IF
              ELSE
                IF(IPRT.GT.0) WRITE(6,*) 
     &          ' Inverse correction switched of'
                CALL REWINO(LU3)
                CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
              END IF
            END IF
*. 1.3 ORTHOGONALIZE TO ALL PREVIOUS VECTORS
           XNRMI =    INPROD(VEC1,VEC1,NVAR)
           CALL REWINO( LU1 )
 
           DO 80 IVEC = 1,NVEC+IADD-1
             CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
   80      CONTINUE
*. 1.4 Normalize vector and check for linear dependency
           SCALE = INPROD(VEC1,VEC1,NVAR)
           IF(ABS(SCALE)/XNRMI .LT. 1.0D-10) THEN
*. Linear dependency
             IADD = IADD - 1
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(6,*) '  Trial vector linear dependent so OUT !!! '
             END IF
           ELSE
             C1NRM = SQRT(SCALE)
             FACTOR = 1.0D0/SQRT(SCALE)
             CALL SCALVE(VEC1,FACTOR,NVAR)
*
             CALL TODSC(VEC1,NVAR,-1  ,LU1)
             IF ( IPRT  .GE.600 ) THEN
               WRITE(6,*) 'ORTHONORMALIZED (D-E)-1 *( HX - EX ) '
               CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
             END IF
           END IF
*
         END IF
  100 CONTINUE
      WRITE(6,'(">>>",I5,25X,(2X,E10.4))')
     &        ITER-1,RNRMMX
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
*
*  2.1: MULTIPLY NEW DIRECTION WITH MATRIX
       CALL REWINO( LU1)
       CALL REWINO( LU2)
       DO 110 IVEC = 1, NVEC
         CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
         CALL FRMDSC(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
  110  CONTINUE
*
      DO 150 IVEC = 1, IADD
         IWRMOD = 1
         XDUM = -1
         LUDUM = 0
         CALL JAC_T_VEC2(L_OR_R,0,IWRMOD,0,0,
     &        VEC1,VEC2,SCRVEC1,SCRVEC2,
     &        NVAR,NVAR0,
     &        XDUM,XNORM,RNORM,
     &        LUCCAMP,LUDUM,LU1,LU2,LUDUM,
     &        LUINTM1,LUINTM2,LUINTM3)
c        CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
c        CALL MV10(VEC1,VEC2)
c        CALL TODSC(VEC2,NVAR,-1  ,LU2)
  150 CONTINUE
*.Augment projected matrix
C     GET_NONSYM_SUBMAT(LUC,LUHC,NVEC,NVECP,SUBMAT,NDIM,
C    &                             VEC1,VEC2,SUBMATP)
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC**2)
      CALL GET_NONSYM_SUBMAT(LU1,LU2,NVEC+IADD,NVEC,APROJ,NVAR,
     &                       VEC1,VEC2,WORK(KAPROJ))
*  DIAGONALIZE PROJECTED MATRIX
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*NVEC)
      CALL EIGGMT3(WORK(KAPROJ),NVEC,WORK(KARVAL),WORK(KAIVAL),
     &             WORK(KARVEC),WORK(KAIVEC),WORK(KZ),WORK(KW),
     &             WORK(KSCR1),1,1)
      
      DO  IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KARVAL-1+IROOT)
      END DO
      CALL COPVEC(WORK(KARVEC),AVEC,NROOT*NVEC)
      IF(I_DO_ROOTHOOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first 
*
*. Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
       END DO
       WRITE(6,*) 
     & ' Norm of projections to previous vector space '
       CALL WRTMAT(XJEP,1,NVEC,1,NVEC)
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL SCALVE(XJEP,ONEM,NVEC)
       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
       IF(NSORT.LT.NVEC) THEN
         WRITE(6,*) ' Warning : Some elements lost in sorting '
         WRITE(6,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
       IF(IPRT.GE.0) THEN
         WRITE(6,*) ' New roots choosen as vectors '
         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
       END IF
*. Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
       END DO
       CALL COPVEC(XJEP,AVEC,NVEC*NVEC)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW) = WORK(KARVAL-1+IOLD)               
       END DO
       DO INEW = 1, NROOT
         EIG(ITER,INEW)=XJEP(INEW)
       END DO
*
       IF(IPRT.GE.3) THEN
         WRITE(6,*) ' Reordered AVEC arrays '
         CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
       END IF
*
      END IF
*     ^ End of root homing procedure
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
           WRITE(6,*) ' Norms of residuals '
           WRITE(6,'(10E13.5)') (RNRM(ITER-1,KROOT),KROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' PROJECTED MATRIX AND EIGEN PAIRS '
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
      END IF
*
**  PERHAPS RESET OR ASSEMBLE CONVERGED EIGENVECTORS
*
  901 CONTINUE
*
      IPULAY = 0
      IF(IPULAY.EQ.1 .AND. MAXVEC.EQ.3 .AND.NVEC.GE.2.
     &   .AND. .NOT.CONVER) THEN
* Save trial vectors : 1 -- current trial vector
*                      2 -- previous trial vector orthogonalized
        CALL REWINO( LU3)
        CALL REWINO( LU1)
*. Current trial vector
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2200 IVEC = 1, NVEC
          CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
         CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2200   CONTINUE
        SCALE = INPROD(VEC1,VEC1,NVAR)
        SCALE  = 1.0D0/SQRT(SCALE)
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC(VEC1,NVAR,-1  ,LU3)
* Previous trial vector orthonormalized
        CALL REWINO(LU1)
        CALL FRMDSC(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        SCALE2 = INPROD(VEC2,VEC2,NVAR)
        SCALE2 = 1.0D0/SQRT(SCALE2)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC(VEC2,NVAR,-1,LU3)
*
        CALL REWINO( LU1)
        CALL REWINO( LU3)
        DO 2411 IVEC = 1,2
          CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC1,NVAR,-1,  LU1)
 2411   CONTINUE
*. Corresponding sigma vectors
        CALL REWINO ( LU3)
        CALL REWINO( LU2)
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2250 IVEC = 1, NVEC
          CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
          CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2250   CONTINUE
*
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC(VEC1,NVAR,-1,  LU3)
* Sigma vector corresponding to second vector on LU1
        CALL REWINO(LU2)
        CALL FRMDSC(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC(VEC2,NVAR,-1,LU3)
*
        CALL REWINO( LU2)
        CALL REWINO( LU3)
        DO 2400 IVEC = 1,2
          CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC2,NVAR,-1  ,LU2)
 2400   CONTINUE
        NVEC = 2
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 2410 IROOT = 1,NVEC
          AVEC((IROOT-1)*NVEC+IROOT) = 1.0D0
 2410   CONTINUE
*.Projected hamiltonian
       CALL REWINO( LU1 )
       DO 2010 IVEC = 1,NVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL REWINO( LU2)
         DO 2008 JVEC = 1, IVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
 2008    CONTINUE
 2010  CONTINUE
      END IF
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
*. Select space spanning lowest NROOT as new subspace    
*. Note that subspace is required to be orthogonal although 
*. the eigenvectors in general not are orthogonal
*
*. Overlap matrix of lowest NROOT eigenvectors
        DO IROOT = 1, NROOT
          DO JROOT = 1, NROOT
            WORK(KAPROJ-1+(IROOT-1)*NROOT+JROOT) = 
     &      INPROD(AVEC((IROOT-1)*NVEC+1),
     &             AVEC((JROOT-1)*NVEC+1),NVEC)
          END DO
        END DO
*
        IF(IPRT.GE.1) THEN
         WRITE(6,*) ' Overlap matrix in new basis '
         CALL WRTMAT(WORK(KAPROJ),NROOT,NROOT,NROOT,NROOT)
        END IF
*. Orthonormal basis for the NROOT lowest roots
C       MGS3(X,S,NDIM,SCR1)
        CALL MGS3(WORK(KARVEC),WORK(KAPROJ),NROOT,WORK(KAIVEC))
*. Orthogonalization matrix is now in WORK(KARVEC)
*  Obtain new basis  - if iteration procedure continues 
        IF(ITER.LT.MAXIT .AND. (.NOT.CONVER)) THEN
          CALL MATML4(WORK(KAIVEC),AVEC,WORK(KARVEC),
     &    NVEC,NROOT,NVEC,NROOT,NROOT,NROOT,0)
          CALL COPVEC(WORK(KAIVEC),AVEC,NROOT*NVEC)
        END IF
*
        CALL REWINO( LU3)
        DO 320 IROOT = 1, NROOT
          CALL REWINO( LU1)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 200 IVEC = 1, NVEC
            CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  200     CONTINUE
*
          SCALE = INPROD(VEC1,VEC1,NVAR)
          SCALE  = 1.0D0/SQRT(SCALE)
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC(VEC1,NVAR,-1  ,LU3)
  320   CONTINUE
        CALL REWINO( LU1)
        CALL REWINO( LU3)
        DO 411 IVEC = 1,NROOT
          CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC1,NVAR,-1,  LU1)
  411   CONTINUE
* CORRESPONDING SIGMA VECTOR
        CALL REWINO ( LU3)
        DO 329 IROOT = 1, NROOT
          CALL REWINO( LU2)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 250 IVEC = 1, NVEC
            CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  250     CONTINUE
*
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC(VEC1,NVAR,-1,  LU3)
  329   CONTINUE
* PLACE C IN LU1 AND HC IN LU2
        CALL REWINO( LU2)
        CALL REWINO( LU3)
        DO 400 IVEC = 1,NROOT
          CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC2,NVAR,-1  ,LU2)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL GET_NONSYM_SUBMAT(LU1,LU2,NROOT,0,APROJ,NVAR,
     &                       VEC1,VEC2,WORK(KAPROJ))
C
      END IF
C
C     IF( ITER .LT. MAXIT .AND. .NOT. CONVER) GOTO 1000
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*. Place first eigenvector in vec1
      CALL REWINO(LU1)
      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
 
 
 
* ( End of loop over iterations )
*
*
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*. Final eigenvalues
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
*
      IF ( IPRT .GT. 1 ) THEN
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL REWINO( LU1)
            CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),
     &                 RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
 1234 CONTINUE
C
*
*. Return final residual norm of roots:
      DO IROOT = 1, NROOT
        RNRM(IROOT,1) = RNRM(ITER,IROOT)
      END DO
*
      CALL RELUNIT(LUINTM1,'DELETE')
      CALL RELUNIT(LUINTM2,'DELETE')
      CALL RELUNIT(LUINTM3,'DELETE')

      CALL ATIM(CPU,WALL)
      CALL PRTIM(6,'time in GEVP solver',CPU-CPU0,WALL-WALL0)
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE DUMSORT(VEC,NDIM,NELMNT,IREO,ISCR)      
*
* Extremely stupid routine for finding the ordering of 
* elements in VEC. Only the lowest NELMNT elements are obtained
*
* On output : IREO: ordered => unordered index
*
* Author : Prefers to be anonymous, May 2000
*
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(NDIM)
*. Output 
      INTEGER IREO(NELMNT)                
*. Scratch through parameter list
      INTEGER ISCR(NDIM)
*
      IZERO = 0
      CALL ISETVC(ISCR,IZERO,NDIM)   
      XMAX = FNDMNX(VEC,NDIM,2)
C?    WRITE(6,*) ' XMAX = ', XMAX
*
      DO I = 1, NELMNT
*. Find the I'th lowest element 
        IELMNT = 0
        XVAL = XMAX
        DO J = 1, NDIM
C?       WRITE(6,*) ' I,J,XVAL,VEC(J) = ', I,J,XVAL,VEC(J)
         IF(VEC(J).LE.XVAL.AND.ISCR(J).EQ.0) THEN
C?         WRITE(6,*) ' Update of lowest element I and J = ', I,J
           XVAL = VEC(J)
           IELMNT = J
         END IF
        END DO
        IREO(I) = IELMNT
        ISCR(IELMNT) = I
      END DO
*
      NTEST = 10
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Elements to be sorted : '
        CALL WRTMAT(VEC,1,NDIM,1,NDIM)
      END  IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Reorder array : new => old order '
        CALL IWRTMA(IREO,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE REO_VEC(IREO,VECIN,NDIM,VECOUT,IWAY)
*
* Reorder vector
*
* IWAY = 1 : VECOUT(I) = VECIN(IREO(I))
* IWAY = 2 : VECOUT(IREO(I)) = VECIN(I)
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VECIN(NDIM)
      INTEGER IREO(NDIM)
*. Output
      DIMENSION VECOUT(NDIM)
*
      IF(IWAY.EQ.1) THEN
        DO I = 1, NDIM
          VECOUT(I) = VECIN(IREO(I))
        END DO
      ELSE
        DO I = 1, NDIM
          VECOUT(IREO(I)) = VECIN(I)
        END DO 
      END IF
*
      RETURN
      END
      SUBROUTINE REO_COL_MAT(IREO,AIN,AOUT,NR,NC,IWAY)
*
* Reorder columns of matrix
*
* IWAY = 1 : AOUT(I,J) = AIN(I,IREO(J))
* IWAY = 2 : AOUT(I,IREO(J)) = AIN(I,J)
*
* Jeppe Olsen, May of 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION AIN(NR,NC) 
      INTEGER IREO(*)
*. Output
      DIMENSION AOUT(NR,NC)
*
      IF(IWAY.EQ.1) THEN
        DO J = 1, NC
          CALL COPVEC(AIN(1,IREO(J)),AOUT(1,J),NR)
        END DO
      ELSE
        DO J = 1, NC
          CALL COPVEC(AIN(1,J),AOUT(1,IREO(J)),NR)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE EIGGMT3(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
     &                   Z,W,SCR,IORD,IEIGVC)
*
* IF IEIGVC = 1 : Calculate eigenvalues and eigenvalues 
*           = 0 : Calculate only eigenvalues
*
* Outer routine for calculating eigenvectors and eigenvalues
* of a general real matrix
*
* Version employing EISPACK path RG
*
* Current implementation is rather wastefull with respect to
* memory but at allows one to work with real arithmetic
* outside this routine
*
* If IORD.EQ.1, the eigenvalues are oredered according to the 
* size of the real part of the eigenvalues
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8 INPROD
      DIMENSION AMAT(NDIM,NDIM),SCR(2*NDIM)
      DIMENSION ARVAL(NDIM),AIVAL(NDIM)
      DIMENSION ARVEC(NDIM,NDIM),AIVEC(NDIM,NDIM)
      DIMENSION Z(NDIM,NDIM),W(NDIM)
*
* Diagonalize
*
      NSCR = 2*NDIM
      IF(IEIGVC.EQ.1) THEN
*. Eigenvalues and eigenvectors 
       CALL RG(NDIM,NDIM,AMAT,ARVAL,AIVAL,1,Z,SCR(1),SCR(1+NDIM),IERR)
      ELSE 
* only eigenvalues :  
       CALL RG(NDIM,NDIM,AMAT,ARVAL,AIVAL,0,DUMMY,SCR(1),SCR(1+NDIM),
     &         IERR)
      END IF
      IF( IERR.NE.0) THEN
        WRITE(6,*) ' Problem in EIGGMTN, no convergence '
        WRITE(6,*) ' I have to stop '
        STOP ' No convergence in EIGGMTN '
      END IF
*
* Extract real and imaginary parts according to Eispack manual p.89
*
      IF(IEIGVC.EQ.1) THEN
      DO 150 K = 1, NDIM
*
        IF(AIVAL(K).NE.0.0D0) GOTO 110
        CALL COPVEC(Z(1,K),ARVEC(1,K),NDIM)
        CALL SETVEC(AIVEC(1,K),0.0D0,NDIM)
        GOTO 150
*
  110   CONTINUE
        IF(AIVAL(K).LT.0.0D0) GOTO 130
        CALL COPVEC(Z(1,K),ARVEC(1,K),NDIM)
        CALL COPVEC(Z(1,K+1),AIVEC(1,K),NDIM)
        GOTO 150
*
  130   CONTINUE
        CALL COPVEC(ARVEC(1,K-1),ARVEC(1,K),NDIM)
        CALL VECSUM(AIVEC(1,K),AIVEC(1,K),AIVEC(1,K-1),
     &              0.0D0,-1.0D0,NDIM)
*
  150 CONTINUE
 
*
* explicit orthogonalization of eigenvectors with
* (degenerate eigenvalues are not orthogonalized by DGEEV)
*
      GOTO 201
      TEST = 1.0D-11
      DO 200 IVEC = 1, NDIM
         RNORM = INPROD(ARVEC(1,IVEC),ARVEC(1,IVEC),NDIM)
     &         + INPROD(AIVEC(1,IVEC),AIVEC(1,IVEC),NDIM)
         FACTOR = 1.0d0/SQRT(RNORM)
         CALL SCALVE(ARVEC(1,IVEC),FACTOR,NDIM)
         CALL SCALVE(AIVEC(1,IVEC),FACTOR,NDIM)
         DO 190 JVEC = IVEC+1,NDIM
           IF(ARVAL(IVEC)-ARVAL(JVEC).LE.TEST) THEN
* orthogonalize jvec to ivec
           OVLAPR = INPROD(ARVEC(1,IVEC),ARVEC(1,JVEC),NDIM)
     &            + INPROD(AIVEC(1,JVEC),AIVEC(1,IVEC),NDIM)
           OVLAPI = INPROD(ARVEC(1,IVEC),AIVEC(1,JVEC),NDIM)
     &            - INPROD(AIVEC(1,IVEC),ARVEC(1,JVEC),NDIM)
           CALL VECSUM(ARVEC(1,JVEC),ARVEC(1,JVEC),ARVEC(1,IVEC),
     &                 1.0D0,-OVLAPR,NDIM )
           CALL VECSUM(AIVEC(1,JVEC),AIVEC(1,JVEC),AIVEC(1,IVEC),
     &                 1.0D0,-OVLAPR,NDIM )
           CALL VECSUM(ARVEC(1,JVEC),ARVEC(1,JVEC),AIVEC(1,IVEC),
     &                 1.0D0,OVLAPI,NDIM )
           CALL VECSUM(AIVEC(1,JVEC),AIVEC(1,JVEC),ARVEC(1,IVEC),
     &                 1.0D0,-OVLAPI,NDIM )
         END IF
  190    CONTINUE
  200 CONTINUE
  201 CONTINUE
 
*
* Normalize eigenvectors
*
      DO 300 L = 1, NDIM
        XNORM = INPROD(ARVEC(1,L),ARVEC(1,L),NDIM)
     &        + INPROD(AIVEC(1,L),AIVEC(1,L),NDIM)
        FACTOR = 1.0D0/SQRT(XNORM)
        CALL SCALVE(ARVEC(1,L),FACTOR,NDIM)
        CALL SCALVE(AIVEC(1,L),FACTOR,NDIM)
  300 CONTINUE
      END IF
*
* Order eigensolutions after size of real part of A
*
      IF(IORD.EQ.1) THEN
*. Get reorder array
        CALL DUMSORT(ARVAL,NDIM,NDIM,W,SCR)
C            DUMSORT(VEC,NDIM,NELMNT,IREO,ISCR)      
*. Reorder eigenvalues
C  REO_VEC(IREO,VECIN,NDIM,VECOUT,IWAY)
        CALL REO_VEC(W,ARVAL,NDIM,SCR,1)
        CALL COPVEC(SCR,ARVAL,NDIM)
        CALL REO_VEC(W,AIVAL,NDIM,SCR,1)
        CALL COPVEC(SCR,AIVAL,NDIM)
*. Reorder eigenvectors
C  REO_COL_MAT(IREO,AIN,AOUT,NR,NC,IWAY)
        IF(IEIGVC.EQ.1) THEN
          CALL REO_COL_MAT(W,ARVEC,Z,NDIM,NDIM,1)
          CALL COPVEC(Z,ARVEC,NDIM**2)
          CALL REO_COL_MAT(W,AIVEC,Z,NDIM,NDIM,1)
          CALL COPVEC(Z,AIVEC,NDIM**2)
         END IF
*
      END IF
*
      NTEST = 0
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Output from EIGGMT '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Real and imaginary parts of eigenvalues '
        CALL WRTMAT_EP(ARVAL,1,NDIM,1,NDIM)
        CALL WRTMAT_EP(AIVAL,1,NDIM,1,NDIM)
      END IF
*
      IF(NTEST.GE.10.AND.IEIGVC.EQ.1) THEN
        WRITE(6,*) ' real part of eigenvectors '
        CALL WRTMAT(ARVEC,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' imaginary part of eigenvectors '
        CALL WRTMAT(AIVEC,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE COPMAT(AIN,AOUT,NIN,NOUT)
C
C COPY MATRIX AIN OF DIMENSION NIN,NIN INTO
C      MATRIX AOUT OF DIMENSAION NOUT,NOUT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AIN(NIN,NIN)
      DIMENSION AOUT(NOUT,NOUT)
C
      DO 100 J = 1, NOUT
       CALL COPVEC(AIN(1,J),AOUT(1,J),NOUT)
  100 CONTINUE
C
      RETURN
      END

c $Id$
