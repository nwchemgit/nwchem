      SUBROUTINE CHK_S_FOR_SING2(S,NDIM,NSING,X,SCR,SCR2,THRESG)
*
* Check overlap matrix S for singularities.
* A singularity is defined by threshold THRES
*
* Complete diagonalization performed
*
* Jeppe Olsen, Sept. 98 for checking IC-PT program
*. Modified so S is destroyed, but SCR needs just to be a vector of DIM NDIM
*
* Note that the matrix x returning the Eigenvectors are only 
* included for backwards compatibility...
* as the eigenvectors are also returned in S
*
*. Comparison with CHK_S_FOR_SING2:  THRESG added
      IMPLICIT REAL*8(A-H,O-Z)
*. Input : Overlap in complete form
      DIMENSION S(*)
*. Output : Ordered so the last eigenvectors are the 
* nonsingular eigenvalues. 
      DIMENSION X(*)
*. Scratch : two vectors of dimension NDIM
      DIMENSION SCR(*),SCR2(*)
*. Diagonalize S
*. USE eispack TRED2-TQL2 routines to diagonalize 
C  DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
      CALL DIAG_SYMMAT_EISPACK(S,SCR,SCR2,NDIM,IRETURN)
      CALL COPVEC(S,X,NDIM**2)
*. Lowest eigenvalues first, count number of singularities 
      NSING = 0
      IF(THRESG.LE.0.0D0) THEN
        THRES = 1.0D0-10
      ELSE
        THRES = THRESG
      END IF
C?    WRITE(6,*) ' Threshold for singularities ', THRES
        DO I = 1, NDIM
          IF(ABS(SCR(I)).LE.THRES) NSING = NSING + 1
        END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues of metric '
        CALL WRTMAT(SCR,1,NDIM,1,NDIM)
        WRITE(6,*) ' Number of singularities ', NSING
      END IF
*
      RETURN
      END
      SUBROUTINE CHK_S_FOR_SING(S,NDIM,NSING,X,SCR,SCR2)
*
* Check overlap matrix S for singularities 
*
* Complete diagonalization performed
*
* Jeppe Olsen, Sept. 98 for checking IC-PT program
*. Modified so S is destroyed, but SCR needs just to be a vector of DIM NDIM
*
* Note that the matrix x returning the Eigenvectors are only 
* included for backwards compatibility...
* as the eigenvectors are also returned in S
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input : Overlap in complete form
      DIMENSION S(*)
*. Output : Ordered so the last eigenvectors are the 
* nonsingular eigenvalues. 
      DIMENSION X(*)
*. Scratch : two vectors of dimension NDIM
      DIMENSION SCR(*),SCR2(*)
*. Diagonalize S
*. USE eispack TRED2-TQL2 routines to diagonalize 
C  DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
      CALL DIAG_SYMMAT_EISPACK(S,SCR,SCR2,NDIM,IRETURN)
      CALL COPVEC(S,X,NDIM**2)
*. Lowest eigenvalues first, count number of singularities 
      NSING = 0
      THRES = 1.0D-10
        DO I = 1, NDIM
          IF(ABS(SCR(I)).LE.THRES) NSING = NSING + 1
        END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues of metric '
        CALL WRTMAT(SCR,1,NDIM,1,NDIM)
        WRITE(6,*) ' Number of singularities ', NSING
      END IF
*
      RETURN
      END
      SUBROUTINE CHK_IADR_DX
*
* Check IADR_DX function 
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*
      WRITE(6,*) '======================'
      WRITE(6,*) 'CHK_IADR_DX in action '
      WRITE(6,*) '======================'
*
      IDX_SYM = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, NSMOB
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
*
       DO IGAS = 1, NGAS
       NI = NOBPTS(IGAS,ISM)
       I_OFF = IOBPTS(IGAS,ISM)
       IF(IGAS.EQ.1) THEN
          IS_OFF = 1
       ELSE
          IS_OFF = IS_OFF + NOBPTS(IGAS-1,ISM)
       END IF
       DO JGAS = 1, NGAS
       NJ = NOBPTS(JGAS,JSM)
       J_OFF = IOBPTS(JGAS,JSM)
       IF(JGAS.EQ.1) THEN
         JS_OFF = 1
       ELSE
         JS_OFF = JS_OFF + NOBPTS(JGAS-1,JSM)
       END IF
       DO KGAS = 1, NGAS
       NK = NOBPTS(KGAS,KSM)
       K_OFF = IOBPTS(KGAS,KSM)
       IF(KGAS.EQ.1) THEN
         KS_OFF = 1
       ELSE
         KS_OFF = KS_OFF + NOBPTS(KGAS-1,KSM)
       END IF
       DO LGAS = 1, NGAS
C?      WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_ACT',
C?   &               IGAS,JGAS,KGAS,LGAS,IJKL_ACT
         NL = NOBPTS(LGAS,LSM)
         L_OFF = IOBPTS(LGAS,LSM)
         IF(LGAS.EQ.1) THEN
           LS_OFF = 1
         ELSE
           LS_OFF = LS_OFF + NOBPTS(LGAS-1,LSM)
         END IF
         IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
C?       WRITE(6,*)  'IGAS, JGAS, KGAS, LGAS, IJKL_ACT'
C?       WRITE(6,*)   IGAS, JGAS, KGAS, LGAS, IJKL_ACT 
C?       WRITE(6,*)  'IS_OFF,JS_OFF,KS_OFF,LS_OFF',
C?   &                   IS_OFF,JS_OFF,KS_OFF,LS_OFF
*
         IF(IJKL_ACT.EQ.1) THEN
*

          DO L = L_OFF,L_OFF+NL-1
          DO K = K_OFF,K_OFF+NK-1
          DO J = J_OFF,J_OFF+NJ-1
          DO I = I_OFF,I_OFF+NI-1
*. Offsets wrt start of sym
            IREL = I - IS_OFF + 1
            JREL = J - JS_OFF + 1
            KREL = K - KS_OFF + 1
            LREL = L - LS_OFF + 1
C?          WRITE(6,*)  ' IREL,JREL,KREL,LREL',IREL,JREL,KREL,LREL
*
C                IADR_DX(ISM,JSM,KSM,LSM,IOB,JOB,KOB,LOB)
          IADR = IADR_DX(ISM,JSM,KSM,LSM,I,J,K,L)
*
          END DO
          END DO
          END DO
          END DO
*         ^ End of loop over orbitals over given TS
        END IF
*       ^ End if allowed block
       END DO
       END DO
       END DO
       END DO
*      ^ End of loop over gasspaces
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      RETURN
      END
      FUNCTION IADR_DX(ISM,JSM,KSM,LSM,IOB,JOB,KOB,LOB)
*
* Address of dx e(ij,kl)
*
* IOB, JOB, KOB, LOB given relative to start of sym
*    
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Jeppe + Erik,  Sept. 98
*
*. Absolute indeces
      IABS = IBSO(ISM)-1+IOB
      JABS = IBSO(JSM)-1+JOB
      KABS = IBSO(KSM)-1+KOB
      LABS = IBSO(LSM)-1+LOB
*. Types
      IGAS = ITPFSO(IABS)
      JGAS = ITPFSO(JABS)
      KGAS = ITPFSO(KABS)
      LGAS = ITPFSO(LABS)
*. Diagonal block ?
      IF(ISM .EQ.KSM .AND.JSM .EQ.LSM .AND. 
     &   IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS     ) THEN 
        IJ_EQ_KL_BLK = 1
      ELSE
        IJ_EQ_KL_BLK = 0
      END IF
*. Relative to start of TS
      IREL = IABS - IOBPTS(IGAS,ISM) + 1
      JREL = JABS - IOBPTS(JGAS,JSM) + 1
      KREL = KABS - IOBPTS(KGAS,KSM) + 1
      LREL = LABS - IOBPTS(LGAS,LSM) + 1
*. 
*. Offset for  block
C  I_OFF_DX(IOFF,ITRNSP,IXGAS,IXSM,JXGAS,JXSM,KXGAS,KXSM,LXGAS,LXSM,IDX_SYM)
      CALL I_OFF_DX(IOFF,ITRNSP,IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &              1)
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*
      IJ = (JREL-1)*NI + IREL
      KL = (LREL-1)*NK + KREL
      IF(IJ_EQ_KL_BLK.EQ.1) THEN
        IF(IJ.GE.KL) THEN
        IJKL_REL = (KL-1)*NI*NJ - KL*(KL-1)/2 + IJ
        ELSE
        IJKL_REL = (IJ-1)*NK*NL - IJ*(IJ-1)/2 + KL
        END IF
C       IJKL_REL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2 + MIN(IJ,KL)
      ELSE
        IF(ITRNSP.EQ.0) THEN
          IJKL_REL = (IJ-1)*NK*NL + KL
        ELSE
          IJKL_REL = (KL-1)*NI*NJ + IJ
        END IF
      END IF
*
      IADR_DX = IOFF -1 + IJKL_REL
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from IADR_DX '
        WRITE(6,'(A,4I4)') ' IOB,JOB,KOB,LOB', IOB,JOB,KOB,LOB
        WRITE(6,'(A,4I4)') ' ISM,JSM,KSM,LSM', ISM,JSM,KSM,LSM
        WRITE(6,'(A,4I4)') ' IGAS,JGAS,KGAS,LGAS', IGAS,JGAS,KGAS,LGAS
        WRITE(6,*) ' IJ KL IOFF IJKL_REL IADR', IJ, KL, IOFF, IJKL_REL,
     &               IADR_DX
      END IF
*
      RETURN
      END
      SUBROUTINE CASPT2_FOCK(CC_AMP,VEC1,VEC2)
*
* Setup CASPT2 Fock matrix
*
* Initial version, Erik + Jeppe,  Sept. 98
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
*. Input : Vector used to communicate with CI codes (KCC1 vector )
      DIMENSION CC_AMP(*)
*
      DIMENSION VEC1(*),VEC2(*)
*
      IDUMMY = 1
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',2,'CASFOC')
*. Read in CC coefficients from EKD
      CALL CC_FROM_ERIK(CC_AMP)
*. Construct first order variation in SD basis
      I_CC_EXC = 1
      I12 = 2
      CALL MV7(VEC1,VEC2,LUC,LUSC36,0,0)
*. and set up Fock matrix
C          PTFOCK(LU0,LUN,N,ISM,ISPC)
      CALL PTFOCK(LUC,LUSC36,1,ICSM,ICSPC)
*  
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',2,'CASFOC')
*
      RETURN
      END
      SUBROUTINE CC_FROM_ERIK(CC_AMP)
*
* Import (E,e) coefficients from EKD
*
* E.K.D + J.O. Sept.22 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION CC_AMP(*)
*. 
      RETURN
      END
      SUBROUTINE COM_IC_MAT(SICMAT,HICMAT,XICV1,VEC1,VEC2,
     &                       LUREF,N_IC_PARM)
*
* Set up complete matrix over internal contrated configurations
*
* Hamiltonian as well as overlap matrix constructed
*
* For the ICCI/ICPT program, Sept. 98
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION HICMAT(*), SICMAT(*)
*. scratch
      DIMENSION XICV1(*)
      DIMENSION VEC1(*),VEC2(*)
*
      IDUMMY = 0
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'COM_IC')
*
      WRITE(6,*) ' H0 and S will be constructed '
      WRITE(6,*) ' Dimension :  ', N_IC_PARM
      DO IC = 1, N_IC_PARM
*
        ZERO = 0.0D0
        CALL SETVEC(XICV1,ZERO,N_IC_PARM)
        XICV1(IC) = 1.0D0
*.  To vectors to be calculated : S*c, H*c
        IMOFF = 1 + (IC-1)*N_IC_PARM
        CALL IC_H_T_REF(LUREF,XICV1,HICMAT(IMOFF),SICMAT(IMOFF),
     &  VEC1,VEC2,1,1,0)
C     IC_H_T_REF(LUREF,VEC,HVEC,SVEC,VEC1,VEC2,IDOHV,IDOSV) 
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Hamiltonian matrix over internal contracted confs'
        CALL WRTMAT(HICMAT,N_IC_PARM,N_IC_PARM,N_IC_PARM,N_IC_PARM)
        WRITE(6,*) ' Overlap     matrix over internal contracted confs'
        CALL WRTMAT(SICMAT,N_IC_PARM,N_IC_PARM,N_IC_PARM,N_IC_PARM)
      END IF
*
      RETURN
      END
      SUBROUTINE SOLVE_PT2_EQ(RHS,C1,VEC1,VEC2,H0,ENOT,S,CC_AMP)
*
* Solve first order pt equations
*
* (H0-E0) C1 = - RHS
*
* Initial version, H0 and S has been constructed,
*
* Jeppe Olsen, Sept. 98
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD, INPRDD
*. Specific Input
      DIMENSION RHS(*)
*. General input
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'oper.inc'
*. Specific input
      DIMENSION H0(*),S(*)
*. Output
      DIMENSION C1(*)
*. Vector used to communicate with CI
      DIMENSION CC_AMP(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*   
      IDUMMY = 1
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'IC_PT2')
*
      NTEST = 1
*  
      NDIM  = N_CC_AMP
      NDIMO = NDIM
      LEN   = N_CC_AMP ** 2
*
      CALL MEMMAN(KLH0E0S,LEN,'ADDL  ',2,'H0E0S ')
      CALL MEMMAN(KLSCR  ,LEN,'ADDL  ',2,'LSCR  ')
      CALL MEMMAN(KLSCR2 ,LEN,'ADDL  ',2,'LSCR2 ')
      CALL MEMMAN(KLSCR3 ,LEN,'ADDL  ',2,'LSCR3 ')
      CALL MEMMAN(KLVEC  ,3*LEN,'ADDL  ',2,'LVEC  ')
      CALL MEMMAN(KLRHS2 ,LEN,'ADDL  ',2,'RHS2  ')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input RHS '
        CALL WRTMAT(RHS,1,NDIM,1,NDIM)
      END IF
*. Check for singularities
C     CHK_S_FOR_SING(S,NDIM,NSING,X,SCR)
      CALL CHK_S_FOR_SING(S,NDIM,NSING,WORK(KLSCR),WORK(KLSCR2))
      KLXNONSING = KLSCR + NSING*NDIM
      NDIM2 = NDIM - NSING
*. Eigenvectors are in KLSCR, first NSING are singularities
*. Obtain H0 - E0 S
      ENOTM = -ENOT
      ONE = 1.0D0
      CALL VECSUM(WORK(KLH0E0S),H0,S,ONE,ENOTM,LEN)
*. Eliminate singularities
C     IF(NSING.NE.0) THEN
*. Transform H0-E0S to nonsingular basis
        CALL TRAN_SYM_BLOC_MAT3(WORK(KLH0E0S),WORK(KLXNONSING),
     &       1,NDIM,NDIM2,WORK(KLSCR2),WORK(KLVEC),0) 
        CALL COPVEC(WORK(KLSCR2),WORK(KLH0E0S),NDIM2**2)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' H0-E0 S in non singular basis '
          CALL WRTMAT(WORK(KLH0E0S),NDIM2,NDIM2,NDIM2,NDIM2)
        END IF
*. Transform RHS
        ZERO = 0.0D0
        ONE = 1.0D0
        CALL MATML7(WORK(KLRHS2),WORK(KLXNONSING), RHS,
     &              NDIM2,1,NDIM,NDIM2,NDIM,1,ZERO,ONE,1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Transformed RHS '
          CALL WRTMAT(WORK(KLRHS2),NDIM2,1,NDIM2,1)
        END IF
        
        NDIM = NDIM2
C     ELSE
C       CALL COPVEC(RHS,WORK(KLRHS2),NDIM)
C     END IF
*. Invert  
      CALL INVERT_BY_DIAG(WORK(KLH0E0S),WORK(KLSCR3),
     &                    WORK(KLSCR2),WORK(KLVEC),NDIM)
*. Multiply rhs with inverse mat to get -1 times first order correction 
      CALL MATVCB(WORK(KLH0E0S),WORK(KLRHS2),WORK(KLVEC),NDIM,NDIM,0)
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLVEC),ONEM,NDIM)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' First order correction  vector in nonsing basis'
        CALL WRTMAT(WORK(KLVEC),1,NDIM,1,NDIM)
      END IF
*. Transform first order correction to (E,e) basis 
      ZERO = 0.0D0
      ONE  = 1.0D0 
      CALL MATML7(WORK(KLSCR2),WORK(KLXNONSING),WORK(KLVEC),
     &            NDIMO,1,NDIMO,NDIM2,NDIM2,1,ZERO,ONE,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' First order correction  vector in (e,E) basis'
        CALL WRTMAT(WORK(KLSCR2),1,NDIMO,1,NDIMO)
      END IF
*. Set up first order correction in det basis
      CALL COPVEC(WORK(KLSCR2),CC_AMP,NDIMO)
      ICC_EXC = 1
      I12 = 2
      CALL MV7(VEC1,VEC2,LUC,LUHC,0,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' First order correction in SD basis '
        CALL WRTVCD(VEC1,LUHC,1,-1)
      END IF
*. H times first order correction 
      ICC_EXC = 0
      I12 = 2
      CALL MV7(VEC1,VEC2,LUHC,LUSC1,0,0)
*. Subtract E0*!0>
      FACTOR = -ENOT
      CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LUSC1,LUREF,LUHC,1,LBLK)
      WRITE(6,*) ' (H-E0)!CASPT> '
      CALL WRTVCD(VEC1,LUHC,1,-1)
*. Inner product with zero-order vector
      E2 = INPRDD(VEC1,VEC2,LUSC1,LUC,1,-1)
      WRITE(6,*) ' Second order energy correction ', E2
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'IC_PT2')
      RETURN
      END
*
      SUBROUTINE LUCIA_ICCI_OLD(ISM,ISPC,CALC)
*
* Internal contracted calculation with ICCI
*    
* Jeppe Olsen, Initiated in May 98
*              Finished before retirement - I hope
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      CHARACTER*6 CALC
      INCLUDE 'cands.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
*
      WRITE(6,*) '*************************'
      WRITE(6,*) '*                       *'
      WRITE(6,*) '* Welcome to LUCIA_ICCI *'
      WRITE(6,*) '*                       *'
      WRITE(6,*) '*************************'
      WRITE(6,*)
      WRITE(6,*)
      IF(CALC(1:4).EQ.'ICCI') THEN
        WRITE(6,*) ' Internal contracted CI calculation '
      ELSE IF (CALC(1:4) .EQ. 'ICPT' ) THEN
        WRITE(6,*) ' Internal contracted PT calculation '
      END IF
*
* We will write the ICCI wave function as
*
*    !0'> = sum_n c_n !n> + P{ sum_{pq} C_{pq} E_{pq} !0> +
*                          sum_{(pq).ge.(rs)} !0> C_{pq,rs} E_{pq}E_{rs}!0>}
*         = sum_{mu} C_{mu} t_{mu} !0>
*
* The sum over dets corresponds to dets included in  w.f. space ISPC-1,
* whereas the excitation generates the part of the w.f. that is in
* space ISPC, but not ISPC -1
*
* The projection operator P projects out determinants belonging to
* space ISPC - 1
*
* P = Sum_{n belonging to ISPC-1} !n><n!
*
*
* 1 : analyze which single excitations are out of space excitations
*
*. a Set up allowed combinations of supergroups for space ISPC-1
*  b Set up allowed combinations of supergroups for space ISPC
*  c Set up matrix of connections between supergroups
*  d From the above, divide single excitations into inspace and out-of-space
*    excitations
*
* 2. Prepare for Internal contracted PT2 calculation
*    a : set up <mu!H!0>
*    b : set up -approximate- diagonal  <mu!H_0!mu>
* 3. Solve first order equations  
*    a :( search for singularities )
*    b : invoke the linear equation solver  
*
      IDUMMY = 0    
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'ICCI  ')
      STOP ' Entering old code '
      WRITE(6,*) ' Symmetry and space ', ISM,ISPC
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
*. Complete operator 
      I12 = 2
*. Coupled cluster flag
      ICC_EXC = 1
*. Divide orbital space into inactive, active, secondary
      CALL CC_AC_SPACES
*. Number of single excitation and double excition amplitudes
      IEXSYM = 1
      CALL FIND_N_CC_AMP
     &(IEXSYM,NSXA,NSXB,NDXAA,NDXBB,NDXAB,NNSXE,NNDXEE)
*
      NSXE = NNSXE
      NDXEE = NNDXEE
*. (NSXE and NDXEE are stored in /CRUN/
*. Allocate space for 3 CC vectors
* : Three type of operators
*  1 : single excitations
*  2 : Double excitations
*  3 : Determinants orthogonal to reference (I guess i am going multiref cc)
* The amplitudes will be stored as above
      NNDET = 0
      N_CC_AMP = NSXE + NDXEE + NNDET
      WRITE(6,*) ' NSXE NDXEE N_CC_AMP', NSXE,NDXEE,N_CC_AMP
      LEN = N_CC_AMP
      CALL MEMMAN(KCC1,2*LEN,'ADDL  ',2,'CC1  ')
      CALL MEMMAN(KCC2,2*LEN,'ADDL  ',2,'CC2  ')
      CALL MEMMAN(KCC3,2*LEN,'ADDL  ',2,'CC3  ')
      CALL MEMMAN(KRHS,  LEN,'ADDL  ',2,'RHS  ')
      CALL MEMMAN(KLRHO1S,NTOOB**2,'ADDL  ',2,'RHO1SA')
*. save one-particle density matrix
      CALL COPVEC(WORK(KRHO1),WORK(KLRHO1S),NTOOB**2)
*. Scratch space for CI - behind the curtain 
       CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)
*. Initialize CC amplitudes in CC1
      CALL INI_CC_AMP(WORK(KCC1),1)
*. Zero order operator
      IF(CALC(1:4).EQ.'ICPT' ) THEN
        ICC_EXC = 1
*Der er her det sker
*. Call subroutine that reads in double coefs.
*. In call use work(KCC1+NSXE) as vector for storing CASPT2 coefs.
*
C       CALL CHK_IADR_DX
C            CASPT2_FOCK(CC_AMP,VEC1,VEC2)
C       CALL CASPT2_FOCK(WORK(KCC1),WORK(KVEC1),WORK(KVEC2))
C       WRITE(6,*) ' Enforced (Jeppe again!) stop in ICPT '
C       STOP       ' Enforced (Jeppe again!) stop in ICPT '
*
*. Obtain zero order operator  F = FI + FA
* ( density is assumed stored in WORK(KRHO1))
*
        CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
        CALL FIFAM(WORK(KFI))
        CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ECORE_H = 0.0D0
*. Set up rhs 
        ZERO = 0.0D0
        CALL SETVEC(WORK(KCC1),ZERO,N_CC_AMP)
        CALL IC_H_T_REF(LUC,WORK(KCC1),WORK(KRHS),WORK(KCC3),
     &                  WORK(KVEC1),WORK(KVEC2),1,1,1)
C?      WRITE(6,*) ' Fresh RHS '
C?      CALL WRT_CC_VEC(WORK(KRHS),6)
*. Initial testing : Set up H and S matrices
        LEN = N_CC_AMP**2
        CALL MEMMAN(KLSMAT,LEN,'ADDL  ',2,'SICMAT')
        CALL MEMMAN(KLHMAT,LEN,'ADDL  ',2,'HICMAT')
*. One-electron operator
        I12 = 1
*. Fifa as one-electron operator 
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
*. Zero order energy
        CALL COPVEC(WORK(KLRHO1S),WORK(KRHO1),NTOOB**2)
        ICC_EXC = 0
        CALL EN_FROM_DENS(EZERO,1,0)
        WRITE(6,*) ' zero-order energy ', EZERO
        CALL COM_IC_MAT(WORK(KLSMAT),WORK(KLHMAT),WORK(KCC1),
     &                  WORK(KVEC1),WORK(KVEC2),LUC,N_CC_AMP)
C            SOLVE_PT2_EQ(RHS,C1,VEC1,VEC2,H0,ENOT,S)
*. Clean up 
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        CALL SOLVE_PT2_EQ(WORK(KRHS),WORK(KCC3),WORK(KVEC1),
     &                    WORK(KVEC2),WORK(KLHMAT),EZERO,WORK(KLSMAT),
     &                    WORK(KCC1))
        CALL CASPT2_FOCK(WORK(KCC1),WORK(KVEC1),WORK(KVEC2))
C       CALL CASPT_CHECK
        WRITE(6,*) ' Enforced (Jeppe again!) stop in ICPT '
        STOP       ' Enforced (Jeppe again!) stop in ICPT '
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'ICCI  ')
      RETURN
      END
      SUBROUTINE IC_H_T_REF(LUREF,VEC,HVEC,SVEC,VEC1,VEC2,IDOHV,IDOSV,
     &                       IADDREF) 
*
* Matrix vector routine for Internal contracted formulation
* Note : reference vector is assumed located on LUREF
*
* IF IADREF = 1, the reference state is added to T |ref>
*

*
* Hvec_i = <ref! (O_IC)_{i}^{\dag} H sum_j (O_IC)_j !ref> Vec_j 
* Svec_i = <ref! (O_IC)_{i}^{\dag}   sum_j (O_IC)_j !ref> Vec_j 
*
* IDOHV = 1 =>  Calc Hvec
* IDOSV = 1 =>  Calc Svec
*
* <0!E(ij)
* <0!e(ijkl)
* <K!
*
* No transformation to biorthonormal basis performed
*
* Jeppe Olsen, September of 98
*
* Initial version     
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Input : Must be the vector used to communicate IC coefficients 
      DIMENSION VEC(*)
*. Output
      DIMENSION HVEC(*),SVEC(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*. 
      INCLUDE 'cc_exc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'oper.inc'
*
C     COMMON/CC_EXC/ICC_EXC
*
      IDUMMY = 0
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'IC_H_T')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' Info from IC_H_T_REF ' 
        WRITE(6,*) ' ====================='
        WRITE(6,*)
        WRITE(6,*) ' Input (E,e) parameters'
        ICC_EXC = 1
        CALL WRT_CC_VEC(VEC,6)
      END IF
*
      I12_SAVE = I12
      LBLK = -1
*. Calculate  sum_j (O_IC)_j !ref> Vec_j 
*. (Notice the IC coefficients are sneaked  in through the back door !
      ICC_EXC = 1
      I12_SAVE = I12
      I12 = 2
      CALL MV7(VEC1,VEC2,LUREF,LUHC,0,0)
      I12 = I12_save
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Slater determinant expansion parameters '
        CALL WRTVCD(VEC1,LUHC,1,-1)
      END IF
*
      IF(IADDREF.EQ.1) THEN
        ONE = 1.0D0
        CALL VECSMD(VEC1,VEC2,ONE,ONE,LUREF,LUHC,LUSC1,1,LBLK)
        CALL COPVCD(LUSC1,LUHC,VEC1,1,LBLK)
      END IF
* result is now on LUHC
* Svec
      IF(IDOSV.EQ.1) THEN
* <IC*ref! (E,e) !ref>
       IPRDEN_ORIG = IPRDEN
       IPRDEN = 0
       XDUM = 0.0D0
       CALL DENSI2(2,WORK(KRHO1),WORK(KRHO2),VEC1,VEC2,
     &             LUHC,LUC,EXPS2,0,XDUM,XDUM,XDUM,XDUM,1)
       IPRDEN = IPRDEN_ORIG
*. Reorganize into IC vector format
       IBIO = 0
       CALL REF_SX(WORK(KRHO1),SVEC(1),2,1,1,IBIO)
       CALL REF_DX_EXP_COMP(WORK(KRHO2),SVEC(1+NSXE),1,1,IBIO)
      END IF
*
      IF(IDOHV.EQ.1) THEN
*. H times !Ic*ref>
       ICC_EXC = 0
* <IC*ref!H (E,e) !ref> on LUSC1
       CALL MV7(VEC1,VEC2,LUHC,LUSC1,0,0)
       IPRDEN_ORIG = IPRDEN
       IPRDEN = 0
       CALL DENSI2(2,WORK(KRHO1),WORK(KRHO2),VEC1,VEC2,
     &             LUSC1,LUC,EXPS2,0,XDUM,XDUM,XDUM,XDUM,1) 
       IPRDEN = IPRDEN_ORIG
*. Reform
       IBIO = 0
       ICC_EXC = 1
       CALL REF_SX(WORK(KRHO1),HVEC(1),2,1,1,IBIO)
       CALL REF_DX_EXP_COMP(WORK(KRHO2),HVEC(1+NSXE),1,1,IBIO)
      END IF
*
      IF(NTEST.GE.100) THEN
CM      WRITE(6,*)
CM      WRITE(6,*) ' ====================='
CM      WRITE(6,*) ' Info from IC_H_T_REF ' 
CM      WRITE(6,*) ' ====================='
CM      WRITE(6,*)
*
        IF(IDOSV.EQ.1) THEN
          WRITE(6,*) ' <ref!(E,e) IC*!ref >'
          CALL WRT_CC_VEC(SVEC,6)
        END IF
*
        IF(IDOHV.EQ.1) THEN
          WRITE(6,*) ' <ref! (E,e) H IC* !ref >'
          CALL WRT_CC_VEC(HVEC,6)
        END IF
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'IC_H_T')
*
      RETURN
      END
      SUBROUTINE CASPT_CHECK
*
* Check CASPT program for He, 2s cas, 1s inact.
*
*The 2 x 2 cas wave function is
*
*!0> = c11  !11> + c12 ( | 12 > + | 21> ) + c22 | 22 >
*
*The order of the SD's are
*
*!33>
*!13> !23>
*!31> !32>
*!11> !12> !21> !22>
*
* Denoting 
*
* !13> + !31> = !1>
* !23> + !32> = !2>
* !33>        = !3>
*and the excitations become
*
*E31 !0> = c11 (!13> + !31> ) + c12 (!32> + | 23>)  (1)
*E32 !0> = c12 (!13> + !31> ) + c22 (!32> + | 23>)  (2)
*e31,11 !0> = c11 (!13> + !31> )  (3)
*e32,11 !0> = c21 (!31> + !13> )  (4)
*e31,21 !0> = c11 (!32> + !23> )  (5)
*e32,21 !0> = c21 (!32> + !23> )  (6)
*e31,12 !0> = c21 (!31> + !13> )  (7)
*e32,12 !0> = c22 (!31> + !13> )  (8)
*e31,22 !0> = c21 (!32> + !23> )  (9)
*e32,22 !0> = c22 (!23> + !32> )  (10)
*1/2e31,31 !0> = c11 !33>  (11)
*e32,31 !0> = 2c21 !33>  (12)
*1/2e32,32 !0> = c22 !33>  (13)
*
* Trans for IC to dets
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION C(600)
      DIMENSION F(600)
      DIMENSION FIC(600)
      DIMENSION SCR(600)
      DIMENSION S(600)
      DIMENSION SIC(600)
*
      C11 = 0.99841279E+00 
      C12 =-0.48896109E-02 
      C22 =-0.55893556E-01
*
      NIC = 13
      NDET = 3
      ZERO = 0.0D0
      CALL SETVEC(C,ZERO,NIC*NDET)
*
      C((3-1)*NDET + 1 ) = C11
      C((4-1)*NDET + 1 ) = C12
      C((5-1)*NDET + 2 ) = C11
      C((6-1)*NDET + 2 ) = C12
      C((7-1)*NDET + 1 ) = C12
      C((8-1)*NDET + 1 ) = C22
      C((9-1)*NDET + 2 ) = C12
      C((10-1)*NDET + 2 ) = C22
      C((11-1)*NDET + 3 ) = C11
      C((12-1)*NDET + 3 ) =2*C12
      C((13-1)*NDET + 3 ) = C22
*
      C((1-1)*NDET + 1 ) = C11
      C((1-1)*NDET + 2 ) = C12
*
      C((2-1)*NDET + 1 ) = C12
      C((2-1)*NDET + 2 ) = C22
 
* The f matrix over orbitals
      FINT11 = -.8905437E+00
      FINT12 = 0.5852368E-03 
      FINT22 = 0.1071525E+01
      FINT13 = 0.9153308E-03 
      FINT23 = 0.2102955E-02 
      FINT33 = 0.6347683E+01
* the f matrix over combinations (as sym packed row order matrix)
      F(1) = 2*(FINT11 + FINT33)
      F(2) = 2* FINT12 
      F(3) = 2*(FINT22 + FINT33)
      F(4) = 2* FINT13
      F(5) = 2* FINT23
      F(6) = 2* FINT33
*. The s matrix over combinations
      CALL SETVEC(S,ZERO,NDET*(NDET+1)/2)
      S(1) = 2.0D0
      S(3) = 2.0D0
      S(6) = 1.0D0
*
*. Transform F to internal contracted form
*
C       TRAN_SYM_BLOC_MAT3
C    &(AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      WRITE(6,*) ' Transformation matrix '
      CALL WRTMAT(C,NDET,NIC,NDET,NIC)
      WRITE(6,*) ' H0 in det(comb) basis '
      CALL PRSYM(F,NDET)    
      CALL TRAN_SYM_BLOC_MAT3(F,C,1,3,13,FIC,SCR,1)
*. Outpack
      CALL TRIPAK(SCR,FIC,2,NIC,NIC)
      WRITE(6,*) ' H0 in IC basis '                
      CALL WRTMAT(SCR,NIC,NIC,NIC,NIC)
*
*. Transform S to internal contracted form
*
C     TRAN_SYM_BLOC_MAT3
C    &(AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL TRAN_SYM_BLOC_MAT3(S,C,1,3,13,SIC,SCR,1)
*. Outpack
      CALL TRIPAK(SCR,SIC,2,NIC,NIC)
      WRITE(6,*) ' S in IC basis '                
      CALL WRTMAT(SCR,NIC,NIC,NIC,NIC)
*
      RETURN
      END
      SUBROUTINE DIAG_SYM_MAT(A,X,SCR,NDIM,ISYM)
*
* Diagonalize symmetric matrix
* (outer shell to eigen)
*
* On output X contains eigenvectors and SCR contains eigenvalues
*
* Jeppe Olsen, Sept 98                 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input  matrix
      DIMENSION A(*)       
*. Output
      DIMENSION X(NDIM**2)
*. Output and scratch 
      DIMENSION SCR(NDIM*(NDIM+1)/2)
*  
      NTEST = 0
*. Reform to packed matrix
      IF(ISYM.EQ.0) THEN
        CALL TRIPAK(A,SCR,1,NDIM,NDIM)
      ELSE 
        CALL COPVEC(A,SCR,NDIM*(NDIM+1)/2)
      END IF
*. Diagonalize
      CALL EIGENL(SCR,X,NDIM,0,1)
*. Pack eigenvalues
      CALL COPDIA(SCR,SCR,NDIM,1)
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT(SCR,NDIM,1,NDIM,1)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvectors of matrix '
        CALL WRTMAT(X,NDIM,NDIM,NDIM,NDIM)     
      END IF
*
      RETURN
      END 



c $Id$
