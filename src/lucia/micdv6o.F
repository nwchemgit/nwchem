      SUBROUTINE MICDV6O(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,
     &                  MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LU8,LUDIA,NROOT,MAXVEC,
     &                  NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,IROOTHOMING,NCNV_RT,RTCNV,IPRECOND,
     &                  CONVER,RNRM_CNV,ISBSPPR_ACT)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version 
*
* Version with 2 or 3 vectors in subspace per root
*
* From MICDV4, Jeppe Olsen, April 1997
*                           Oct 97: root homing added
*                           May 99: New preconditioner switch 
*                           Feb.13: Extensive recoding for improved stability and readability..
*                                   Lu8 added (for compatibility)
*
* Input :
* =======
*        MV7: Sigma routine 
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
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
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       LOGICAL RTCNV(NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
       DIMENSION RNRM_CNV(NROOT)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*MAXVEC                               
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
*      IROOTHOMING : Do roothoming, i.e. select the
*      eigenvectors in iteration n+1 as the approximations
*      with largest overlap with the previous space
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER
       REAL*8 INPRDD, INPROD
*. Notice XJEP is also used for ROOTHOMING, should be allocated
* outside (for roothoming :dim = 3*MAXVEC )
       DIMENSION XJEP(10000)
       INTEGER   IXJEP(10000)
*. Dimension of Work2: 3*MAXVEC**2
       DIMENSION WORK2(10000)
*. For working with routines accessing several files
       DIMENSION LU_FILE1(5),NVC_FILE1(5)
       DIMENSION LU_FILE2(5),NVC_FILE2(5)
*
       EXTERNAL MV7
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
       IPRT = 10
       IPRT = MAX(IPRTXX,IPRT)
       IF(IPRT.GE.0) WRITE(6,*) ' MICDV6 in action, NROOT, MAXVEC =',
     &                            NROOT,MAXVEC
       WRITE(6,*) ' TESTY, ISBSPPR_ACT = ', ISBSPPR_ACT
*
       IOLSTM = 0
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV6 '
         STOP 20
       END IF
*
C?     WRITE(6,*) ' MICDV6 : NROOT, IPRECOND ', NROOT, IPRECOND
C?     WRITE(6,*) ' LU1 LU2 LU3 LU4 LU5 LU6',
C?   &              LU1,LU2,LU3,LU4,LU5,LU6 
       IF(IROOTHOMING.EQ.1) THEN
         WRITE(6,*) ' Root homing performed '
       END IF
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = -1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       CALL GFLUSH(6)
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       DO IVEC = 1,NINVEC
         CALL REWINO(LU5)
         CALL REWINO(LU6)
         CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU5,LU6,0,0)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU6)
         CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
         DO JVEC = 1, IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
         END DO
*        ^ End of loop over JVEC
       END DO
*      ^ End of loop over IVEC
*
       IF(IPRT.GE.10000) THEN
         CALL REWINO(LU1)
         CALL REWINO(LU2)
         DO IVEC = 1, NINVEC
           WRITE(6,*)
           WRITE(6,*) ' ==============================='
           WRITE(6,*) ' Initial C and S vector ', IVEC
           WRITE(6,*) ' ==============================='
           WRITE(6,*)
           CALL WRTVCD(VEC1,LU1,0,LBLK)
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
     
       IF( IPRT .GE.5 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,*) ' Eigenvalues of initial iteration (with shift)'
         WRITE(6,'(5F21.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(6,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
       CALL REWINO(LU3)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINO(LU3)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)
*
       IF(IPRT.GE.10000) THEN
         WRITE(6,*) ' Initial set of eigenvectors '
         CALL REWINO(LU1)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(6,*) ' Initial set of sigma vectors '
         CALL REWINO(LU2)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       ZERO = 0.0D0
       CALL SETVEC(APROJ,ZERO,NROOT*(NROOT+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO
*
       NVEC = NROOT 
CJAN20 IF (MAXIT .EQ. 1 ) GOTO  901
C      IF (MAXIT .EQ. 1 ) GOTO  1001
*
* ======================
*. Loop over iterations
* ======================
*
*
 1000 CONTINUE
       ITER = ITER + 1
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*) ' Info from iteration .... ', ITER
       END IF
       CALL GFLUSH(6)
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       CALL REWINO(LU1)
       CALL REWINO(LU2)
*
       DO 100 IROOT = 1, NROOT
         IF(IPRT.GE.100) WRITE(6,*) ' Loop 100, IROOT = ', IROOT
*
*
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINO(LU7)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
*
         CALL REWINO(LU7)
         CALL REWINO(LU5)
         FACTOR = - EIGAPR
         CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU2,LU7,LU5,0,LBLK)
*
         IF ( IPRT  .GE. 10000 ) THEN
           WRITE(6,*) '  ( HX - EX ) for IROOT =', IROOT
           CALL WRTVCD(VEC1,LU5,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         RNORM = SQRT( ABS(INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)) )
         RNRM(ITER-1,IROOT) = RNORM
C?       WRITE(6,*) ' MICDV6: ITER = ', ITER
*
         IF(IPRT.GE.1000) WRITE(6,*) ' RNORM = ', RNORM
         IF(RNORM.LT. TEST .OR. 
     &      (ITER.GT.2.AND.
     &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            IF(IROOT.LE.NCNV_RT) CONVER = .FALSE.
         END IF
*
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU6
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(6,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU5)
           CALL REWINO(LU6)
           CALL H0M1TD(LU6,LUDIA,LU5,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND,ISBSPPR_ACT)
           IF ( IPRT  .GE. 1000) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU6,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU5
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND,ISBSPPR_ACT)
*. Get X back on LU7
             CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
             CALL REWINO(LU7)
             CALL COPVCD(LU1,LU7,VEC1,0,LBLK)

*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             GAMMA = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
             IF(IPRT.GE.10) WRITE(6,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINO(LU5)
             CALL REWINO(LU7)
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU7,LU5,0,LBLK))
              IF(IPRT.GE.1000) write(6,*) ' VNORM = ', VNORM
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                FACTOR = -(DELTA/GAMMA)
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU6,LU5,LU7,1,LBLK)
                CALL COPVCD(LU7,LU6,VEC1,1,LBLK)
                XNORM2 = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
                XNORM_INI = SQRT(ABS(XNORM2))
*
                IF(IPRT.GE.10000) THEN
                  WRITE(6,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU6,1,LBLK)
                END IF
*
              END IF
            END IF
*
*. 1.3 Orthogonalize to all previous vectors
*
*
* MGS- more I/O intensive and stable
*
C              ADD_ORTN_2SUBSPC(LU1,LU2,N1,N2,LUNEW,IMNEW,
C    &                          VEC1,VEC2,LUSC1,LUSC2,THRES)
            THRES = 1.0D-3 
            CALL ADD_ORTN_2SUBSPC(LU1,LU3,NROOT,NVEC+IADD-1-NROOT,
     &           LU6,IMNEW,VEC1,VEC2,LU7,LU5,THRES)
            I_SUPER_ORT = 0
            IF(I_SUPER_ORT.EQ.1) THEN
*. A second orthonormalization 
              THRES = 1.0D-3
              CALL ADD_ORTN_2SUBSPC(LU1,LU3,NROOT,NVEC+IADD-1-NROOT,
     &             LU6,IMNEW,VEC1,VEC2,LU7,LU5,THRES)
            END IF
*
            IF(IMNEW.EQ.1) THEN
*. Vector is new, copy to LU3- which was left at EOF by add
              CALL REWINO(LU6)
              CALL COPVCD(LU6,LU3,VEC1,0,LBLK)
            ELSE
*. Nothing new 
              IADD = IADD - 1
              IF(IPRT.GE.3) THEN
               WRITE(6,*) ' Correction vector removed for root',IROOT
              END IF
            END IF
         END IF ! Root is not converged
  100  CONTINUE
*. IADD new vectors have been obtained, NROOT are on LU1, IADD + NVEC-NROOT on LU3
       IF(IPRT.GE.10) 
     &  WRITE(6,*) ' Number of orthogonal corrections added ',  IADD
       NADD1 = IADD
       NVEC3 = NVEC-NROOT+IADD
*. Notice: NVEC is not updated yet, it is still the number of vectors from prev. it
*
       I_CHECK_OVLAP = 1
       IF(I_CHECK_OVLAP.EQ.1) THEN
*. Overlap LU1 LU2
         WRITE(6,*) ' Overlap roots'
C             GET_OVLAP_FOR_VECTORS(LU1,N1,LU2,N2,LUSCR,OVLAP,VEC1,VEC2)
         CALL GET_OVLAP_FOR_VECTORS(LU1,NROOT,LU1,NROOT,LU7,
     &        XJEP,VEC1,VEC2)
         CALL CHECK_UNIT_MAT(XJEP,NROOT,XMAX_DIA,XMAX_OFD,1)
C             CHECK_UNIT_MAT(UNI,NDIM,XMAX_DIFF_DIAG,XMAX_DIFF_OFFD,ISYM)
 
         WRITE(6,*) ' Overlap roots and corrections '
         CALL GET_OVLAP_FOR_VECTORS(LU1,NROOT,LU3,NVEC3,LU7,
     &        XJEP,VEC1,VEC2)
         XMAX = XMNMX(XJEP,NROOT*NVEC3,3)
         WRITE(6,*) ' Largest element of <roots!Corrections> ', XMAX
         WRITE(6,*) ' Overlap corrections '
         CALL GET_OVLAP_FOR_VECTORS(LU3,NVEC3,LU3,NVEC3,LU7,
     &        XJEP,VEC1,VEC2)
         CALL CHECK_UNIT_MAT(XJEP,NVEC3,XMAX_DIA,XMAX_OFD,1)
       END IF
       IF(CONVER) THEN
*. Well, we converged     
C         ITER = ITER-1
          GOTO  1001
       END IF
       IF( ITER.GT. MAXIT) THEN
          ITER = MAXIT
          GOTO 1001
       END IF
*
* ====================================================
*  2 : Optimal combination of new and old directions
* ====================================================
*
*  2.1: Multiply new directions with matrix
*
       IF(IPRT.GE.10000) THEN
         WRITE(6,*) ' Vectors on LU3'
         WRITE(6,*) ' Number of vectors = ', NVEC3
         CALL REWINO(LU3)
         DO IVEC = 1, NVEC3
           CALL WRTVCD(VEC1,LU3,0,LBLK)
         END DO 
       END IF
*
       CALL SKPVCD(LU3,NVEC-NROOT,VEC1,1,LBLK)
       CALL SKPVCD(LU4,NVEC-NROOT,VEC1,1,LBLK)
       DO IVEC = 1, NADD1
         CALL REWINE(LU5,LBLK)
         CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU5,LU6,0,0)
         CALL REWINE(LU6,LBLK)
         CALL COPVCD(LU6,LU4,VEC1,0,LBLK)
*.2.2 Augment projected matrix
         CALL REWINE( LU1,LBLK)
         DO JVEC = 1, NROOT
           CALL REWINE(LU6,LBLK)
           IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
         END DO
*
         CALL REWINE(LU3,LBLK)
         DO JVEC = NROOT+1, NVEC+IVEC
          CALL REWINE(LU6,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
         END DO
       END DO
*      /\ End do over new trial vectors
*. 2.3 Diagonalize projected matrix
       NVEC = NVEC + IADD
       CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
       CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
*. Test : transform projected matrix
C TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
C      CALL TRAN_SYM_BLOC_MAT(APROJ,AVEC,1,NVEC,XJEP(1000),XJEP(1))
C      WRITE(6,*) ' Explicitly transformed matrix'
C      CALL PRSYM(XJEP(1000),NVEC)

       IF(IPICO.NE.0) THEN
         E0VAR = WORK(KAPROJ)
         C0VAR = AVEC(1)
         C1VAR = AVEC(2)
         C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
         AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
         AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
         E0PERT = AVEC(1)**2*APROJ(1)
     &          + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &          + AVEC(2)**2*APROJ(3)
         WORK(KAPROJ) = E0PERT
         IF(IPRT.GE.10) THEN
          WRITE(6,*) ' Var and Pert solution, energy and coefficients'
          WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
          WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
         END IF
       END IF
*
       IF(IROOTHOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first 
*
*. Norm of wavefunction in previous space
        DO IVEC = 1, NVEC
          XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                 AVEC(1+(IVEC-1)*NROOT),NROOT)
        END DO
        IF(IPRT.GT.10) WRITE(6,*) 
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
        IF(IPRT.GE.3) THEN
          WRITE(6,*) ' New roots choosen as vectors '
          CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
        END IF
*. Reorder
        DO INEW = 1, NVEC
          IOLD = IXJEP(INEW)
          CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
        END DO
        CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
        DO INEW = 1, NVEC
          IOLD = IXJEP(INEW)
          XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
        END DO
        DO INEW = 1, NVEC
          WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
        END DO
*
        IF(IPRT.GE.3) THEN
          WRITE(6,*) ' Reordered WORK and AVEC arrays '
          CALL PRSYM(WORK,NVEC)
          CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
        END IF
*
       END IF
*      ^ End of root homing procedure
*
       DO IROOT = 1, NROOT
         EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F21.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(6,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
       IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,F22.13)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
       END IF
*
**  Reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset      
*
*
* case 1 : Only NROOT vectors can be stored
*          save current eigenvector approximations
* Case 2 : Atleast 2*NROOT eigenvectors can be saved
*          Current eigenvactor approximations+
*          vectors allowing generation of previous approxs.
*
*
C       IF(IPRT.GE.1000) THEN 
        IF(IPRT.GE.000) THEN 
          write(6,*) ' I am going to reset '
          write(6,*) ' nroot, nvec ', nroot,nvec
        END IF
C
        IF(IPRT.GE.10000) THEN
          WRITE(6,*) ' Initial vectors on LU1'
          CALL REWINO(LU1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
          WRITE(6,*) ' Initial vectors on LU3'
          CALL REWINO(LU3)
          DO IVEC = 1, NVEC3
             CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF
*
*. First orthonormalization, based on subspace matrices
*
C       GET_CNEWCOLD_BAS(CN,CNO,NVEC,NROOT,SCR,NVECUT,THRES)
        THRES = 1.0D-5
        CALL GET_CNEWCOLD_BAS(AVEC,XJEP,NVEC,NROOT,WORK2,NVECUT,
     &       THRES)
*
*. 3.1 : Current wave function approximations, collect on LU7
*
*
* Novel reset (Feb. 2013)
*
*. Obtain the corrections and their sigma vectors
*
C            CNEW_AS_COLD_PLUS_CORR(CNEW,NROOT,NVEC,DELTA)
        CALL CNEW_AS_COLD_PLUS_CORR(AVEC,NROOT,NVEC,WORK)
*
*. Work now contains expansion of NROOT corrections, assemble these 
*
*
        LU_FILE1(1) = LU1
        LU_FILE1(2) = LU3
        NVC_FILE1(1) = NROOT
        NVC_FILE1(2) = NVEC3
C            TRAVCMF(NFILE,LU_FILE,NVC_FILE,NIN,NOUT,LUOUT,X,
C    &           LUSC1, LUSC2,VEC1,VEC2)
        CALL TRAVCMF(2,LU_FILE1,NVC_FILE1,NVEC,NROOT,LU5,WORK,
     &       LU6,LU7,VEC1,VEC2)
C            COPNVCD(LUIN,LUOUT,NVEC,SEGMNT,IREW,LBLK)
        CALL COPNVCD(LU5,LU3,NROOT,VEC1,1,LBLK)
*
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' Corrections on LU3 '
          CALL REWINO(LU3)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF
*
*. 3.2 : corresponding sigma vectors
*
        LU_FILE1(1) = LU2
        LU_FILE1(2) = LU4
        NVC_FILE1(1) = NROOT
        NVC_FILE1(2) = NVEC-NROOT
C            TRAVCMF(NFILE,LU_FILE,NVC_FILE,NIN,NOUT,LUOUT,X,
C    &           LUSC1, LUSC2,VEC1,VEC2)
        CALL TRAVCMF(2,LU_FILE1,NVC_FILE1,NVEC,NROOT,LU5,WORK,
     &       LU6,LU7,VEC1,VEC2)
C            COPNVCD(LUIN,LUOUT,NVEC,SEGMNT,IREW,LBLK)
        CALL COPNVCD(LU5,LU4,NROOT,VEC1,1,LBLK)
*
        WRITE(6,*) ' Correction vectors reset '
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' Updated correction vectors '
          WRITE(6,*) ' =========================== '
          WRITE(6,*)
          CALL REWINO(LU4)
          DO JVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU4,0,LBLK)
             WRITE(6,*)
          END DO
          WRITE(6,*) ' Updated Sigma for correction vectors '
          WRITE(6,*) ' ==================================== '
          CALL REWINO(LU4)
          DO JVEC = 1, NROOT
             CALL WRTVCD(VEC2,LU4,0,LBLK)
             WRITE(6,*)
          END DO
         END IF  ! Print
*
*. Obtain the new wave function approximations and their sigma vectors on LU1, LU2, respectively
*
        CALL REWINO(LU1)
        CALL REWINO(LU3)
        CALL REWINO(LU7)
        DO IROOT = 1, NROOT
C              VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
          CALL VECSMD(VEC1,VEC2,ONE,ONE,LU1,LU3,LU7,0,LBLK)
        END DO
        CALL REWINO(LU1)
        CALL REWINO(LU7)
        DO IROOT = 1, NROOT
          CALL COPVCD(LU7,LU1,VEC1,0,LBLK)
        END DO
*
        CALL REWINO(LU2)
        CALL REWINO(LU4)
        CALL REWINO(LU7)
        DO IROOT = 1, NROOT
C              VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
          CALL VECSMD(VEC1,VEC2,ONE,ONE,LU2,LU4,LU7,0,LBLK)
        END DO
        CALL REWINO(LU2)
        CALL REWINO(LU7)
        DO IROOT = 1, NROOT
          CALL COPVCD(LU7,LU2,VEC1,0,LBLK)
        END DO
        WRITE(6,*) '  Solution vectors reset '
*
* We now have the corrections vectors on LU2/LU4 and the current wave-function
* approximations on LU1/LU3. Orthonormalize the correction vectors
*

        IF(3*NROOT.LE.MAXVEC) THEN
          WRITE(6,*) ' 3*NROOT < MAXVEC section entered '
*
* We have on LU1/LU2 : NROOT orthonormal eigenvectors
*            LU3/LU4: NADD1 auxilary vectors that are not orthogonal to each other or to
*                     the vectors on LU1/LU2
*
*. Orthogonalize the vectors on LU3 to those at LU1
*
C     ORTVECTORS_TO_ORTNVECTORS(LUVEC1,NVEC1,LUVEC2,NVEC2,
C    &           THRES,VEC1,VEC2,LUVEC3,NVEC3,LUSC1,LUSC2,LUSC3,X,XC)
* (LUSC1 or LUSC2 may be = to LUOUT)
          THRES2 = 1.0D-3
          CALL ORTVECTORS_TO_ORTNVECTORS(LU3,NROOT,LU1,NROOT,
     &         THRES2,VEC1,VEC2,LU5,NVECEFF,LU5,LU6,LU7,WORK,
     &         WORK(1+NVEC**2),IXJEP)
          IF(NVECEFF.NE.NROOT) THEN
            WRITE(6,*) 
     &      ' Warning: Vector eliminated in subspace construction'
          END IF
          NADD = NVECEFF
          NADD2 = NVECEFF
          WRITE(6,*) ' NADD2 = ', NADD2
*. And place orthogonalized vectors on LU3
          CALL COPNVCD(LU5,LU3,NADD2,VEC1,1,LBLK)
C              COPNVCD(LUIN,LUOUT,NVEC,SEGMNT,IREW,LBLK)
*. Check that vectors on LU3 and LU1 are orthogonal
          WRITE(6,*) ' Overlap <roots!additional vectors(1)>'
          CALL GET_OVLAP_FOR_VECTORS(LU1,NROOT,LU3,NADD2,LU7,
     &          XJEP,VEC1,VEC2)
         XMAX = XMNMX(XJEP,NROOT*NADD2,3)
         WRITE(6,*) ' Largest element of <roots!add vectors(1)> ', XMAX
*
*. The sigma vectors corresponding to the orthogonalized corrections
*
*. The part in the ROOT basis on LU5
C     TRAVCD(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
          CALL REWINO(LU5)
          CALL TRAVCD(VEC1,VEC2,WORK,NROOT,NADD2,LU2,LU5,0,
     &                LBLK,LU6,LU7)
          WRITE(6,*) ' Sigma vectors of orth. corrections, (1) finito'
*. And add to LU4 and save on LU6
          CALL REWINO(LU5)
          CALL REWINO(LU4)
          CALL REWINO(LU6)
          IVEC3 = 1
          NWRIT = 0
          DO IVEC = 1, NROOT
            IF(IXJEP(IVEC3).EQ.IVEC) THEN
C                  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
              CALL VECSMD(VEC1,VEC2,ONE,ONE,LU5,LU4,LU6,0,LBLK)
              WRITE(6,*) ' vector used, IVEC, IVEC3 = ', IVEC, IVEC3
              NWRIT = NWRIT + 1
              IF(IVEC3.LT.NROOT) IVEC3 = IVEC3 +1
            ELSE
C                  SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
              CALL SKPVCD(LU4,1,VEC1,0,LBLK)
              WRITE(6,*) ' vector skipped, IVEC, IVEC3 = ', IVEC, IVEC3
            END IF
          END DO
          WRITE(6,*) ' Number of written (a) ', NWRIT
          CALL COPNVCD(LU6,LU4,NADD2,VEC1,1,LBLK)
          WRITE(6,*) ' Sigma vectors constructed for ort to roots '
*
*. Orthonormalize the correction vectors on LU3
*
*. Overlap before MGS, for check
*
          THRES1 = 1.D-3
          CALL MGSCD(LU3,NADD,LU3,NVECUT,WORK,VEC1,VEC2,
     &         LU5,LU6,LU7,THRES1,WORK(NVEC**2+1))
C     MGSCD(LUIN,NVECIN,LUOUT,NVECUT,X,VEC1,VEC2,
C    &                LUSCR1,LUSCR2,LUSCR3,THRES,XC)
          WRITE(6,*) ' Home from MGSCD '
          IF(NVECUT.NE.NADD) THEN
           WRITE(6,*) ' Number of vectors reduced '
          END IF
          NADDO = NADD
          NADD = NVECUT
*. Test overlap of new vectors
          CALL GET_OVLAP_FOR_VECTORS(LU3,NADD,LU3,NADD,LU7,
     &         XJEP,VEC1,VEC2)
C          GET_OVLAP_FOR_VECTORS(LU1,N1,LU2,N2,LUSCR,OVLAP,
C    &                                VEC1,VEC2)
          CALL CHECK_UNIT_MAT(XJEP,NADD,XMAX_DIA,XMAX_OFD,1)
*. And overlap with original vectors
          WRITE(6,*) ' Overlap <roots!additional vectors(2)>'
          CALL GET_OVLAP_FOR_VECTORS(LU1,NROOT,LU3,NADD,LU7,
     &          XJEP,VEC1,VEC2)
         XMAX = XMNMX(XJEP,NROOT*NADD,3)
         WRITE(6,*) ' Largest element of <roots!add vectors(2)> ', XMAX
*
* And obtain the corresponding Sigma vectors
*
C     TRAVCD(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
          CALL REWINO(LU6)
          CALL TRAVCD(VEC1,VEC2,WORK,NADDO,NADD,LU4,LU6,
     &           1,LBLK,LU5,LU7)
*
          NNVEC = NROOT + NADD
        END IF
*       /\ End if more than NROOT vectors to be reset
        NVEC = NNVEC 
        WRITE(6,*) ' End of new reset '
*
      IF(IPRT.GE.1000) THEN
        WRITE(6,*) ' Updated eigenvector approximations '
        WRITE(6,*) ' ================================== '
        WRITE(6,*)
        CALL REWINO(LU1)
        DO JVEC = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
           WRITE(6,*)
        END DO
*
        WRITE(6,*) ' Updated Sigma for eigenvector approximations '
        WRITE(6,*) ' ============================================ '
        CALL REWINO(LU2)
        DO JVEC = 1, NADD
           CALL WRTVCD(VEC2,LU2,0,LBLK)
           WRITE(6,*)
        END DO
*
        WRITE(6,*) ' Updated Sigma for Additional vectors '
        WRITE(6,*) ' ==================================== '
        CALL REWINO(LU4)
        DO JVEC = 1, NADD
           CALL WRTVCD(VEC2,LU4,0,LBLK)
           WRITE(6,*)
        END DO
       END IF  ! Print
*
* We now have:
* NROOT current approximations to the eigenvectors and their sigmavectors on LU1 and LU3, respectively
* NVEC-NROOT additional vectors and their sigmavectors on LU2 and LU4, respectively
*
*. New subspace
*
*. Calculate subspace Hamiltonian from actual vectors 
*. on disc
        I_NEW_OR_OLD = 1
*
        IF(I_NEW_OR_OLD.EQ.2) THEN
        IF(IPRT.GE.1000) write(6,*) ' Subspace hamiltonian' 
        CALL REWINO(LU1)
        CALL REWINO(LU3)
        DO IVEC = 1, NVEC
*
         CALL REWINO(LU5)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINO(LU2)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
         END DO
*
         CALL REWINO(LU4)
         DO JVEC = NROOT+1,IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
         END DO
        END DO ! Loop over IVEC
        IF(IPRT.GE.10) THEN
          write(6,*) ' Reset hamiltonian'
          call prsym(aproj,nvec)
        END IF
*. Test : Orthogonality of new vectors
        CALL REWINO(LU1)
        CALL REWINO(LU3)
        DO IVEC = 1, NVEC
*
         CALL REWINO(LU5)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINO(LU1)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU1,0,LBLK)
         END DO
         CALL REWINO(LU3)
         DO JVEC = NROOT+1,IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU3,0,LBLK)
         END DO
        END DO ! Loop over IVEC
        IF(IPRT.GE.10) THEN
         write(6,*) ' Overlap matrix    '
         call prsym(xjep,nvec)
        END IF
*
        ELSE
*
*. Alternative calculation of subspace Hamiltonian '
*
        NVC_FILE1(1) = NROOT
        LU_FILE1(1) = LU1
        NVC_FILE1(2) = NADD
        LU_FILE1(2) = LU3
*
        NVC_FILE2(1) = NROOT
        LU_FILE2(1) = LU2
        NVC_FILE2(2) = NADD
        LU_FILE2(2) = LU4
*
        CALL SUBSPC_MF(2,LU_FILE1,NVC_FILE1,
     &                 2,LU_FILE2,NVC_FILE2,
     &                 VEC1,VEC2,APROJ,1,LU6)
        IF(IPRT.GE.10) THEN
          write(6,*) ' Reset hamiltonian (new route)'
          call prsym(aproj,nvec)
        END IF
*. Alternative calculation of subspace Overlap '
        NVC_FILE1(1) = NROOT
        LU_FILE1(1) = LU1
        NVC_FILE1(2) = NADD
        LU_FILE1(2) = LU3
        CALL SUBSPC_MF(2,LU_FILE1,NVC_FILE1,
     &                 2,LU_FILE1,NVC_FILE1,
     &                 VEC1,VEC2,XJEP,1,LU6)
C      SUBSPC_MF(NFL_L,LUFL_L,NVCFL_L,NFL_R,LUFL_R,NVCFL_R,
C    &                     VEC1,VEC2,SUBMAT,ISYM,LUSCR)
        IF(IPRT.GE.10) THEN
          write(6,*) ' Overlap matrix (alt. route) '
          call prsym(xjep,nvec)
        END IF
        CALL CHECK_UNIT_MAT(XJEP,NVEC,XMAX_DIA,XMAX_OFD,1)
        END IF ! Switch between reset types
*. End of resetting business
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000 ! End of loop over iterations
 1001 CONTINUE
*
* ( End of loop over iterations )
*
C?    WRITE(6,*) ' ITER, MAXIT = ', ITER, MAXIT
      IF(CONVER.AND.ITER.LE.MAXIT) THEN
        ITER_LAST = ITER-1
      ELSE
        ITER_LAST = MAXIT
      END IF
C?    WRITE(6,*) ' ITER, MAXIT, ITER_LAST, CONVER = ',
C?   &             ITER, MAXIT, ITER_LAST, CONVER
      DO 1601 IROOT = 1, NROOT
         FINEIG(IROOT) = EIG(ITER_LAST,IROOT)+EIGSHF
         RNRM_CNV(IROOT) = RNRM(ITER_LAST,IROOT)
C?    WRITE(6,*) ' ITER, MAXIT, ITER_LAST, RNRM_CNV(1) = ', 
C?   &             ITER, MAXIT, ITER_LAST, RNRM_CNV(1)
 1601 CONTINUE
C?    WRITE(6,*) ' MICDV6(still), FINEIG(1) = ', FINEIG(1)
 
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
*
      DO IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(ITER,IROOT) + EIGSHF
      END DO
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
C?    WRITE(6,*) ' MICDV6(c), FINEIG(1) = ', FINEIG(1)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.1000) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
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
      ELSE
CM      DO 1601 IROOT = 1, NROOT
CM         FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
CM         RNRM_CNV(IROOT) = RNRM(ITER,IROOT)
CM1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
c $Id$
