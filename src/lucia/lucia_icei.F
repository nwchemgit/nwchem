*
* Note on the transformation
* For a given type of external state and a given symmetry of 
* internal states, the transformation from the elementary EI
* operators to the orthonormal basis goes as
*
* 1: a orthogonal transformation X1 that diagonalizes metric is 
* set up: X1(T) S X1 = Sigma, X1(T) X(1) = 1, dimension: LINT*LINT
* 2: The nonsingular part of X1 is scaled so it gives eigenvectors 
*    of unit norm: X1s = X1 * sigma(-1/2), X1s(T) S X1S = 1, dimension:
*    LINT*LORTN
* 3: The zeroorder Hamiltonian is set up in the X1s basis and
* diagonalized by a orthogonal matrix:
* H0(tilde) = X1s(T) H0 X1s, X2(T) H0(tilde) X2 = epsilon
*
* The transformation form the elementary basis to the orthonormal 
* basis reads therefore X = X1 sigma(-1/2) X2, and its inverse is 
* X(-1) = X2(T) sigma(1/2) X1(T). To obtain inverse X(-1) readily
* the matrices X1, X2 and the vector sigma are stored.
*
* Missing: Calculate diagonal
* <0!0(INT_J)O(EXT(J)! H0!O(EXT_I)O(INT_I)|0>
*. Assuming that H0 = H0(INT) + H0(EXT) gives
* <0!0(INT_J)O(EXT(J)! H0!O(EXT_I)O(INT_I)|0> =
* 
      SUBROUTINE PROJ_TO_NONRED(VECIN,VECOUT,ITSYM,VECSCR)
*
*. Project vector to nonredundant space
*
*. VECOUT = X1 Sigma^-1 X1^T S VECIN
*
*. Jeppe Olsen, Oct 17, 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'lucinp.inc'
*
      REAL*8
     &INPROD
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT(*)
*. Scratch
      DIMENSION VECSCR(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'PROJNR')
      CALL LUCIAQENTER('PRJNR')
*
      NTEST = 100
     
      CALL PROJ_TO_NONRED_SLAVE(VECIN,VECOUT,ITSYM,VECSCR,
     &WORK(KL_X1_INT_EI_FOR_SE),WORK(KL_SG_INT_EI_FOR_SE),
     &WORK(KL_S_INT_EI_FOR_SE),
     &WORK(KL_NDIM_EX_ST),WORK(KL_NDIM_IN_SE),
     &WORK(KL_N_ORTN_FOR_SE),
     &WORK(KL_IREO_EI_ST),
     &N_EXTOP_TP,NSMOB,N_ZERO_EI,NDIM_EI,I_INCLUDE_UNI)
*
      XNRM_IN = SQRT(INPROD(VECIN,VECIN,NDIM_EI-1))
      XNRM_OUT = SQRT(INPROD(VECOUT,VECOUT,NDIM_EI-1))
      XNRM_DIFF = XNORM_DIFF(VECIN,VECOUT,NDIM_EI-1)
*
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Input and output vectors from PROJ_TO_NONRED'
       CALL WRT_2VEC(VECIN,VECOUT,NDIM_EI)
      END IF
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' PROJ_TO_NONRED: XNRM_IN, XNRM_OUT,XNRM_DIFF =',
     &                               XNRM_IN, XNRM_OUT,XNRM_DIFF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'PROJNR')
      CALL LUCIAQEXIT('PRJNR')
     
*
      RETURN
      END
      FUNCTION XNORM_DIFF(VEC1,VEC2,NDIM)
*
* Norm of difference of two vectors
*
*. Jeppe Olsen, Oct. 18, 2009
*
      INCLUDE 'implicit.inc'
      DIMENSION VEC1(NDIM),VEC2(*)
*
      XNORM = 0.0D0
      DO I = 1, NDIM
       XNORM = XNORM + (VEC1(I)-VEC2(I))**2
      END DO
      XNORM = SQRT(XNORM)
*
      XNORM_DIFF = XNORM
*
      RETURN
      END
*
      SUBROUTINE PROJ_TO_NONRED_SLAVE(VECIN,VECOUT,ITSYM,VECSCR,
     &X1,SG,S,
     &NDIM_EX_ST,NDIM_IN_ST,NDIM_ORT_ST,
     &IREO_EI_ST,
     &N_EXTP,NSMOB,N_ZERO_EI,NDIM_EI,I_INCLUDE_UNI)
*
*. Project a vector to nonredundant basis
*
*    Vecout = X1 Sigma^-1 X1^T S Vecin
*
* Input and output vectors in CAAB order
*
*. Jeppe Olsen, Oct18, 2009
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*.General input
*
      DIMENSION X1(*),SG(*),S(*)
      INTEGER IREO_EI_ST(*)
*
      DIMENSION NDIM_EX_ST(NSMOB,N_EXTP),
     &          NDIM_ORT_ST(NSMOB,N_EXTP),
     &          NDIM_IN_ST(NSMOB,N_EXTP)
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT(*)
*. Scratch
      DIMENSION VECSCR(*) 
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'PROJNS')
*.1: Reorder to EI order:  CAAB(in VECIN) => EI-order(in VECSCR) 
      CALL SGATVEC(VECSCR,VECIN,IREO_EI_ST,NDIM_EI)
      IF(I_INCLUDE_UNI.EQ.1) THEN 
        VECSCR(NDIM_EI+1) = VECIN(NDIM_EI+1)
      END IF
*.2 Perform transformation, save in VECSCR
      IOFF_SG=1
      IOFF_S=1
      IOFF_X1=1
      IOFF_ORT=1
      IOFF_VEC = 1
*
      DO I_EXTP = 1, N_EXTP
       DO I_EXSM = 1, NSMOB
        I_INSM = MULTD2H(I_EXSM,ITSYM)
*
        N_EX = NDIM_EX_ST(I_EXSM,I_EXTP)
        N_IN = NDIM_IN_ST(I_INSM,I_EXTP)
        N_ORT = NDIM_ORT_ST(I_INSM,I_EXTP)
*
        CALL  PROJ_EI_BL(X1(IOFF_X1),S(IOFF_S),SG(IOFF_SG),
     &        VECSCR(IOFF_VEC),VECOUT(IOFF_VEC),N_EX,N_IN,N_ORT) 
        CALL COPVEC(VECOUT(IOFF_VEC),VECSCR(IOFF_VEC),N_EX*N_IN)
*
        IOFF_SG = IOFF_SG + N_ORT
        IOFF_S = IOFF_S + N_IN*(N_IN+1)/2
        IOFF_X1 = IOFF_X1 + N_ORT*N_IN
        IOFF_VEC = IOFF_VEC + N_EX*N_IN
*
       END DO
      END DO
*. Reorder to CAAB order
      CALL SSCAVEC(VECOUT,VECSCR,IREO_EI_ST,NDIM_EI)
      IF(I_INCLUDE_UNI.EQ.1) THEN
        VECOUT(NDIM_EI+1) = VECSCR(NDIM_EI+1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'PROJNS')
      RETURN
      END
      SUBROUTINE PROJ_EI_BL(X1,S,SG,VECIN,VECOUT,N_EX,N_IN,N_ORT)
*
* Project a block of a vector to non-redundant basis
* Input vector is destroyed in the process
*
*    Vecout = X1 Sigma^-1 X1^T S Vecin
*
*. Jeppe Olsen, oct. 18, 2009
*
      INCLUDE 'implicit.inc'
*. General input
      DIMENSION X1(N_ORT,N_IN),S(N_IN*(N_IN+1)/2),SG(N_ORT)
*. input and scratch
C     DIMENSION VECIN(N_IN,N_EX)
      DIMENSION VECIN(N_IN*N_EX)
*. output
      DIMENSION VECOUT(N_IN,N_EX)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Input vector to PROJ_EI_BL'
       CALL WRTMAT(VECIN,N_IN,N_EX,N_IN,N_EX)
       WRITE(6,*) 'N_IN, N_ORT, N_EX = ', N_IN, N_ORT, N_EX
      END IF
*1: S Vecin in Vecout
      DO I_EX = 1, N_EX
C       CALL SYMTVC(S,VECIN(1,I_EX),VECOUT(1,I_EX),N_IN)
        CALL SYMTVC(S,VECIN((I_EX-1)*N_IN+1),VECOUT(1,I_EX),N_IN)
C       SYMTVC(A,VECIN,VECOUT,NDIM)
      END DO

*2:  X1^T S Vecin in Vecin
      FACTOR_AB = 1.0D0
      FACTOR_C  = 0.0D0
*
      CALL MATML7(VECIN,X1,VECOUT,N_ORT,N_EX,N_IN,N_ORT,N_IN,N_EX,
     &            FACTOR_C,FACTOR_AB,1)
*3:  Sigma^-1 X1^T S Vecin in Vecin
      DO J = 1, N_EX
       DO I = 1, N_ORT
C        VECIN(I,J) = VECIN(I,J)/SG(I)
         VECIN((J-1)*N_ORT+I) = VECIN((J-1)*N_ORT+I)/SG(I)
       END  DO
      END  DO
*4:  X1 Sigma^-1 X1^T S Vecin in Vecout
      CALL MATML7(VECOUT,X1,VECIN,N_IN,N_EX,N_IN,N_ORT,N_ORT,N_EX,
     &            FACTOR_C,FACTOR_AB,0)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output block from PROJ_EI_BL'
       CALL WRTMAT(VECOUT,N_IN,N_EX,N_IN,N_EX)
      END IF
*
      RETURN
      END
      SUBROUTINE SYMMAT_TV(SYMMAT,VECIN,VECOUT,NDIM)
*
* A symmetric matrix is stored rowwise, lowerhalf
* multiply with vector VECIN
* 
*. NO cache-considerations have been invoked
*
*. Jeppe Olsen, oct. 18, 2009
*
      INCLUDE 'implicit.inc'
*.input
      DIMENSION SYMMAT(*)
      DIMENSION VECIN(*)
*.Output
      DIMENSION VECOUT(*)
*
      ZERO = 0.0D0
      CALL SETVEC(VECOUT,ZERO,NDIM)
      IJ = 0
      DO  I = 1, NDIM
        DO J = 1, I
          IJ = IJ + 1
          VECOUT(J) = VECOUT(J) + SYMMAT(IJ)*VECIN(I)
          VECOUT(I) = VECOUT(I) + SYMMAT(IJ)*VECIN(J)
        END DO
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output vector from SYMMAT_TV'
       CALL WRTMAT(VECOUT,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE MODDIAG(H0DIAG,NDIM,XMIN)
*
*. Jacobian Diagonal of H0 is given.
*. Check this and replace all values smaller than XMIN
* by XMIN
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
      DIMENSION H0DIAG(NDIM)
*
      XMIN_FOUND = ABS(H0DIAG(1))
      NSHIFT = 0
*
      DO I = 1, NDIM
        IF(ABS(H0DIAG(I)).LT.XMIN_FOUND) THEN
          XMIN_FOUND = H0DIAG(I)
        END IF
        IF(ABS(H0DIAG(I)).LT.XMIN) THEN
          H0DIAG(I) = XMIN
          NSHIFT = NSHIFT + 1
        END IF
      END DO
*
      NTEST = 10
      IF(NTEST.GE.1) THEN
       WRITE(6,*) ' Check of J(H0DIAG)'
       WRITE(6,*) ' Imposed lower value =', XMIN
       WRITE(6,*) ' Number of elements shifted =', NSHIFT
       WRITE(6,*) ' Lowest value found =', XMIN_FOUND
      END IF
*
      RETURN
      END
      SUBROUTINE TEST_GENMAT(MSTV,NVAR,NVAR_INT,I_DO_DIAG)
*
* A method using a matrix vector routine MSTV is 
* tested by constructing complete matrix and metric and maybe
* diagonalizing (If I_DO_DIAG = 1)
*
*. NVAT_INT: Sometimes a different number of variables
*. are used internally in matrix vector. The vectors
*. should have this dimension
*
* Jeppe Olsen, Hotel room in Dusseldorf, 2009
*
      INCLUDE 'wrkspc.inc'
      EXTERNAL MSTV
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'TST_GN')
*
      WRITE(6,*) ' NVAR, NVAR_INT =',NVAR,NVAR_INT
      CALL MEMMAN(KLMATH,NVAR**2 ,'ADDL  ',2,'MATH  ')
      CALL MEMMAN(KLMATS,NVAR**2 ,'ADDL  ',2,'MATS  ')
      CALL MEMMAN(KLVEC1,NVAR_INT,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,NVAR_INT,'ADDL  ',2,'VEC2  ')
      CALL MEMMAN(KLVEC3,NVAR_INT,'ADDL  ',2,'VEC3  ')
*
*. Space for diagonalization(I use a simple routine for this)
      IF(I_DO_DIAG.EQ.1) THEN
        CALL MEMMAN(KLSCRM1,NVAR**2,'ADDL  ',2,'MAT1  ')
        CALL MEMMAN(KLSCRM2,NVAR**2,'ADDL  ',2,'MAT2  ')
      ELSE
        KLSCRM1 = 1
        KLSCRM2 = 1
      END IF
*
      CALL TEST_GENMAT_INNER(MSTV,NVAR,I_DO_DIAG,
     &                       WORK(KLMATH),WORK(KLMATS),
     &                       WORK(KLVEC1),WORK(KLVEC2),WORK(KLVEC3),
     &                       WORK(KLSCRM1),WORK(KLSCRM2))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'TST_GN')
*
      RETURN
      END
      SUBROUTINE TEST_GENMAT_INNER(MSTV,NVAR,I_DO_DIAG,
     &            HMAT,SMAT,VEC1,VEC2,VEC3,SCRM1,SCRM2)
*
* A method using a matrix vector routine MSTV is 
* tested by constructing complete matrix and metric and maybe
* diagonalizing - Inner routine ( Well, sounds better than slave
* routine)
*
* Jeppe Olsen, Hotel room in Dusseldorf, 2009
*
      INCLUDE 'wrkspc.inc'
      EXTERNAL MSTV
      DIMENSION HMAT(NVAR,NVAR),SMAT(NVAR,NVAR)
      DIMENSION VEC1(NVAR),VEC2(NVAR),VEC3(NVAR)
      DIMENSION SCRM1(*),SCRM2(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'TS_GNI')
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' TEST_GENMAT_INNER reporting '
      END IF
*. Test: just parameter 136
C?    ZERO = 0.0D0
C?    CALL SETVEC(VEC1,ZERO,NVAR)
C?    VEC1(136) = 1.0D0
C?    CALL MSTV(VEC1,VEC2,VEC3,1,1)
C?    STOP 'Enforced stop after call to test MSTV'
*
      DO I = 1, NVAR
       ZERO = 0.0D0
       CALL SETVEC(VEC1,ZERO,NVAR)
       ONE = 1.0D0
       VEC1(I) = ONE
       WRITE(6,*) ' Constructing HV, SV for element = ', I
       WRITE(6,*) ' --------------------------------------'
       CALL MSTV(VEC1,VEC2,VEC3,1,1)
       CALL COPVEC(VEC2,HMAT(1,I),NVAR)
       CALL COPVEC(VEC3,SMAT(1,I),NVAR)
      END DO
*. Compare to unit matrix
      CALL COMPARE_TO_UNI(SMAT,NVAR)
C     SUBROUTINE COMPARE_TO_UNI(A,NDIM)
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' H-matrix: '
       CALL WRTMAT(HMAT,NVAR,NVAR,NVAR,NVAR)
       WRITE(6,*) ' S-matrix: '
       CALL WRTMAT(SMAT,NVAR,NVAR,NVAR,NVAR)
      END IF
*
      IF (I_DO_DIAG.EQ.1) THEN
*. Diagonalize 
C       GENDIA(HIN,SIN,VOUT,EIGENV,PVEC,NDIM,ISORT)
        CALL GENDIA(HMAT,SMAT,SCRM1,VEC1,SCRM2,NVAR,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'TS_GNI')
      RETURN
      END
*
      SUBROUTINE PRINT_ZERO_TRMAT
*
* Print transformation matrices for obtaining orthogonal
* internal zero-order states
*
*. Jeppe Olsen, the train to dusseldorf, aug14, 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cei.inc'
*
      WRITE(6,*) ' ----------------------------------------------'
      WRITE(6,*) ' Transformation matrices for orthonormal states'
      WRITE(6,*) ' ----------------------------------------------'
      WRITE(6,*) ' PRINT_ZERO.(a): NSMOB, N_EXTOP_TP=',NSMOB,N_EXTOP_TP

      CALL PRINT_ZERO_TRMAT_SLAVE(NSMOB,N_EXTOP_TP,
     &     WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &     WORK(KL_X1_INT_EI_FOR_SE),WORK(KL_X2_INT_EI_FOR_SE),
     &     WORK(KL_SG_INT_EI_FOR_SE),
     &     WORK(KL_IBX1_INT_EI_FOR_SE),WORK(KL_IBX2_INT_EI_FOR_SE),
     &     WORK(KL_IBSG_INT_EI_FOR_SE))
*
      RETURN
      END
      SUBROUTINE PRINT_ZERO_TRMAT_SLAVE(NSMOB,N_EXTP,
     &            N_ORTN_FOR_SE,N_INT_FOR_SE,
     &            X1,X2,SG,
     &            IBX1_INT_EI_FOR_SE,IBX2_INT_EI_FOR_SE,
     &            IBSG_INT_EI_FOR_SE)
*. Print matrices relevant for generation of orthonormal
*. zero-order stated
*
*. Jeppe Olsen, Aug.14, 2009
*
      INCLUDE 'implicit.inc'
      DIMENSION X1(*),X2(*),SG(*)
*
      INTEGER N_ORTN_FOR_SE(NSMOB,N_EXTP),N_INT_FOR_SE(NSMOB,N_EXTP)
*
      INTEGER IBX1_INT_EI_FOR_SE(NSMOB,N_EXTP)
      INTEGER IBX2_INT_EI_FOR_SE(NSMOB,N_EXTP)
      INTEGER IBSG_INT_EI_FOR_SE(NSMOB,N_EXTP)
*
      NTEST = 0
      WRITE(6,*) ' PRINT_ZERO.. : N_EXTP, NSMOB = ', N_EXTP, NSMOB
      DO I_EXTP = 1, N_EXTP
       DO I_EXSM = 1, NSMOB
        WRITE(6,*) ' Info for external type and sym =', I_EXTP,I_EXSM
        N_INT = N_INT_FOR_SE(I_EXSM,I_EXTP)
        N_ORTN= N_ORTN_FOR_SE(I_EXSM,I_EXTP)
*
        IOFF_X1 = IBX1_INT_EI_FOR_SE(I_EXSM,I_EXTP)
        IOFF_X2 = IBX2_INT_EI_FOR_SE(I_EXSM,I_EXTP)
        IOFF_SG = IBSG_INT_EI_FOR_SE(I_EXSM,I_EXTP)
*. Test output for printing routines - I must be getting old
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' N_INT, N_ORTN = ',  N_INT, N_ORTN
          WRITE(6,*) ' IOFF_X1, IOFF_X2 ,IOFF_SG =',
     &                 IOFF_X1, IOFF_X2 ,IOFF_SG 
        END IF
*
        WRITE(6,*) ' Block of X1, X2, SG '
        CALL WRTMAT(X1(IOFF_X1),N_INT,N_ORTN,N_INT,N_ORTN)
        WRITE(6,*)
        CALL WRTMAT(X2(IOFF_X2),N_ORTN,N_ORTN,N_ORTN,N_ORTN)
        WRITE(6,*)
        CALL WRTMAT(SG(IOFF_SG),1,N_ORTN,1,N_ORTN)
        WRITE(6,*)
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE TRANS_CAAB_ORTN(T_CAAB,T_ORTN,ITSYM,ICO,ILR,SCR,
     &                          ICOCON)
*
*. Transform a vector between standard CAAB order and EI order 
*. with orthornormal internal states
*
* ICO = 1: CAAB => Ortn
* ICO = 2: Ortn => CAAB
*
* ICOCON = 1 => Covariant transformation
* ICOCON = 2 => Contravariant transformation
*
*. Jeppe Olsen, July 29, 2009, sidder på Rønnevej 12 med Jette
*                              gættende kryds og tværs...
*
*. NOTE: at the moment no signs are used. When transforming tp
* spin-adapted operators, sign changes going from CAAB tp EI form 
* must be considered.
*
      INCLUDE 'wrkspc.inc'
      INCLUDE  'cei.inc'
*. Input/output
      DIMENSION T_ORTN(*),T_CAAB(*)
*. Scratch: Dimension of full CAAB expansion
      DIMENSION SCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) WRITE(6,*) ' Starting with transformation'
*
*
      IF(ICO.EQ.1) THEN
*. CAAB => Ortn is done in two steps 
*  CAAB(in T_CAAB) => EI-order(in SCR) => Ortn(in T_ORTN)
       CALL SGATVEC(SCR,T_CAAB,WORK(KL_IREO_EI_ST),NDIM_EI)
       IF(I_INCLUDE_UNI.EQ.1) THEN 
         SCR(NDIM_EI+1) = T_CAAB(NDIM_EI+1)
       END IF
       CALL  TRANS_EI_ORTN(SCR,T_ORTN,ITSYM,1,ILR,ICOCON)
      ELSE
*. Ortn => CAAB is done in two steps
*. Ortn(in T_ORTN) => EI-order(in SCR) => CAAB (in T_CAAB)
       CALL TRANS_EI_ORTN(SCR,T_ORTN,ITSYM,2,ILR,ICOCON)
       CALL SSCAVEC(T_CAAB,SCR,WORK(KL_IREO_EI_ST),NDIM_EI)
*
       IF(I_INCLUDE_UNI.EQ.1) THEN
         T_CAAB(NDIM_EI+1) = SCR(NDIM_EI+1)
       END IF
      END IF
*
      IF(NTEST.GE.100) THEN
        IF(ICO.EQ.1) THEN
          WRITE(6,*) ' CAAB => ORTN transformation '
        ELSE
          WRITE(6,*) ' ORTN => CAAB transformation '
        END IF
        IF(ICOCON.EQ.1) THEN
          WRITE(6,*) ' Covariant transformation'
        ELSE
          WRITE(6,*) ' Contravariant transformation'
        END IF
*
        WRITE(6,*) ' The T vector in zero-order EI basis: '
        CALL PRINT_T_EI(T_ORTN,2,ITSYM)
COLD    IF(I_INCLUDE_UNI.EQ.1) THEN
COLD      WRITE(6,*) ' Element corresponding to unit operator',
COLD &    T_ORTN(N_ZERO_EI+1)
COLD    END IF
        WRITE(6,*) ' The T vector in CAAB form'
        CALL WRTMAT(T_CAAB,1,NDIM_EI,1,NDIM_EI)
        IF(I_INCLUDE_UNI.EQ.1) THEN
          WRITE(6,*) ' Element corresponding to unit operator',
     &    T_CAAB(NDIM_EI+1)
        END IF
*
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' Finished with transformation'
*
      RETURN
      END 
      SUBROUTINE PRINT_T_EI(T,IEO,ITSYM)
*
* Print a T(I,E) vector given in elementary(IEO=1) or 
* orthonormal(IEO=2) form
*
*. Jeppe Olsen, July 29, 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'lucinp.inc'
*. Specific input
      DIMENSION T(*)
*
      WRITE(6,*)
      WRITE(6,*) ' ------------------------------------'
      WRITE(6,*) ' T(I,E) vector with symmetry ', ITSYM
      WRITE(6,*) ' ------------------------------------'
      WRITE(6,*)
      IF(IEO.EQ.1) THEN
        CALL PRINT_T_EI_SLAVE(T,ITSYM,N_EXTOP_TP,
     &       WORK(KL_NDIM_EX_ST),WORK(KL_NDIM_IN_SE),NSMOB)
        IF(I_INCLUDE_UNI.EQ.1) THEN
          WRITE(6,*) ' Element corresponding to unit-operator',
     &    T(NDIM_EI+1)
        END IF
      ELSE
        CALL PRINT_T_EI_SLAVE(T,ITSYM,N_EXTOP_TP,
     &       WORK(KL_NDIM_EX_ST),WORK(KL_N_ORTN_FOR_SE),NSMOB)
        IF(I_INCLUDE_UNI.EQ.1) THEN
          WRITE(6,*) ' Element corresponding to unit-operator',
     &    T(N_ZERO_EI+1)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE PRINT_T_EI_SLAVE(T,ITSYM,N_EXTP,
     &           NDIM_EX_ST,NDIM_IN_ST,NSMOB)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Specific input
      DIMENSION T(*)
      DIMENSION NDIM_EX_ST(NSMOB,N_EXTP),NDIM_IN_ST(NSMOB,N_EXTP)
*
      IOFF = 1
      DO I_EXTP = 1, N_EXTP
       DO I_EXSM = 1, NSMOB
        I_INSM = MULTD2H(I_EXSM,ITSYM)
        WRITE(6,*) ' Block with external type and sym :', I_EXTP, I_EXSM
        N_EX = NDIM_EX_ST(I_EXSM,I_EXTP)
        N_IN = NDIM_IN_ST(I_INSM,I_EXTP)
C?      WRITE(6,*) ' N_EX, N_IN = ', N_EX, N_IN
        CALL WRTMAT(T(IOFF),N_IN,N_EX,N_IN,N_EX)
        IOFF = IOFF + N_EX*N_IN
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GET_DIAG_H0_EI(DIAG)
*
* Construct diagonal of H0 over orthonormal states
*
*. Jeppe Olsen, March 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'lucinp.inc'
*. Output
      DIMENSION DIAG(*)
*
      NTEST = 100
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_DIA')
C?    WRITE(6,*) ' Entering GET_DIAG.., WORK(KINT1) = ', WORK(KINT1)
*. Generate H0 as FI + FA - starting from scatch to be on the safe side..
      CALL COPVEC(WORK(KH),WORK(KHINA),NINT1)
*. Inactive Fock matric
      CALL FISM(WORK(KHINA),ECC)
*. Active Fock matrix
      CALL FAM(WORK(KFIFA))
*. and add 
      ONE = 1.0D0
      CALL VECSUM(WORK(KFIFA),WORK(KFIFA),WORK(KHINA),
     &            ONE,ONE,NINT1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' FIFA from GET_DIAG_H0_EI: '
        CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*. External part of zero-order energy: External part is assumed to 
*. be doubly occupied hole space and unoccupied particle space
      E0REF_EXT = 0.0D0
*. Obtain diagonal of H0
      CALL MEMMAN(KLH0DIAS,NTOOB,'ADDL  ',2,'H0DIAS')
      CALL MEMMAN(KLH0DIA,NTOOB,'ADDL  ',2,'H0DIA ')
      CALL GET_DIAG_BLOC_MAT(WORK(KFIFA),WORK(KLH0DIAS),NSMOB,NTOOBS,1)
C          GET_DIAG_BLOC_MAT(A,ADIAG,NBLOCK,LBLOCK,ISYM)
      IF(NTEST.GE.100) THEN
      WRITE(6,*)' Diagonal of FIFA in sym-order'
      CALL WRTMAT(WORK(KLH0DIAS),1,NTOOB,1,NTOOB)
      END IF
*. Was obtained in symmetry ordered basis, type used below so
      DO I= 1,NTOOB    
       WORK(KLH0DIA-1+IREOST(I)) = WORK(KLH0DIAS-1+I)
      END DO
        
*. 4 Blocks for occupation of  external strings
      CALL MEMMAN(KL_IST_EX_CA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STE_CA')
      CALL MEMMAN(KL_IST_EX_CB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STE_CB')
      CALL MEMMAN(KL_IST_EX_AA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STE_AA')
      CALL MEMMAN(KL_IST_EX_AB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STE_AB')
*
C     GET_DIAG_H0_EI_SLAVE(DIAG,
C    &           ISPOBEX_TP,N_EXTP,E0REF_EXT,
C    &           N_ORTN_FOR_SE,
C    &           I_IN_TP,IB_INTP,IB_EXTP,H0DIAG,
C    &           IST_EX_CA, IST_EX_CB, IST_EX_AA,IST_EX_AB)
C?    WRITE(6,*) ' I_INT_OFF before callto GET_DI..SLAVE',I_INT_OFF
      CALL GET_DIAG_H0_EI_SLAVE(DIAG,WORK(KLSOBEX),N_EXTOP_TP,E0REF_EXT,
     &     WORK(KL_N_ORTN_FOR_SE),I_IN_TP,
     &     I_INT_OFF,I_EXT_OFF,WORK(KLH0DIA),
     &     WORK(KL_IST_EX_CA),WORK(KL_IST_EX_CB),
     &     WORK(KL_IST_EX_AA),WORK(KL_IST_EX_AB))
*. Reinstall one-electron integrals
      CALL COPVEC(WORK(KINT1O),WORK(KFIFA),NINT1)
C?    WRITE(6,*) ' Leaving GET_DIAG.., WORK(KINT1) = ', WORK(KINT1)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT_DIA')
*
      RETURN
      END
      SUBROUTINE GET_DIAG_H0_EI_SLAVE(DIAG,
     &           ISPOBEX_TP,N_EXTP,E0REF_EXT,
     &           N_ORTN_FOR_SE,
     &           I_IN_TP,IB_INTP,IB_EXTP,H0DIAG,
     &           IST_EX_CA, IST_EX_CB, IST_EX_AA,IST_EX_AB)
*
* Generate diagonal of H0 over zero-order states using EI approach 1
*
* Jeppe Olsen, March 2009, trying to do a bit of science
*
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc' 
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
*. E0REF_EXT is external parts of E0 of reference state
*. Input: 
      INTEGER ISPOBEX_TP(NGAS,4,*)
*. Number of orthonormal internal states per symmetry and external type
      INTEGER N_ORTN_FOR_SE(NSMOB,N_EXTP)
*. Diagonal part of H0 in orbital basis
      DIMENSION H0DIAG(*)
*. Output
      DIMENSION DIAG(*)
*. Scratch
*. Blocks for holding group of external strings with given sym
      DIMENSION IST_EX_CA(*),IST_EX_CB(*), IST_EX_AA(*),IST_EX_AB(*)
*. Local scratch
      INTEGER IGRP(MXPNGAS)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_DIA')
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' -------------- '
        WRITE(6,*) ' GET_DIAG_H0_EI '
        WRITE(6,*) ' -------------- '
        WRITE(6,*)
C?      WRITE(6,*) ' E0REF_EXT at start of GET_DIAG.. ', E0REF_EXT
      END IF
*. Symmetry of O+(e)O+(i)
      ISYM_EI  = 1
*. Largest number of internal zero-order states for given symmetry
C          IMNXVC(IVEC,NDIM,MXMN,IVAL,IPLACE)
      CALL IMNXVC(N_ORTN_FOR_SE,NSMOB*N_EXTP,1,NMAX_ORTN_FOR_SE,IPLACE)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NMAX_ORTN_FOR_SE = ', NMAX_ORTN_FOR_SE
        CALL IWRTMA(N_ORTN_FOR_SE,NSMOB,N_EXTP,NSMOB,N_EXTP)
      END IF
*. an array for internal state energies
      CALL MEMMAN(KLE0INT,NMAX_ORTN_FOR_SE,'ADDL  ',2,'E0_INT')
*
      IEI_ORTN = 0
      DO J_EXTP = 1, N_EXTP
        J_EXTP_ABS = IB_EXTP-1+J_EXTP
*
        NEL_EX_CA = IELSUM(ISPOBEX_TP(1,1,J_EXTP_ABS),NGAS)
        NEL_EX_CB = IELSUM(ISPOBEX_TP(1,2,J_EXTP_ABS),NGAS)
        NEL_EX_AA = IELSUM(ISPOBEX_TP(1,3,J_EXTP_ABS),NGAS)
        NEL_EX_AB = IELSUM(ISPOBEX_TP(1,4,J_EXTP_ABS),NGAS)
*
*. Loop over strings in EI order
        DO I_EX_SM = 1, NSMOB
          I_IN_SM = MULTD2H(I_EX_SM,ISYM_EI)
          IF(NTEST.GE.20) WRITE(6,*) ' J_EXTP, I_EX_SM,I_IN_SM =',
     &                                 J_EXTP, I_EX_SM,I_IN_SM
*. Obtain diagonal of internal energy contributions
C         GET_BLOCK_OF_HS_IN_INTERNAL_SPACE(
C    &    IEXTP,IINTSM,I_HS,HSBLCK,I_INT_TP,I_ONLY_DIA)
          CALL GET_BLOCK_OF_HS_IN_INTERNAL_SPACE(
     &    J_EXTP,I_IN_SM,1,WORK(KLE0INT),I_INT_TP,1)
          IF(NTEST.GE.20) THEN
            L_INTOP = N_ORTN_FOR_SE(I_IN_SM,J_EXTP)
            WRITE(6,*) ' H0 energies of internal states'
            CALL WRTMAT(WORK(KLE0INT),1,L_INTOP,1,L_INTOP)
          END IF
*. Loop over symmetries of external operators
          DO I_EX_C_SM = 1, NSMOB
           I_EX_A_SM = MULTD2H(I_EX_C_SM,I_EX_SM)
           DO I_EX_CA_SM = 1, NSMOB
            I_EX_CB_SM = MULTD2H(I_EX_CA_SM,I_EX_C_SM)
             DO I_EX_AA_SM = 1, NSMOB
              I_EX_AB_SM = MULTD2H(I_EX_AA_SM,I_EX_A_SM)
*
              IF(NTEST.GE.1000) THEN
                WRITE(6,'(A,4I4)') 
     &          'I_EX_CA_SM, I_EX_CB_SM, I_EX_AA_SM, I_EX_AB_SM = ',
     &           I_EX_CA_SM, I_EX_CB_SM, I_EX_AA_SM, I_EX_AB_SM
              END IF
*. Occupation of external strings
*. CA in IST_EX_CA
              CALL OCC_TO_GRP(ISPOBEX_TP(1,1,J_EXTP_ABS),IGRP,1)
              CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &             I_EX_CA_SM,NEL_EX_CA,NSTR_EX_CA,IST_EX_CA,NTOOB,0,
     &             IDUM,IDUM)
*. CB in IST_EX_CB
              CALL OCC_TO_GRP(ISPOBEX_TP(1,2,J_EXTP_ABS),IGRP,1)
              CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &             I_EX_CB_SM,NEL_EX_CB,NSTR_EX_CB,IST_EX_CB,NTOOB,0,
     &             IDUM,IDUM)
*. AA in IST_EX_AA
              CALL OCC_TO_GRP(ISPOBEX_TP(1,3,J_EXTP_ABS),IGRP,1)
              CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &             I_EX_AA_SM,NEL_EX_AA,NSTR_EX_AA,IST_EX_AA,NTOOB,0,
     &             IDUM,IDUM)
*. AB in IST_EX_AB
              CALL OCC_TO_GRP(ISPOBEX_TP(1,4,J_EXTP_ABS),IGRP,1)
              CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &             I_EX_AB_SM,NEL_EX_AB,NSTR_EX_AB,IST_EX_AB,NTOOB,0,
     &             IDUM,IDUM)
C?      WRITE(6,*) ' End of external '
*. Loop over External strings (CA, CB, AA, AB)
              DO I_EX_AB = 1, NSTR_EX_AB
*. Contribution from I_EX_AB to diagonal
C             E1_FOR_STRING(HDIAG,ISTRING,NEL)
              D_EX_AB =
     &        E1_FOR_STRING(H0DIAG,IST_EX_AB(1+(NEL_EX_AB)*(I_EX_AB-1)),
     &                      NEL_EX_AB)
              DO I_EX_AA = 1, NSTR_EX_AA
              D_EX_AA =
     &        E1_FOR_STRING(H0DIAG,IST_EX_AA(1+(NEL_EX_AA)*(I_EX_AA-1)),
     &                      NEL_EX_AA)
              DO I_EX_CB = 1, NSTR_EX_CB
              D_EX_CB =
     &        E1_FOR_STRING(H0DIAG,IST_EX_CB(1+(NEL_EX_CB)*(I_EX_CB-1)),
     &                      NEL_EX_CB)
              DO I_EX_CA = 1, NSTR_EX_CA
              D_EX_CA =
     &        E1_FOR_STRING(H0DIAG,IST_EX_CA(1+(NEL_EX_CA)*(I_EX_CA-1)),
     &                      NEL_EX_CA)
                E_EXT = E0REF_EXT + D_EX_CB +D_EX_CA -D_EX_AB - D_EX_AA 
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' D_EX_CB,D_EX_CA,D_EX_AB,D_EX_AA =',
     &                         D_EX_CB,D_EX_CA,D_EX_AB,D_EX_AA
                END IF
CM*. Obtain diagonal of internal energy contributions
C?              CALL GET_BLOCK_OF_HS_IN_INTERNAL_SPACE(
C?   &          J_EXTP,I_IN_SM,1,WORK(KLE0INT),I_INT_TP,1)
*
                L_INTOP = N_ORTN_FOR_SE(I_IN_SM,J_EXTP)
                DO I_INTOP = 1, L_INTOP
                 IEI_ORTN = IEI_ORTN + 1
                 IF(NTEST.GE.1000) THEN
                   WRITE(6,*) ' I, Internal, external term =',
     &             WORK(KLE0INT-1+I_INTOP), E_EXT
                 END IF
                 DIAG(IEI_ORTN) = WORK(KLE0INT-1+I_INTOP) + E_EXT
                 
                END DO
*.              ^ End of loop over internal states
              END DO
              END DO
              END DO
              END DO
*             ^ End of loop over external CA,CB,AA,AB strings of given sym
            END DO
           END DO
          END DO
*.        ^ End of loop over symmetry of external CA, CB, AA, AB strings
        END DO
*.      ^ End of loop over symmetry of external operators
      END DO
*.    ^ End of loop over external types
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Diagonal of zero-order energy of ortn states'
        CALL WRTMAT(DIAG,1,IEI_ORTN,1,IEI_ORTN)
      END  IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',1,'GT_DIA')
*
      RETURN
      END 
      SUBROUTINE LUCIA_IC_EI
     &           (IREFSPCE,ITREFSPC,ICTYP,EREF,I_DO_CUMULANTS,
     &            EFINAL,CONVER,VNFINAL)
*
*
* Master routine for internally contracted calculation 
*
* Specialized for reference states that allows division into 
* internal and external parts
*
* Jeppe Olsen, September 06
*
* The spaces are assumed organized as
* 1 : Reference space on which excitations are performed (IREFSPC)
* 2 : Space that defines excitations                    (ITREFSPC)
* 3 : Space where calculation is performed               (ITREFSPC)
*
* Deviations from this will cause trouble
* 
      INCLUDE 'wrkspc.inc'
      REAL*8 
     &INPROD
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'newccp.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. A bit of local scratch
      DIMENSION ICASCR(MXPNGAS)
      CHARACTER*6 ICTYP
*
      EXTERNAL MTV_FUSK, STV_FUSK
      EXTERNAL H_S_EFF_ICCI_TV,H_S_EXT_ICCI_TV
      EXTERNAL HOME_SD_INV_T_ICCI
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICCI  ')
*
      I_SPIN_ADAPT = 0
      IREFSPC = IREFSPCE
      MSCOMB_CC = 0
*
C?    WRITE(6,*) ' I_INC_AA = ', I_INC_AA
      IF(I_INC_AA.EQ.0) THEN
        NOAAEX =  1
      ELSE
        NOAAEX = 0
      END IF
*. Form of internal states
* 1 => diagonalize metric
* 2 => Diagonalize zero-order Hamiltonian matrix
* 3 => Diagonalize zero-order Jacobian => gives different left and right
*      zero-order stated
      I_IN_TP = 2
*. This is using splitting into internals and externals so
      I_DO_EI = 1
* Internal inactive -> inactive, sec -> sec will be eliminated
      IF(ICEXC_INT.NE.0) THEN
        WRITE(6,*) 
     &  ' inactive -> inactive, sec -> sec excitations turned off'
        ICEXC_INT = 0
      END IF
*
      NTEST = 2
      IPRNT = NTEST
      IF(NTEST.GE.1) THEN
         WRITE(6,*)
         WRITE(6,*) ' Internal contracted section entered '
         WRITE(6,*) ' ==================================== '
         WRITE(6,*)
         WRITE(6,*) ' Version exploiting external/internal division'
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,'(A,A)') ' Form of calculation  ', ICTYP
         WRITE(6,*)
         WRITE(6,*) '  Symmetri of reference vector ' , IREFSM 
         WRITE(6,*) '  Space of Reference vector ', IREFSPC
         WRITE(6,*) '  Space of Internal contracted vector ', ITREFSPC-1
         WRITE(6,*)
         WRITE(6,*) ' Parameters defining operator manifold:'
         WRITE(6,*) '       Max operator rank   ', ICOP_RANK_MAX
         WRITE(6,*) '       Max number of active indeces ', 
     &   ICEXC_MAX_ACT
         WRITE(6,*) '       Max number of external indeces ', 
     &   ICEXC_MAX_EXT
         IF(ICEXC_INT.EQ.1) THEN
           WRITE(6,*) '       ',
     &   'Internal (ina->ina, sec->sec) excitations allowed'
         ELSE
           WRITE(6,*) '       ',
     &   'Internal (ina->ina, sec->sec) excitations not allowed'
         END IF
         IF(NOAAEX.EQ.1) THEN
           WRITE(6,*) '       ',
     &     'Pure active excitations are not included'
         ELSE
           WRITE(6,*) '       ',
     &     ' Pure active excitations are included'
         END IF
         WRITE(6,*)
         WRITE(6,*) 
     &   ' Threshold for nonsingular  eigenvalues of metric',
     &     THRES_SINGU
*
         IF(IRESTRT_IC.EQ.1) THEN
           WRITE(6,*) ' Restarted calculation : '
           WRITE(6,*) '      IC coefficients read from LUSC54'
           WRITE(6,*) '      CI for reference read from LUEXC '
         END IF
         WRITE(6,*) ' Zero-order states obtained by:'
         IF(I_IN_TP.EQ.1) THEN
           WRITE(6,*) ' diagonalizing metric '
         ELSE IF(I_IN_TP.EQ.2) THEN
           WRITE(6,*) ' Diagonalizing zero-order Hamiltonian matrix'
         ELSE IF (I_IN_TP.EQ.3) THEN
           WRITE(6,*) ' Diagonalizing zero-order Jacobian matrix'
         END IF
      END IF
*
*. Divide orbital spaces into inactive, active, secondary using 
*. space 1
      CALL CC_AC_SPACES(1,IREFTYP)
*. Obtain the orbital excitations
* (copied more or less from LUCIA_GENCC)
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*
COLD  ICSPC = ICISPC
COLD  ISSPC = ICISPC
*
C?    WRITE(6,*) ' Zero-order Hamiltonian with zero-order density '
C. IPHGAS1 should be used to divide into H,P,V, but IPHGAS is used, so swap
      CALL ISWPVE(IPHGAS(1),IPHGAS1(1),NGAS)
      CALL COPVEC(WORK(KINT1O),WORK(KFIFA),NINT1)
      CALL FIFAM(WORK(KFIFA))
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' FIFA: '
        CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*. External part of zero-order energy: External part is assumed to 
*. be doubly occupied hole space and unoccupied particle space
C     GET_E0REF_EXT(FI,IPHGASX)
C     E0_REF_EXT = GET_E0REF_EXT(WORK(KFIFA),IPHGAS)
*. And clean up
      CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)

*
*
* ========================
* info for reference space
* ========================
*
*. Make sure that there is just a single occupation space
      CALL OCCLSE(1,NOCCLS_REF,IOCCLS,NEL,IREFSPC,0,0,NOBPT)
      IF(NOCCLS_REF.NE.1) THEN
        WRITE(6,*) ' Problem in LUCIA_IC_NEW : '
        WRITE(6,*)
     &  ' Reference space is not a single occupation space'
        STOP
     &  ' Reference space is not a single occupation space'
      END IF
*. and the reference occupation space
      CALL MEMMAN(KLOCCLS_REF,NGAS,'ADDL  ',1,'OCC_RF')
      CALL OCCLSE(2,NOCCLS_REF,WORK(KLOCCLS_REF),NEL,IREFSPC,0,0,NOBPT)
*
* ====================================
* Info for space defining excitations
* ====================================
*
*. Number
C     IT2REFSPC = ITREFSPC - 1
      IT2REFSPC = ITREFSPC 
      CALL OCCLSE(1,NOCCLS,IOCCLS,NEL,IT2REFSPC,0,0,NOBPT)
*. And the occupation classes
      CALL MEMMAN(KLOCCLS,NOCCLS*NGAS,'ADDL  ',1,'OCCLS ')
      CALL OCCLSE(2,NOCCLS,WORK(KLOCCLS),NEL,IT2REFSPC,0,0,NOBPT)
*. Number of occupation classes for T-operators
      NTOCCLS = NOCCLS
*. It could be an idea to check that reference space is included
* ========================
* Orbital excitation types
* ========================
*
*. Number of excitation types
      IFLAG = 1
      IDUM = 1
*
      MX_NCREA = ICOP_RANK_MAX
      MX_NANNI = ICOP_RANK_MAX
      MX_AAEXC = ICEXC_MAX_ACT
      I_OOCC = 0
*. Pure AA excitations (without external part?)
      IF(I_INC_AA.EQ.0) THEN
        NOAAEX =  1
      ELSE
        NOAAEX = 0
      END IF
C?    WRITE(6,*) ' NOAAEX =', NOAAEX
      IFLAG = 1
      CALL TP_OBEX3(NOCCLS,NEL,NGAS,WORK(IDUM),
     &             WORK(IDUM),WORK(IDUM),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(IDUM),MX_AAEXC,IFLAG,
     &             I_OOCC,NOBEX_TP,NOAAEX,ICEXC_MAX_EXT,IPRCC)
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' NOBEX_TP,MX_EXC_LEVEL = ', NOBEX_TP,MX_EXC_LEVEL
*. And the actual orbital excitations
*.  An orbital excition operator is defined by
*   1 : Number of creation operators
*   2 : Number of annihilation operators
*   3 : The actual creation and annihilation operators
*. The number of orbital excitations is increased by one to include
*. excitations within the reference space
      NOBEX_TPE = NOBEX_TP+1
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TPE,'ADDL  ',1,'LCOBEX')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TPE,'ADDL  ',1,'LAOBEX')
      CALL MEMMAN(KOBEX_TP ,NOBEX_TPE*2*NGAS,'ADDL  ',1,'IOBE_X')
*. Excitation type => Original occupation class
      CALL MEMMAN(KEX_TO_OC,NOBEX_TPE,'ADDL  ',1,'EX__OC')
      IFLAG = 0
*. (Unit operator is created even if only NOBEX_TP is transferred
      CALL TP_OBEX3(NOCCLS,NEL,NGAS,WORK(KOBEX_TP),
     &     WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),
     &     WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &     MX_EXC_LEVEL,WORK(KEX_TO_OC),MX_AAEXC,IFLAG,
     &     I_OOCC,NOBEX_TP,NOAAEX,ICEXC_MAX_EXT,IPRCC)

*
* =======================
* Spinorbital excitations  
* =======================
*
*. Spin combinations of CC excitations : Currently we assume that
*. The T-operator is a singlet, can 'easily' be changed
*
*. Notice : The first time in OBEX_TO_SPOBEX we always use MSCOMB_CC = 0.
*. This may lead to the allocation of too much space for
*. spinorbital excitations, but MSCOMB_CC = 1, requires access
*. to WORK(KLSOBEX) which has not been defined
*
*
* Combined external-internal excitations
*. Largest spin-orbital excitation level
      IF(MXSPOX.NE.0) THEN
        MXSPOX_L = MXSPOX
      ELSE
        MXSPOX_L = MX_EXC_LEVEL
      END IF
      IF(NTEST.GE.10) 
     &WRITE(6,*) ' MXSPOX, MXSPOX_L, MX_EXC_LEVEL = ',
     &             MXSPOX, MXSPOX_L, MX_EXC_LEVEL
      IZERO = 0
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP,NGAS,
     &     NOBPT,0,IZERO ,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     MXSPOX_L,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC)
*. Extended number of spin-orbital excitations : Include
*. unit operator as last spinorbital excitation operator
      NSPOBEX_TPE = NSPOBEX_TP + 1
      IF(IPRNT.GE.1) THEN
        WRITE(6,*) ' Number of spinorbital excitations(with unit)',
     &  NSPOBEX_TPE
      END IF
      IF(IPRNT.GE.10) WRITE(6,*) ' NSPOBEX_TP = ', NSPOBEX_TP
*. Allocate space for EI, E, I (external-internal, external, internal ) 
*. spinorbital excitations. As the number of E and I are not known 
*. we set their dimension to NSPOBEX_TP ( an upper limit)
*. And the actual spinorbital excitation operators
      CALL MEMMAN(KLSOBEX,4*NGAS*3*NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,3*NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TP+1,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TP+1,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TP+1,'ADDL  ',1,'IBSOXF')
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     NOBPT,0,0,IAAEXC_TYP,IACT_SPC,
     &     IPRCC,WORK(KLSOX_TO_OX),MXSPOX_L,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC)
      NSPOBEX_TPE = NSPOBEX_TP + 1
*. Add unit-operator as last spinorbital excitation
      IZERO = 0
      CALL ISTVC3(WORK(KLSOBEX),(NSPOBEX_TPE-1)*4*NGAS+1,IZERO,4*NGAS)
      IF(IPRNT.GE.5) THEN
        WRITE(6,*) ' Extended list of spin-orbital excitations : '
        CALL WRT_SPOX_TP(WORK(KLSOBEX),NSPOBEX_TPE)
      END IF
*
*. Construct the various external and internal operators of the 
*. various operatortpys and set up mappings from IE operator 
*. to I,E operators
*. Question Jeppe : Should the zero-operator also be splitted. 
*. Yes, I am doing this from today -aug20
      I_EXT_OFF =  NSPOBEX_TPE + 1
*. Offset in KLSOBEX to internal part is obtained in SPLIT_IE_SPOBEXTP
      I_INT_OFF =  0
*. The various internal operators for the same external operators 
*. will be collected. Mappings for this
*. Number of internal types per external type
      CALL MEMMAN(KL_N_INT_FOR_EXT,NSPOBEX_TPE,'ADDL  ',1,'N_INT_')
*. Offsets for internals for given external
      CALL MEMMAN(KL_IB_INT_FOR_EXT,NSPOBEX_TPE,'ADDL  ',1,'IB_INT')
*. And the actual internals for each external 
      CALL MEMMAN(KL_I_INT_FOR_EXT,NSPOBEX_TPE,'ADDL  ',1,'I_INT_')
*   
      CALL SPLIT_IE_SPOBEXTP(WORK(KLSOBEX),NSPOBEX_TPE, N_EXTOP_TP,
     &     N_INTOP_TP,I_EXT_OFF, I_INT_OFF,WORK(KL_N_INT_FOR_EXT),
     &     WORK(KL_IB_INT_FOR_EXT), WORK(KL_I_INT_FOR_EXT),NGAS,
     &     IHPVGAS )
*. Obtain reorder array EI-order => standard order
       CALL MEMMAN(KL_I_EI_TP_REO,NSPOBEX_TPE,'ADDL  ',1,'EITPRE')
C?     WRITE(6,*) ' NSPOBEX_TPE before EITP.. =', NSPOBEX_TPE
       CALL EITP_TO_SPOXTP_REO(WORK(KL_I_EI_TP_REO),
     &      WORK(KLSOBEX),NSPOBEX_TPE,
     &      I_EXT_OFF,I_INT_OFF,NGAS,N_EXTOP_TP,N_INTOP_TP,
     &      WORK(KL_N_INT_FOR_EXT), WORK(KL_IB_INT_FOR_EXT),
     &      WORK(KL_I_INT_FOR_EXT) )
C?     WRITE(6,*) ' I_INT_OFF after sec call to EITP_TO.. ', I_INT_OFF
C     EITP_TO_SPOXTP_REO(I_EI_TP_REO,ISPOBEX_TP,NSPOBEX_TP,
C    &           IB_EXTP, IB_INTP,NGAS,N_EXTP,N_INTP,
C    &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX)
*. We have now stored information about the external types 
*. in WORK(KLSOBEX) from I_EXT_OFF and for external types in 
*. I_INT_OFF. We do however not increase NSPOBEX_TP,
*. so in the following there are hopefully invisible
*. 
*. Obtain info on the dimension of the various internal and external types
*
C      DIMENSION_EI_EXP(ISPOBEX_TP,IB_EXTP,IB_INTP,N_EXTP,
C    &            N_INTP,N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,
C    &            NDIM_EX_ST,NDIM_IN_ST,NDIM_EI,NDIM_IN_SE, NSMOB)
      CALL MEMMAN(KL_NDIM_EX_ST,NSMOB*N_EXTOP_TP,'ADDL  ',1,'NDIM_E')
      CALL MEMMAN(KL_NDIM_IN_ST,NSMOB*N_INTOP_TP,'ADDL  ',1,'NDIM_I')
      CALL MEMMAN(KL_NDIM_IN_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,'NDIMSE')
      CALL DIMENSION_EI_EXP(WORK(KLSOBEX),I_EXT_OFF,I_INT_OFF,
     &     N_EXTOP_TP,N_INTOP_TP,WORK(KL_N_INT_FOR_EXT),
     &     WORK(KL_IB_INT_FOR_EXT), WORK(KL_I_INT_FOR_EXT),
     &     WORK(KL_NDIM_EX_ST),WORK(KL_NDIM_IN_ST),NDIM_EI,
     &     WORK(KL_NDIM_IN_SE),NSMOB)
*
      CALL ISTVC3(WORK(KLSOX_TO_OX),NSPOBEX_TPE,NOBEX_TP+1,1)
*. Mapping spinorbital excitations => occupation classes
      CALL MEMMAN(KIBSOX_FOR_OCCLS,NOCCLS,'ADDL  ',1,'IBSXOC')
      CALL MEMMAN(KNSOX_FOR_OCCLS,NOCCLS,'ADDL  ',1,' NSXOC')
      CALL MEMMAN(KISOX_FOR_OCCLS,NSPOBEX_TPE,'ADDL  ',1,' ISXOC')
C       SPOBEX_FOR_OCCLS(
C    &           IEXTP_TO_OCCLS,NOCCLS,ISOX_TO_OX,NSOX,
C    &           NSOX_FOR_OCCLS,ISOX_FOR_OCCLS,IBSOX_FOR_OCCLS)
      CALL SPOBEX_FOR_OCCLS(WORK(KEX_TO_OC),NOCCLS,WORK(KLSOX_TO_OX),
     &     NSPOBEX_TPE,WORK(KNSOX_FOR_OCCLS),WORK(KISOX_FOR_OCCLS),
     &     WORK(KIBSOX_FOR_OCCLS))
*
*. Frozen spin-orbital excitation types
      CALL MEMMAN(KLSPOBEX_FRZ, NSPOBEX_TPE,'ADDL  ',1,'SPOBFR')
      CALL FRZ_SPOBEX(WORK(KLSPOBEX_FRZ),WORK(KLCOBEX_TP),NSPOBEX_TP,
     &                WORK(KLSOX_TO_OX),IFRZ_CC_AR,NFRZ_CC)
      IZERO = 0
      CALL ISTVC3(WORK(KLSPOBEX_FRZ),NSPOBEX_TPE,IZERO,1)
*. Spin-orbital excitation types related by spin-flip
      CALL MEMMAN(KLSPOBEX_AB,NSPOBEX_TPE,'ADDL  ',1,'SPOBAB')
      CALL SPOBEXTP_PAIRS(NSPOBEX_TPE,WORK(KLSOBEX),NGAS,
     &                    WORK(KLSPOBEX_AB))
C          SPOBEXTP_PAIRS(NSPOBEX_TP,ISPOBEX,NGAS,ISPOBEX_PAIRS)
C     SELECT_AB_TYPES(NSPOBEX_TP,ISPOBEX_TP,
C    &           ISPOBEX_PAIRS,NGAS)
      CALL SELECT_AB_TYPES(NSPOBEX_TPE,WORK(KLSOBEX),
     &                     WORK(KLSPOBEX_AB),NGAS)
*. Alpha- and beta-excitations constituting the spinorbital excitations
*. Number
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMMAN(KLAOBEX,LENA,'ADDL  ',2,'IAOBEX')
      CALL MEMMAN(KLBOBEX,LENB,'ADDL  ',2,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     0,NAOBEX_TP,NBOBEX_TP,WORK(KLAOBEX),WORK(KLBOBEX))
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(IPRNT.GE.10) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
* Max Dimension of spinorbital excitation operators
*
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TPE,'ADDL  ',1,'SPOBAC')
*. ALl spinorbital excitations are initially active
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TPE)
*
      MX_ST_TSOSO_MX = 0
      MX_ST_TSOSO_BLK_MX = 0
      MX_TBLK_MX = 0
      MX_TBLK_AS_MX = 0
      LEN_T_VEC_MX = 0
      DO ICCAMP_SM = 1, NSMST
*. Dimension without zero-particle operator
        CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ICCAMP_SM,
     &       MX_ST_TSOSOL,MX_ST_TSOSO_BLKL,MX_TBLKL,
     &       WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VECL,
     &       MSCOMB_CC,MX_TBLK_AS,
     &       WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &       NTCONF,IPRCC)
*
        MX_ST_TSOSO_MX = MAX(MX_ST_TSOSO_MX,MX_ST_TSOSOL)
        MX_ST_TSOSO_BLK_MX = MAX(MX_ST_TSOSO_BLK_MX,MX_ST_TSOSO_BLKL)
        MX_TBLK_MX = MAX(MX_TBLK_MX,MX_TBLKL)
        MX_TBLK_AS_MX = MAX(MX_TBLK_AS_MX,MX_TBLK_AS)
        LEN_T_VEC_MX = MAX(LEN_T_VEC_MX, LEN_T_VECL)
*
      END DO
      IF(IPRNT.GE.10) WRITE(6,*) ' MX_TBLK_AS_MX = ', MX_TBLK_AS_MX
      IF(IPRNT.GE.10) WRITE(6,*) ' LEN_T_VEC_MX = ', LEN_T_VEC_MX
*. And dimensions for symmetry 1
      ITOP_SM = 1
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_SBSTR,
     &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRCC)
*. NCAAB and N_CC_AMP does not include zero-particle operator
      N_CC_AMP = LEN_T_VEC
      NCAAB = N_CC_AMP
      IF(IPRNT.GE.5)
     &WRITE(6,*) ' LUCIA_GENCC : N_CC_AMP = ', N_CC_AMP
      IF(IPRNT.GE.5) THEN
        WRITE(6,*) ' Number of amplitudes per operator type: '
        CALL IWRTMA(WORK(KLLSOBEX),NSPOBEX_TP,1,NSPOBEX_TP,1)
      END IF
*. Hard wire info for unit operator stored as last spinorbital excitation
C  ISTVC2(IVEC,IBASE,IFACT,NDIM)
      IONE = 1
      CALL ISTVC3(WORK(KLLSOBEX),NSPOBEX_TPE,IONE,1)
      N_CC_AMPP1 = N_CC_AMP + 1
      CALL ISTVC3(WORK(KLIBSOBEX),NSPOBEX_TPE,N_CC_AMPP1,1)
*. Obtain mapping between EI order and standard order 
      CALL MEMMAN(KL_IREO_EI_ST,N_CC_AMP+1,'ADDL  ',1,'IREOST')
C?    WRITE(6,*) ' I_EXT_OFF, I_INT_OFF  = ', I_EXT_OFF,I_INT_OFF
      CALL EI_REORDER_ARRAYS(WORK(KL_IREO_EI_ST),WORK(KLSOBEX),
     &     I_EXT_OFF,I_INT_OFF,WORK(KLLSOBEX),WORK(KLIBSOBEX),
     &     NSPOBEX_TPE,N_EXTOP_TP,N_INTOP_TP,
     &     WORK(KL_N_INT_FOR_EXT),WORK(KL_IB_INT_FOR_EXT),
     &     WORK(KL_I_INT_FOR_EXT),1)
C     EI_REORDER_ARRAYS(IREO_EI_ST,
C    &           ISPOBEX_TP,IB_EXTP,IB_INTP,
C    &           L_SPOBEX_TP,IB_SPOBEX_TP,NSPOBEX_TP,
C    &           N_EXTP,N_INTP,
C    &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,ISYM)
*. Determine dimensions of the various internal and external 
*. types
* =============
* Scratch space
* =============
*
       CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
       KVEC1P = KVEC1
       KVEC2P = KVEC2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Standard list of CAAB operators'
        CALL PRINT_CAAB_LIST(1)
      END IF
*
*
* --------------------------
* Obtain the internal states
* --------------------------
*
      CALL MEMMAN(KL_N_ORTN_FOR_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,'ORT_SE')
      CALL MEMMAN(KL_N_INT_FOR_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,'INT_SE')
*. Dimension of matrices with given symmetry and external types
C?    WRITE(6,*) ' NDIM_IN_SE array : '
C?    CALL IWRTMA(WORK(KL_NDIM_IN_SE),
C?   &     NSMOB,N_EXTOP_TP,NSMOB,N_EXTOP_TP)
*
      I_SE_DIM2 = ISQELSUM(WORK(KL_NDIM_IN_SE),NSMOB*N_EXTOP_TP,0)
      WRITE(6,*) ' Dimension of transformation matrices X1, X2',
     &            I_SE_DIM2
      I_SE_DIM2S= ISQELSUM(WORK(KL_NDIM_IN_SE),NSMOB*N_EXTOP_TP,1)
      WRITE(6,*) ' Dimension of internal overlap matrix ',
     &            I_SE_DIM2S
      I_SE_DIM = IELSUM(WORK(KL_NDIM_IN_SE),NSMOB*N_EXTOP_TP)
*. Space for transformation matrices: For each block matrices 
*. X1, X2 and a vector sigma
      CALL MEMMAN(KL_X1_INT_EI_FOR_SE,I_SE_DIM2,'ADDL  ',2,'X1INEI')
      CALL MEMMAN(KL_X2_INT_EI_FOR_SE,I_SE_DIM2,'ADDL  ',2,'X2INEI')
      CALL MEMMAN(KL_SG_INT_EI_FOR_SE,I_SE_DIM ,'ADDL  ',2,'SGINEI')
      CALL MEMMAN(KL_S_INT_EI_FOR_SE,I_SE_DIM2S,'ADDL  ',2,'SINEI ')
      CALL MEMMAN(KL_IBX1_INT_EI_FOR_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,
     &            'IBX1IN')
      CALL MEMMAN(KL_IBX2_INT_EI_FOR_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,
     &            'IBX2IN')
      CALL MEMMAN(KL_IBSG_INT_EI_FOR_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,
     &            'IBSGIN')
      CALL MEMMAN(KL_IBS_INT_EI_FOR_SE,NSMOB*N_EXTOP_TP,'ADDL  ',1,
     &            'IBSIN ')
*. For left hand side transformation to internal state basis
*. if Jacobian is diagonalized
*. matrix diagonalizing metric is common for L and R as is sigma, 
*. so only an extra matrix for the last diagonalization is needed
      IF(I_IN_TP.LE.2) THEN
        KL_X2L_INT_EI_FOR_SE = KL_X2_INT_EI_FOR_SE
      ELSE
        CALL MEMMAN(KL_X2L_INT_EI_FOR_SE,I_SE_DIM2,'ADDL  ',2,'XRINEI')
      END IF
*
*.---------------------------------------------
*. Obtain internal states by diagonalizing H0
*.---------------------------------------------
C?    WRITE(6,*) ' Call to IDIM_TCC just before zero-order states'
C?    CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,1,
C?   &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
C?   &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
C?   &     MSCOMB_CC,MX_SBSTR,
C?   &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
C?   &     NTCONF,IPRCC)
*
C?    WRITE(6,*) ' Before GET_INT, I_INT_OFF,I_EXT_OFF ',
C?   &                             I_INT_OFF,I_EXT_OFF
      CALL GET_INTERNAL_STATES(N_EXTOP_TP,N_INTOP_TP,
     &     WORK(KLSOBEX),WORK(KL_N_INT_FOR_EXT),WORK(KL_IB_INT_FOR_EXT),
     &     WORK(KL_I_INT_FOR_EXT),WORK(KL_NDIM_IN_SE),
     &     WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &     WORK(KL_X1_INT_EI_FOR_SE), WORK(KL_X2_INT_EI_FOR_SE),
     &     WORK(KL_SG_INT_EI_FOR_SE),WORK(KL_S_INT_EI_FOR_SE),
     &     WORK(KL_IBX1_INT_EI_FOR_SE), WORK(KL_IBX2_INT_EI_FOR_SE),
     &     WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KL_IBS_INT_EI_FOR_SE),
     &     WORK(KL_X2L_INT_EI_FOR_SE),
     &     I_IN_TP,I_INT_OFF,I_EXT_OFF)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of internal states per sym and ext type'
        CALL IWRTMA(WORK(KL_N_INT_FOR_SE),
     &  NSMOB,N_EXTOP_TP, NSMOB,N_EXTOP_TP)
        WRITE(6,*) 
     &  ' Number of orthn internal states per sym and ext type'
        CALL IWRTMA(WORK(KL_N_ORTN_FOR_SE),
     &       NSMOB,N_EXTOP_TP, NSMOB,N_EXTOP_TP)
      END IF
*. Largest number of internal states for given sym and external type
C IMNNMX(IVEC,NDIM,MINMAX)
      N_INT_MAX = IMNMX(WORK(KL_N_INT_FOR_SE),N_EXTOP_TP*NSMOB,2)
*. Largest number of zero-order states of given sym and external type
      N_ORTN_MAX = IMNMX(WORK(KL_N_ORTN_FOR_SE),N_EXTOP_TP*NSMOB,2)
      WRITE(6,*) ' N_INT_MAX, N_ORTN_MAX = ', N_INT_MAX, N_ORTN_MAX
*. Largest transformation block 
      N_XEO_MAX = N_INT_MAX*N_ORTN_MAX
      IF(NTEST.GE.5) WRITE(6,*) ' Largest (EL,ORTN) block = ', N_XEO_MAX
*
*. Number of zero-order states - does now include the unit-operator
      N_ZERO_EI = N_ZERO_ORDER_STATES(WORK(KL_N_ORTN_FOR_SE),
     &            WORK(KL_NDIM_EX_ST),N_EXTOP_TP,1)
      WRITE(6,*) ' Number of zero-order states with sym 1 = ', N_ZERO_EI
C                 N_ZERO_ORDER_STATES(NORTN_FOR_SE,NDIM_EX_ST,N_EXTP,ITOTSYM)
C?    WRITE(6,*) ' Call to IDIM_TCC just after zero-order states'
C?    CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,1,
C?   &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
C?   &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
C?   &     MSCOMB_CC,MX_SBSTR,
C?   &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
C?   &     NTCONF,IPRCC)
* =============
* Scratch space
* =============
*
*. Scratch space for CI - behind the curtain
COLD   CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
*. Pointers to KVEC1 and KVEC2, transferred through GLBBAS
COLD   KVEC1P = KVEC1
COLD   KVEC2P = KVEC2
C?     WRITE(6,*) ' KVEC3 after GET_3BLKS.. ', KVEC3
*. and two CC vectors , extra element for unexcited SD at end of vectors
       N_SD_INT = 1
       LENNY = LEN_T_VEC_MX + N_SD_INT
       if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then
C!       lenny = max(lenny,nooexc(1)+nooexc(2))
         STOP ' Jeppe copied out (nooexc not defined)'
       end if
       CALL MEMMAN(KCC1,LENNY,'ADDL  ',2,'CC1_VE')
       CALL MEMMAN(KCC2,LENNY,'ADDL  ',2,'CC2_VE')
*. And the CC diagonal
       CALL MEMMAN(KDIA,LENNY,'ADDL  ',2,'CC_DIA')
*
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(NTEST.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
      IF(NTEST.GE.100) CALL PRINT_ZERO_TRMAT
*
      IF(ICTYP(1:4).EQ.'ICCI') THEN
*
*                    ==============================
*                    Internal contracted CI section 
*                    ==============================
*
        CALL LUCIA_ICCI(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                  EFINAL,CONVER,VNFINAL)
*
      ELSE IF(ICTYP(1:4).EQ.'ICPT') THEN
*
*                    ==========================================
*                    Internal contracted Perturbation expansion 
*                    ==========================================
*
        CALL LUCIA_ICPT(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                  EFINAL,CONVER,VNFINAL)
*
      ELSE IF(ICTYP(1:4).EQ.'ICCC') THEN
*
*                    ======================================
*                    Internal contracted Coupled Cluster 
*                    =======================================
*
        CALL LUCIA_ICCC(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                  EFINAL,CONVER,VNFINAL)
      END IF
*
*.
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICCI  ')
*
      RETURN
      END
      SUBROUTINE SPLIT_IE_SPOBEXTP(ISOBEX_TP, NSPOBEX_TP,
     &           N_EXTOP_TP, N_INTOP_TP,IB_EXTOP_TP,IB_INTOP_TP,
     &           N_INT_FOR_EXT, IB_INT_FOR_EXT, I_INT_FOR_EXT,
     &           NGAS, IHPVGAS)
*. A set of spinorbital excitations is given. 
*. Split these into internal and external parts
*. and store these in ISOBEX_TP
*. IHPVGAS is used to identify Internal/external parts
*
*. Jeppe Olsen, October 2006
*
*. Input :
* ISOBEX_TP : The spinorbital types in CAAB form
* NSPOBEX_TP : Number of spinorbitaltypes
* IB_EXTOP_TP : offset for storing external operators in ISOBEX_TP
*. output
* N_EXTOP_TP : number of external spinorbitalexcitation types
* N_INTOP_TP : number of internal spinorbitalexcitation types
* N_INT_FOR_EXT(IEXT) : number of internal types for external type IEXT
* I_INT_FOR_EXT() : gives the internal types for each external
* IB_INT_FOR_EXT(IEXT) : Offset for given external types in I_INT_FOR_EXT
*
* Note the assymmetry : IB_EXTOP_TP is input whereas IB_INTOP_TP is output
* 
* 
      INCLUDE 'wrkspc.inc'
*. Input (output added)
      INTEGER ISOBEX_TP(NGAS,4,*), IHPVGAS(NGAS)
*. Output
      INTEGER N_INT_FOR_EXT(*), IB_INT_FOR_EXT(*), I_INT_FOR_EXT(*)
*. Scratch : Internal and external parts
      INTEGER I_INT_TP(4*MXPNGAS)
      INTEGER I_EXT_TP(4*MXPNGAS)
*
      NTEST = 100
*
*. Obtain the various external parts
*
      IZERO = 0
      CALL ISETVC(N_INT_FOR_EXT,IZERO,NSPOBEX_TP)
      N_EXTOP_TP = 0
      DO I_EI_OP = 1, NSPOBEX_TP
*. Split operator into ext and int parts
        CALL SPLIT_EIOP(ISOBEX_TP(1,1,I_EI_OP),I_EXT_TP, I_INT_TP,
     &                  NGAS,IHPVGAS)
*. Is this a new external operator ?
        NEW = 1
        IDENT = 0
        DO J_EXTOP_TP = 1, N_EXTOP_TP
C        COMPARE_TWO_INTARRAYS(IA,IB,NAB,IDENT)
         CALL COMPARE_TWO_INTARRAYS(I_EXT_TP,
     &        ISOBEX_TP(1,1,IB_EXTOP_TP-1+J_EXTOP_TP),4*NGAS,IDENT)
         IF(IDENT.EQ.1) THEN
*. Match with previous type
           N_INT_FOR_EXT(J_EXTOP_TP) =  N_INT_FOR_EXT(J_EXTOP_TP) + 1
           NEW = 0
         END IF
        END DO
        IF(NEW.EQ.1) THEN
          N_EXTOP_TP = N_EXTOP_TP + 1
          N_INT_FOR_EXT(N_EXTOP_TP) = 1
          CALL ICOPVE(I_EXT_TP,
     &         ISOBEX_TP(1,1,IB_EXTOP_TP-1+N_EXTOP_TP),4*NGAS)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of external spinorbitalexcitation types ',
     &               N_EXTOP_TP
        WRITE(6,*) ' And the external spinorbitalexcitation types : '
        CALL WRT_SPOX_TP(ISOBEX_TP(1,1,IB_EXTOP_TP),N_EXTOP_TP)
      END IF
*. Offsets for the various internal parts for each external part
      IB_INT_FOR_EXT(1) = 1
      DO J_EXTOP_TP = 2, N_EXTOP_TP
        IB_INT_FOR_EXT(J_EXTOP_TP) =  IB_INT_FOR_EXT(J_EXTOP_TP-1)
     &                             +  N_INT_FOR_EXT(J_EXTOP_TP-1)
      END DO
C?    WRITE(6,*) ' N_INT_FOR_EXT, IB_INT_FOR_EXT '
C?    CALL IWRTMA(N_INT_FOR_EXT,1,N_EXTOP_TP,1,N_EXTOP_TP)
C?    CALL IWRTMA(IB_INT_FOR_EXT,1,N_EXTOP_TP,1,N_EXTOP_TP)
*. Offsets for internal parts
      IB_INTOP_TP = IB_EXTOP_TP+N_EXTOP_TP
C?    WRITE(6,*) ' IB_INTOP_TP, IB_EXTOP_TP, N_EXTOP_TP = ',
C?   &             IB_INTOP_TP, IB_EXTOP_TP, N_EXTOP_TP
*. And generate the various internal parts
      N_INTOP_TP = 0
      CALL ISETVC(N_INT_FOR_EXT,IZERO,NSPOBEX_TP)
      DO I_EI_OP = 1, NSPOBEX_TP
C?      WRITE(6,*) ' Info for I_EI_OP = ', I_EI_OP
*. Split operator into ext and int parts
        CALL SPLIT_EIOP(ISOBEX_TP(1,1,I_EI_OP),I_EXT_TP, I_INT_TP,
     &                  NGAS,IHPVGAS)
*. Type of this external operator
        JJ_EXTOP_TP = -1
        DO J_EXTOP_TP = 1, N_EXTOP_TP
         CALL COMPARE_TWO_INTARRAYS(I_EXT_TP,
     &        ISOBEX_TP(1,1,IB_EXTOP_TP-1+J_EXTOP_TP),4*NGAS,IDENT)
         IF(IDENT.EQ.1) JJ_EXTOP_TP = J_EXTOP_TP
        END DO
C?      WRITE(6,*) 'JJ_EXTOP_TP = ', JJ_EXTOP_TP
*. Is internal type new ?
        NEW = 1
        IDENT = 0
        DO J_INTOP_TP = 1, N_INTOP_TP
         CALL COMPARE_TWO_INTARRAYS(I_INT_TP,
     &        ISOBEX_TP(1,1,IB_INTOP_TP-1+J_INTOP_TP),4*NGAS,IDENT)
         IF(IDENT.EQ.1) THEN
*. Match with previous type
           NEW = 0
           JJ_INTOP_TP = J_INTOP_TP
         END IF
C?       WRITE(6,*) ' IDENT = ', IDENT
        END DO
C?      WRITE(6,*) ' NEW = ', NEW
        IF(NEW.EQ.1) THEN
          N_INTOP_TP = N_INTOP_TP + 1
          JJ_INTOP_TP = N_INTOP_TP
          CALL ICOPVE(I_INT_TP,ISOBEX_TP(1,1,IB_INTOP_TP-1+N_INTOP_TP),
     &                4*NGAS)
        END IF
        N_INT_FOR_EXT(JJ_EXTOP_TP) =  N_INT_FOR_EXT(JJ_EXTOP_TP) + 1
        I_INT_FOR_EXT(IB_INT_FOR_EXT(JJ_EXTOP_TP)-1
     &                +N_INT_FOR_EXT(JJ_EXTOP_TP)) = JJ_INTOP_TP
C?      WRITE(6,*) ' IB_INT_FOR_EXT(JJ_EXTOP_TP)-1 + ... ',
C?   &  IB_INT_FOR_EXT(JJ_EXTOP_TP)-1 + N_INT_FOR_EXT(JJ_EXTOP_TP)
C?      WRITE(6,*) ' JJ_EXTOP_TP = ', JJ_EXTOP_TP
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of internal excitationtypes ', N_INTOP_TP
        WRITE(6,*) ' And the internal excitationtypes : '
        CALL WRT_SPOX_TP(ISOBEX_TP(1,1,IB_INTOP_TP),N_INTOP_TP)
        WRITE(6,*) ' And the various EI types for each E-type'
        DO J_EXTOP_TP = 1, N_EXTOP_TP
         WRITE(6,*) ' External operator type ', J_EXTOP_TP
         WRITE(6,*) ' Number of internal types for this ', 
     &   N_INT_FOR_EXT(J_EXTOP_TP)
         WRITE(6,*) ' and the internal operator types '
         LEN = N_INT_FOR_EXT(J_EXTOP_TP)
         CALL IWRTMA(I_INT_FOR_EXT(IB_INT_FOR_EXT(J_EXTOP_TP)),
     &               1,LEN,1,LEN)
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE SPLIT_EIOP(I_EIOP,I_EOP,I_IOP,NGAS,IHPVGAS)
*
* Split operator IEOP into external and internal parts 
* according to IHPVGAS
*
*. Jeppe Olsen, October 2006
*
      INCLUDE 'wrkspc.inc'
*. Input
      INTEGER I_EIOP(NGAS,4), IHPVGAS(NGAS)
*. Output
      INTEGER I_IOP(NGAS,4),I_EOP(NGAS,4)
*
      IZERO = 0
      CALL ISETVC(I_EOP,IZERO,4*NGAS)
      CALL ISETVC(I_IOP,IZERO,4*NGAS)
*
      DO IGAS = 1, NGAS
       IF(IHPVGAS(IGAS).LE.2) THEN
*. External : Secondary or inactive
          DO ICAAB = 1, 4
            I_EOP(IGAS,ICAAB) = I_EIOP(IGAS,ICAAB)
          END DO
        ELSE 
          DO ICAAB = 1, 4
            I_IOP(IGAS,ICAAB) = I_EIOP(IGAS,ICAAB)
          END DO
        END IF
       END DO
*
       NTEST = 00
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Splitting EI operator : '
         WRITE(6,*) ' Input EI and output E, I '
         CALL WRT_SPOX_TP(I_EIOP,1)
         CALL WRT_SPOX_TP(I_EOP,1)
         CALL WRT_SPOX_TP(I_IOP,1)
       END IF
*
       RETURN
       END 
       SUBROUTINE DIMENSION_EI_EXP(ISPOBEX_TP,IB_EXTP,IB_INTP,N_EXTP,
     &            N_INTP,N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,
     /            NDIM_EX_ST,NDIM_IN_ST,NDIM_EI,NDIM_IN_SE,NSYM)
*
* NOTE: Symmetry of EI is assumed to be 1 (totsym)
*
* Obtain dimension for EI expansion
* ISPOBEX_TP contains all spinorbital excitations with external types
* starting at IB_EXTP and internal types starting at IB_INTP.
* The combinations of internal and external types are specified 
* by  N_IN_FOR_EX, IB_IN_FOR_EX, I_IN_FOR_EX.
*
*. On output
*     NDIM_IN_ST(ISYM,INTP) : Number of internal strings of sym ISYM and type INTP
*     NDIM_EX_ST(ISYM,INTP) : Number of external strings of sym ISYM and type INTP
*     NDIM_IN_SE(ISYM,IEXTP) : Number of internal strings per symmetry and extenal type
* N_DIM_EI : Dimension of complete expansion
*
*. Jeppe Olsen, October 2006
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cprnt.inc'
*. Input
      INTEGER ISPOBEX_TP(NGAS,4,*)
      INTEGER  N_IN_FOR_EX(N_EXTP),IB_IN_FOR_EX(N_EXTP),
     &         I_IN_FOR_EX(N_EXTP)
*. Local scratch
      INTEGER ISCR1(MXPSTT),ISCR2(MXPSTT)
*. Output
      INTEGER NDIM_EX_ST(NSYM,N_EXTP), NDIM_IN_ST(NSYM,N_INTP)
      INTEGER NDIM_IN_SE(NSYM,N_EXTP)
*
      NTEST = 100
      NTEST = MAX(NTEST,IPRCC)
*
*. External parts 
      MX_ST_TSOSO_EX = 0
      MX_ST_TSOSO_BLK_EX = 0
      MX_TBLK_EX = 0
      MX_TBLK_AS_EX = 0
*
      DO ISYM = 1, NSYM
C     IDIM_TCC(ITSOSO_TP,NTSOSO_TP,ISYM,
C    &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
C    &           LTSOSO_TP,IBTSOSO_TP,IDIM_T,MSCOMB_CC,
C    &           MX_TBLK_AS,ISPOX_FOR_OCCLS,NOCCLS,IBSPOX_FOR_OCCLS,
C    &           NTCONF,IPRCC)
*. External parts
        CALL IDIM_TCC(ISPOBEX_TP(1,1,IB_EXTP),N_EXTP,ISYM,
     &       MX_ST_TSOSO_EXL,MX_ST_TSOSO_BLK_EXL,MX_TBLK_EXL,
     &       ISCR1, ISCR2, IDIM_T_EXL, 0,
     &       MX_TBLK_AS_EXL,IDUM,0,IDIM,IDUM,IPRCC)
        MX_ST_TSOSO_EX = MAX(MX_ST_TSOSO_EX,MX_ST_TSOSO_EXL)
        DO I_EXTP = 1, N_EXTP
          NDIM_EX_ST(ISYM,I_EXTP) = ISCR1(I_EXTP)
        END DO
        MX_ST_TSOSO_BLK_EX = MAX(MX_ST_TSOSO_BLK_EX,MX_ST_TSOSO_BLK_EXL)
        MX_TBLK_EX = MAX(MX_TBLK_EX,MX_TBLK_EXL)
        MX_TBLK_AS_EX = MAX(MX_TBLK_AS_EX,MX_TBLK_AS_EXL)
      END DO
*. Internal parts 
      MX_ST_TSOSO_IN = 0
      MX_ST_TSOSO_BLK_IN = 0
      MX_TBLK_IN = 0
      MX_TBLK_AS_IN = 0
*
      DO ISYM = 1, NSYM
        CALL IDIM_TCC(ISPOBEX_TP(1,1,IB_INTP),N_INTP,ISYM,
     &       MX_ST_TSOSO_INL,MX_ST_TSOSO_BLK_INL,MX_TBLK_INL,
     &       ISCR1, ISCR2, IDIM_T_INL, 0,
     &       MX_TBLK_AS_INL,IDUM,0,IDIM,IDUM,IPRCC)
        MX_ST_TSOSO_IN = MAX(MX_ST_TSOSO_IN,MX_ST_TSOSO_INL)
        DO I_INTP = 1, N_INTP
          NDIM_IN_ST(ISYM,I_INTP) = ISCR1(I_INTP)
        END DO
        MX_ST_TSOSO_BLK_IN = MAX(MX_ST_TSOSO_BLK_IN,MX_ST_TSOSO_BLK_INL)
        MX_TBLK_IN = MAX(MX_TBLK_IN,MX_TBLK_INL)
        MX_TBLK_AS_IN = MAX(MX_TBLK_AS_IN,MX_TBLK_AS_INL)
      END DO
*. Obtain number of internals per symmetry and external.
      DO I_EXTP =  1, N_EXTP
       IB = IB_IN_FOR_EX(I_EXTP)
       N_IN = N_IN_FOR_EX(I_EXTP)
       DO ISYM = 1, NSYM
        NDIM_IN_SE(ISYM,I_EXTP) = 0
        DO II_INTP = 1, N_IN
         I_INTP = I_IN_FOR_EX(IB-1+II_INTP)
         NDIM_IN_SE(ISYM,I_EXTP) =  
     &   NDIM_IN_SE(ISYM,I_EXTP) + NDIM_IN_ST(ISYM,I_INTP)
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Dimension of externals: '
        WRITE(6,*) ' ======================== '
        DO I_EXTP = 1, N_EXTP
         WRITE(6,*) ' External type  : ', I_EXTP
         CALL IWRTMA(NDIM_EX_ST(1,I_EXTP),1,NSYM,1,NSYM)
        END DO
        WRITE(6,*) ' Dimension of internals: '
        WRITE(6,*) ' ======================== '
        DO I_INTP = 1, N_INTP
         WRITE(6,*) ' Internal type  : ', I_INTP
         CALL IWRTMA(NDIM_IN_ST(1,I_INTP),1,NSYM,1,NSYM)
        END DO
        WRITE(6,*) ' Number of internals per symmetry and external:'
        WRITE(6,*) ' ============================================== '
        CALL IWRTMA(NDIM_IN_SE,NSYM,N_EXTP,NSYM,N_EXTP)
      END IF
*. And then for the various combinations of internals and externals for sym1 
*. (for comparison with standard test)
      NDIM_EI = 0
      DO I_EXTP = 1, N_EXTP
       DO II_INTP = 1, N_IN_FOR_EX(I_EXTP)
         I_INTP = I_IN_FOR_EX(IB_IN_FOR_EX(I_EXTP)-1+II_INTP)
         DO ISYM = 1, NSYM
           NDIM_EI = NDIM_EI 
     &   + NDIM_EX_ST(ISYM,I_EXTP)*NDIM_IN_ST(ISYM,I_INTP)
         END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of EI combinations with sym 1 = ', NDIM_EI
      END IF
*
      RETURN
      END
      SUBROUTINE EI_REORDER_ARRAYS(IREO_EI_ST,
     &           ISPOBEX_TP,IB_EXTP,IB_INTP,
     &           L_SPOBEX_TP,IB_SPOBEX_TP,NSPOBEX_TPL,
     &           N_EXTP,N_INTP,
     &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,ITSYM)
*
* Generate reorder array for going between EI and standard CAAB type
* orders
* Jeppe Olsen, October 2006, modified March 2009: Changed into 
* T(Internal, external ordering), with all internal and external 
* strings corresponding to given type of external string and 
* symmetry of internal string(sic) are a single matrix block
*. Total symmetry of operator is ITSYM
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'ctcc.inc'
*. Some of the dimensions.
      INTEGER N_IN_FOR_EX(N_EXTP)
      
*
      IDUM = 0
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NELMAX = MAX(NAEL,NBEL)
*. largest number of external types for given internal type
      MX_INTP = IMNMX(N_IN_FOR_EX,N_EXTP,2)
C?    WRITE(6,*) ' MX_INTP = ', MX_INTP
      IF(MX_INTP.GT.MXP_NINTP_FOR_EX) THEN
        WRITE(6,*) ' Problem with fixed dimension '
        WRITE(6,*) ' Observed in routine EI_REORDER_ARRAYS'
        WRITE(6,*)  ' MXP_NINTP_FOR_EX: actual and required value',
     &               MX_INTP,MXP_NINTP_FOR_EX
        WRITE(6,*) ' Increase value of MXP_NINTP_FOR_EX'
        STOP       ' Increase value of MXP_NINTP_FOR_EX'
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',1,'EI_REO')
*. 
      CALL MEMMAN(KL_IST_SCR,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'I_STR_')
      LEN_Z = NELMAX*NTOOB
      CALL MEMMAN(KL_IZ_CA,LEN_Z*MX_INTP,'ADDL  ',1,'IZ_CA ')
      CALL MEMMAN(KL_IZ_CB,LEN_Z*MX_INTP,'ADDL  ',1,'IZ_CB ')
      CALL MEMMAN(KL_IZ_AA,LEN_Z*MX_INTP,'ADDL  ',1,'IZ_AA ')
      CALL MEMMAN(KL_IZ_AB,LEN_Z*MX_INTP,'ADDL  ',1,'IZ_AB ')
      L_IZ_SCR = (NELMAX+1)*(NTOOB+1) + 2 * NTOOB
      CALL MEMMAN(KL_IZ_SCR,L_IZ_SCR,'ADDL  ',1,'IZ_SCR')
*
      LEN_REO= NSMOB*MX_ST_TSOSO_MX
      CALL MEMMAN(KL_IREO_CA,MX_INTP*LEN_REO,'ADDL  ',1,'IRE_CA')
      CALL MEMMAN(KL_IREO_CB,MX_INTP*LEN_REO,'ADDL  ',1,'IRE_CB')
      CALL MEMMAN(KL_IREO_AA,MX_INTP*LEN_REO,'ADDL  ',1,'IRE_AA')
      CALL MEMMAN(KL_IREO_AB,MX_INTP*LEN_REO,'ADDL  ',1,'IRE_AB')
*.
      LEN_STBLK = NELMAX*MX_ST_TSOSO_BLK_MX
      CALL MEMMAN(KL_IST_EX_CA,LEN_STBLK*NSMOB,'ADDL  ',1,'STE_CA')
      CALL MEMMAN(KL_IST_EX_CB,LEN_STBLK*NSMOB,'ADDL  ',1,'STE_CB')
      CALL MEMMAN(KL_IST_EX_AA,LEN_STBLK*NSMOB,'ADDL  ',1,'STE_AA')
      CALL MEMMAN(KL_IST_EX_AB,LEN_STBLK*NSMOB,'ADDL  ',1,'STE_AB')
*
      CALL MEMMAN(KL_IST_IN_CA,LEN_STBLK*NSMOB*MX_INTP,'ADDL  ',1,
     &            'STI_CA')
      CALL MEMMAN(KL_IST_IN_CB,LEN_STBLK*NSMOB*MX_INTP,'ADDL  ',1,
     &            'STI_CB')
      CALL MEMMAN(KL_IST_IN_AA,LEN_STBLK*NSMOB*MX_INTP,'ADDL  ',1,
     &            'STI_AA')
      CALL MEMMAN(KL_IST_IN_AB,LEN_STBLK*NSMOB*MX_INTP,'ADDL  ',1,
     &            'STI_AB')
*. Obtain array for reordering types from EI to standard
C?    WRITE(6,*) ' IB_EXTP, IB_INTP =', IB_EXTP, IB_INTP
      CALL MEMMAN(KL_EI_TP_REO,NSPOBEX_TP+1,'ADDL  ',1,'EI_TPR')
      CALL EITP_TO_SPOXTP_REO(WORK(KL_EI_TP_REO),ISPOBEX_TP,
     &      NSPOBEX_TPL,IB_EXTP,IB_INTP,NGAS,N_EXTP,N_INTP,
     /      N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX) 
C     EITP_TO_SPOXTP_REO(I_EI_TP_REO,ISPOBEX_TP,NSPOBEX_TP,
C    &           IB_EXTP, IB_INTP,NGAS,N_EXTP,N_INTP,
C    &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX)
      CALL MEMCHK2('AFT_TP')
*
      LCHECK = 123456789
      CALL EI_REORDER_ARRAYS_S(IREO_EI_ST,
     &     ISPOBEX_TP,IB_EXTP,IB_INTP,NGAS,
     &     L_SPOBEX_TP,IB_SPOBEX_TP,NOBPT,NTOOB,
     &     NSPOBEX_TPL,N_EXTP,N_INTP,
     &     N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,
     &     WORK(KL_EI_TP_REO),ITSYM,NSMOB,
     &     WORK(KL_IST_SCR),WORK(KL_IZ_CA),WORK(KL_IZ_CB),
     /     WORK(KL_IZ_AA),WORK(KL_IZ_AB),WORK(KL_IZ_SCR),
     /     WORK(KL_IREO_CA),WORK(KL_IREO_CB),
     &     WORK(KL_IREO_AA),WORK(KL_IREO_AB),
     &     WORK(KL_IST_EX_CA),WORK(KL_IST_EX_CB),
     &     WORK(KL_IST_EX_AA),WORK(KL_IST_EX_AB),
     &     WORK(KL_IST_IN_CA),WORK(KL_IST_IN_CB),
     &     WORK(KL_IST_IN_AA),WORK(KL_IST_IN_AB),LEN_Z,LEN_REO,
     &     LEN_STBLK,LCHECK)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'EI_REO')
      RETURN
      END 
      SUBROUTINE EI_REORDER_ARRAYS_S(IREO_EI_ST,
     &           ISPOBEX_TP,IB_EXTP, IB_INTP,NGAS,
     &           L_SPOBEX_TP,IB_SPOBEX_TP,NOBPT,NTOOB,
     &           NSPOBEX_TP,N_EXTP,N_INTP,
     &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,
     &           I_EI_TP_REO,ITSYM,NSMOB,
     &           IST_SCR,IZ_CA,IZ_CB,IZ_AA,IZ_AB,IZ_SCR,
     &           IREO_CA,IREO_CB,IREO_AA,IREO_AB,
     &           IST_EX_CA,IST_EX_CB,IST_EX_AA,IST_EX_AB,
     &           IST_IN_CA,IST_IN_CB,IST_IN_AA,IST_IN_AB,LEN_Z,LEN_REO,
     &           LEN_STBLK,LCHECK)
*. Obtain ordering : EI-ordering => standard-ordering
*
*. Jeppe Olsen, October 2006, modified march 2009
*
* Obtain reorder array going from standard CAAB ordering into
* EI ordering. The EI order is T(internal, external)
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER ISPOBEX_TP(NGAS,4,*)
      INTEGER NOBPT(NGAS)
*. Length and offset for each spinorbitalexcitation
      INTEGER L_SPOBEX_TP(NSPOBEX_TP),IB_SPOBEX_TP(NSPOBEX_TP)
      INTEGER N_IN_FOR_EX(N_EXTP),IB_IN_FOR_EX(N_EXTP),I_IN_FOR_EX(*)
      INTEGER I_EI_TP_REO(NSPOBEX_TP)
*. Output : reorder array with signs included
      INTEGER IREO_EI_ST(*)
*. Scratch through argument list
*. Space for occupations of strings 
      INTEGER IST_SCR(*)
*. Z arrays for lexical ordering and scratch array for constructing then
      INTEGER IZ_CA(LEN_Z,*),IZ_CB(LEN_Z,*),
     &        IZ_AA(LEN_Z,*),IZ_AB(LEN_Z,*),IZ_SCR(*)

*. Reorder arrays
      INTEGER IREO_CA(LEN_REO,*), IREO_CB(LEN_REO,*),
     &        IREO_AA(LEN_REO,*), IREO_AB(LEN_REO,*)
*. Occupation of external and internal strings
      INTEGER IST_EX_CA(LEN_STBLK,NSMOB), IST_EX_CB(LEN_STBLK,NSMOB)
      INTEGER IST_EX_AA(LEN_STBLK,NSMOB), IST_EX_AB(LEN_STBLK,NSMOB)
      INTEGER IST_IN_CA(LEN_STBLK,NSMOB,*), IST_IN_CB(LEN_STBLK,NSMOB,*)
      INTEGER IST_IN_AA(LEN_STBLK,NSMOB,*), IST_IN_AB(LEN_STBLK,NSMOB,*)
*. Local scratch ( DIMENSIONS ???? )
      INTEGER I_EI_CA(MXPORB),I_EI_CB(MXPORB),I_EI_AA(MXPORB),
     &        I_EI_AB(MXPORB)
      INTEGER IGRP(MXPNGAS)
*
      INTEGER NEL_IN_CA(MXP_NINTP_FOR_EX),NEL_IN_CB(MXP_NINTP_FOR_EX)
      INTEGER NEL_IN_AA(MXP_NINTP_FOR_EX),NEL_IN_AB(MXP_NINTP_FOR_EX)
      INTEGER NEL_CA(MXP_NINTP_FOR_EX),NEL_CB(MXP_NINTP_FOR_EX)
      INTEGER NEL_AA(MXP_NINTP_FOR_EX),NEL_AB(MXP_NINTP_FOR_EX)
      INTEGER I_ST_TP(MXP_NINTP_FOR_EX)
      INTEGER IB_ST_TP(MXP_NINTP_FOR_EX)
      INTEGER LEN_ST(MXP_NINTP_FOR_EX)
      INTEGER IB_ST_SSS(8,8,8,MXP_NINTP_FOR_EX)
*
      INTEGER NST_CA(8,MXP_NINTP_FOR_EX),
     &        NST_CB(8,MXP_NINTP_FOR_EX),
     &        NST_AA(8,MXP_NINTP_FOR_EX),
     &        NST_AB(8,MXP_NINTP_FOR_EX)
*
      INTEGER NST_EX_CA(8),NST_EX_CB(8),NST_EX_AA(8),NST_EX_AB(8)
*
      INTEGER NST_IN_CA(8,MXP_NINTP_FOR_EX),
     &        NST_IN_CB(8,MXP_NINTP_FOR_EX),
     &        NST_IN_AA(8,MXP_NINTP_FOR_EX),
     &        NST_IN_AB(8,MXP_NINTP_FOR_EX) 

*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ------------------------------------------------'
        WRITE(6,*) ' Information from subroutine EI_REORDER_ARRAYS_S '
        WRITE(6,*) ' ------------------------------------------------'
        WRITE(6,*)
      END IF
C?    WRITE(6,*) 'LEN_STBLK = ', LEN_STBLK
C?    WRITE(6,*) ' I_EI_TP_REO at start of EI_REORDER' 
C?    CALL IWRTMA(I_EI_TP_REO,1,NSPOBEX_TP,1,NSPOBEX_TP)
*
      I_EI_TP = 0
      I_EI_OP = 0
      DO J_EXTP = 1, N_EXTP
        J_EXTP_ABS = IB_EXTP-1+J_EXTP
*
        NEL_EX_CA = IELSUM(ISPOBEX_TP(1,1,J_EXTP_ABS),NGAS)
        NEL_EX_CB = IELSUM(ISPOBEX_TP(1,2,J_EXTP_ABS),NGAS)
        NEL_EX_AA = IELSUM(ISPOBEX_TP(1,3,J_EXTP_ABS),NGAS)
        NEL_EX_AB = IELSUM(ISPOBEX_TP(1,4,J_EXTP_ABS),NGAS)
*
        N_INTP_J = N_IN_FOR_EX(J_EXTP) 
        IF(NTEST.GE.100) THEN
          WRITE(6,*)
          WRITE(6,*) ' Information for external type = ', J_EXTP
          WRITE(6,*) ' Number of internal types =', N_INTP_J
          WRITE(6,*) ' I_EI_OP at start ', I_EI_OP
          WRITE(6,*)
        END IF
        IB_INTP_J = IB_IN_FOR_EX(J_EXTP)
*
*. Construct information for the various internal and combined ei types
*. that will be used for this external type
        DO JJ_INTP = 1, N_INTP_J
          I_EI_TP = I_EI_TP + 1
          J_INTP = I_IN_FOR_EX(IB_INTP_J-1+JJ_INTP)
          J_INTP_ABS = IB_INTP-1+J_INTP
*. Corresponding type in standard order
          I_ST_TP(JJ_INTP) = I_EI_TP_REO(I_EI_TP)
          IF(NTEST.GE.10000)
     &    WRITE(6,*) ' JJ_INTP, I_EI_TP, I_ST, I_EI ',
     &    JJ_INTP,I_EI_TP,I_EI_TP_REO(I_EI_TP),I_ST_TP(JJ_INTP)
          IB_ST_TP(JJ_INTP) = IB_SPOBEX_TP(I_ST_TP(JJ_INTP))
C?        WRITE(6,*) ' IB_ST_TP = ', IB_ST_TP
*
          IF(NTEST.GE.10000) THEN
            WRITE(6,'(A,4I4)') ' J_INTP, IB_INTP, J_INTP_ABS =', 
     &                           J_INTP, IB_INTP, J_INTP_ABS
          END IF
*. Number of operators in internal part
          NEL_IN_CA(JJ_INTP) = IELSUM(ISPOBEX_TP(1,1,J_INTP_ABS),NGAS)
          NEL_IN_CB(JJ_INTP) = IELSUM(ISPOBEX_TP(1,2,J_INTP_ABS),NGAS)
          NEL_IN_AA(JJ_INTP) = IELSUM(ISPOBEX_TP(1,3,J_INTP_ABS),NGAS)
          NEL_IN_AB(JJ_INTP) = IELSUM(ISPOBEX_TP(1,4,J_INTP_ABS),NGAS)
*. Number of operators in combined  external/internal part
          NEL_CA(JJ_INTP) =IELSUM(ISPOBEX_TP(1,1,I_ST_TP(JJ_INTP)),NGAS)
          NEL_CB(JJ_INTP) =IELSUM(ISPOBEX_TP(1,2,I_ST_TP(JJ_INTP)),NGAS)
          NEL_AA(JJ_INTP) =IELSUM(ISPOBEX_TP(1,3,I_ST_TP(JJ_INTP)),NGAS)
          NEL_AB(JJ_INTP) =IELSUM(ISPOBEX_TP(1,4,I_ST_TP(JJ_INTP)),NGAS)
*
          IF(NTEST.GE.10000) THEN
            WRITE(6,'(A,4I4)') ' << E, I, EI, ST types >> ', 
     &      J_EXTP, J_INTP, I_EI_TP, I_ST_TP(JJ_INTP)
          END IF
          IF(NTEST.GE.10000) THEN
            WRITE(6,'(A,4I4)') 'NEL_CA, NEL_CB, NEL_AA, NEL_AB = ',
     &      NEL_CA(JJ_INTP), NEL_CB(JJ_INTP), 
     &      NEL_AA(JJ_INTP), NEL_AB(JJ_INTP)
          END IF
*. Construct reorder arrays for the combined CA, CB, AA, AB strings
C  WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,ISCR,NTEST)
          CALL WEIGHT_SPGP(IZ_CA(1,JJ_INTP),NGAS,
     &         ISPOBEX_TP(1,1,I_ST_TP(JJ_INTP)),NOBPT,IZ_SCR,0)
          CALL WEIGHT_SPGP(IZ_CB(1,JJ_INTP),NGAS,
     &         ISPOBEX_TP(1,2,I_ST_TP(JJ_INTP)),NOBPT,IZ_SCR,0)
          CALL WEIGHT_SPGP(IZ_AA(1,JJ_INTP),NGAS,
     &         ISPOBEX_TP(1,3,I_ST_TP(JJ_INTP)),NOBPT,IZ_SCR,0)
          CALL WEIGHT_SPGP(IZ_AB(1,JJ_INTP),NGAS,
     &         ISPOBEX_TP(1,4,I_ST_TP(JJ_INTP)),NOBPT,IZ_SCR,0)
          DO ISYM = 1, NSMOB
*. CA
            CALL OCC_TO_GRP(ISPOBEX_TP(1,1,I_ST_TP(JJ_INTP)),IGRP,1)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,ISYM,NEL_CA(JJ_INTP),
     &           NST_CA_L,IST_SCR,NTOOB,1,IZ_CA(1,JJ_INTP),
     &           IREO_CA(1,JJ_INTP))
            NST_CA(ISYM,JJ_INTP) = NST_CA_L
*. CB
            CALL OCC_TO_GRP(ISPOBEX_TP(1,2,I_ST_TP(JJ_INTP)),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,ISYM,NEL_CB(JJ_INTP),
     &           NST_CB_L,IST_SCR,NTOOB,1,IZ_CB(1,JJ_INTP),
     &           IREO_CB(1,JJ_INTP))
            NST_CB(ISYM,JJ_INTP) = NST_CB_L

*. AA
            CALL OCC_TO_GRP(ISPOBEX_TP(1,3,I_ST_TP(JJ_INTP)),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,ISYM,NEL_AA(JJ_INTP),
     &           NST_AA_L,IST_SCR,NTOOB,1,IZ_AA(1,JJ_INTP),
     &           IREO_AA(1,JJ_INTP))
            NST_AA(ISYM,JJ_INTP) = NST_AA_L
*. AB
            CALL OCC_TO_GRP(ISPOBEX_TP(1,4,I_ST_TP(JJ_INTP)),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,ISYM,NEL_AB(JJ_INTP),
     &           NST_AB_L,IST_SCR,NTOOB,1,IZ_AB(1,JJ_INTP),
     &           IREO_AB(1,JJ_INTP))
            NST_AB(ISYM,JJ_INTP) = NST_AB_L
*
*. And the internal strings
*. CA
            CALL OCC_TO_GRP(ISPOBEX_TP(1,1,J_INTP_ABS),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &           ISYM,NEL_IN_CA(JJ_INTP),NST_IN_CA(ISYM,JJ_INTP),
     &           IST_IN_CA(1,ISYM,JJ_INTP),NTOOB,0,
     &           IDUM,IDUM)
*. CB
            CALL OCC_TO_GRP(ISPOBEX_TP(1,2,J_INTP_ABS),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &           ISYM,NEL_IN_CB(JJ_INTP),NST_IN_CB(ISYM,JJ_INTP),
     &           IST_IN_CB(1,ISYM,JJ_INTP),NTOOB,0,
     &           IDUM,IDUM)
*. AA
            CALL OCC_TO_GRP(ISPOBEX_TP(1,3,J_INTP_ABS),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &           ISYM,NEL_IN_AA(JJ_INTP),NST_IN_AA(ISYM,JJ_INTP),
     &           IST_IN_AA(1,ISYM,JJ_INTP),NTOOB,0,
     &           IDUM,IDUM)
C?          WRITE(6,*) ' The internal AA strings right after delivery'
C?          CALL IWRTMA(IST_IN_AA(1,ISYM,JJ_INTP),NEL_IN_AA(JJ_INTP),
C?   &         NST_IN_AA(1,JJ_INTP) ,NEL_IN_AA(JJ_INTP),
C?   &         NST_IN_AA(1,J_INTPJ))
*. AB 
            CALL OCC_TO_GRP(ISPOBEX_TP(1,4,J_INTP_ABS),IGRP,1)
            CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &           ISYM,NEL_IN_AB(JJ_INTP),NST_IN_AB(ISYM,JJ_INTP),
     &           IST_IN_AB(1,ISYM,JJ_INTP),NTOOB,0,
     &           IDUM,IDUM)
*

          END DO
*         ^ End of loops over symmetries
*. Array for accessing CA,CB,AA,AB blocks of given sym
          IF(NTEST.GE.10000) THEN
            WRITE(6,*) ' NST_CA, NST_CB, NST_AA, NST_AB:'
            CALL IWRTMA(NST_CA(1,JJ_INTP),1,NSMOB,1,NSMOB)
            CALL IWRTMA(NST_CB(1,JJ_INTP),1,NSMOB,1,NSMOB)
            CALL IWRTMA(NST_AA(1,JJ_INTP),1,NSMOB,1,NSMOB)
            CALL IWRTMA(NST_AB(1,JJ_INTP),1,NSMOB,1,NSMOB)
          END IF
*
          CALL Z_TCC_OFF2(IB_ST_SSS(1,1,1,JJ_INTP),LEN_ST(JJ_INTP),
     &         NST_CA(1,JJ_INTP),NST_CB(1,JJ_INTP),
     &         NST_AA(1,JJ_INTP),NST_AB(1,JJ_INTP),ITSYM,NSMOB)
        END DO
*       ^ End of loop over internal types for given external type
C?      WRITE(6,*) ' The NEL_IN_AA array in a test'
C?      CALL IWRTMA(NEL_IN_AA,1, N_INTP_J,1, N_INTP_J)
C?      WRITE(6,*) ' The NST_IN_AA-array in a test'
C?      CALL IWRTMA(NST_IN_AA,1, N_INTP_J,8, N_INTP_J)
C?      WRITE(6,*) ' The IN_AA strings in sym 1 '
C?      DO JX = 1, N_INTP_J
C?        CALL IWRTMA(IST_IN_AA(1,1,JX),NEL_IN_AA(JX),
C?   &         NST_IN_AA(1,JX) ,NEL_IN_AA(JX),
C?   &         NST_IN_AA(1,JX))
C?      END DO
*. Occupation of external strings
*. CA
        DO ISYM = 1, NSMOB
           CALL OCC_TO_GRP(ISPOBEX_TP(1,1,J_EXTP_ABS),IGRP,1)
           CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &          ISYM,NEL_EX_CA,NST_EX_CA(ISYM),IST_EX_CA(1,ISYM),
     &          NTOOB,0,IDUM,IDUM)
*. CB
           CALL OCC_TO_GRP(ISPOBEX_TP(1,2,J_EXTP_ABS),IGRP,1)
           CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &          ISYM,NEL_EX_CB,NST_EX_CB(ISYM),IST_EX_CB(1,ISYM),
     &          NTOOB,0,IDUM,IDUM)
*. AA
           CALL OCC_TO_GRP(ISPOBEX_TP(1,3,J_EXTP_ABS),IGRP,1)
           CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &          ISYM,NEL_EX_AA,NST_EX_AA(ISYM),IST_EX_AA(1,ISYM),
     &          NTOOB,0,IDUM,IDUM)
*. AB
          CALL OCC_TO_GRP(ISPOBEX_TP(1,4,J_EXTP_ABS),IGRP,1)
          CALL GETSTR2_TOTSM_SPGP(IGRP,NGAS,
     &         ISYM,NEL_EX_AB,NST_EX_AB(ISYM),IST_EX_AB(1,ISYM),
     &          NTOOB,0,IDUM,IDUM)
        END DO
C?      WRITE(6,*) ' End of external '
C?      WRITE(6,*) ' End of combined '
*. Loop over the individual strings of this External type
*  and obtain strings in EI-order
*. Loop over strings in EI order
        DO I_EX_SM = 1, NSMOB
         I_IN_SM = MULTD2H(I_EX_SM,ITSYM)
*. Loop over symmetries of external operators
         DO I_EX_C_SM = 1, NSMOB
         I_EX_A_SM = MULTD2H(I_EX_C_SM,I_EX_SM)
         DO I_EX_CA_SM = 1, NSMOB
         I_EX_CB_SM = MULTD2H(I_EX_CA_SM,I_EX_C_SM)
         DO I_EX_AA_SM = 1, NSMOB
         I_EX_AB_SM = MULTD2H(I_EX_AA_SM,I_EX_A_SM)
*. Loop over External strings (CA, CB, AA, AB)
          DO I_EX_AB = 1, NST_EX_AB(I_EX_AB_SM)
          DO I_EX_AA = 1, NST_EX_AA(I_EX_AA_SM)
          DO I_EX_CB = 1, NST_EX_CB(I_EX_CB_SM)
          DO I_EX_CA = 1, NST_EX_CA(I_EX_CA_SM)
           IF(NTEST.GE.10000) THEN
            WRITE(6,*) ' I_EX_AB, I_EX_AA, I_EX_CB, I_EX_CA = ',
     &                   I_EX_AB, I_EX_AA, I_EX_CB, I_EX_CA
           END IF
*. Loop over the internal types for this external types
           DO JJ_INTP = 1, N_INTP_J
            J_INTP = I_IN_FOR_EX(IB_INTP_J-1+JJ_INTP)
            J_INTP_ABS = IB_INTP-1+J_INTP
*. Loop over symmetries of internal operators 
            DO I_IN_C_SM = 1, NSMOB
            I_IN_A_SM = MULTD2H(I_IN_C_SM,I_IN_SM)
            DO I_IN_CA_SM = 1, NSMOB
            I_IN_CB_SM = MULTD2H(I_IN_C_SM,I_IN_CA_SM)
            DO I_IN_AA_SM = 1, NSMOB
             I_IN_AB_SM = MULTD2H(I_IN_A_SM,I_IN_AA_SM)
             IF(NTEST.GE.10000) THEN
              WRITE(6,'(A,4I4)') 
     &        'I_IN_A_SM, I_IN_CB_SM,I_IN_AB_SM,I_IN_AB_SM = ',
     &         I_IN_A_SM, I_IN_CB_SM,I_IN_AB_SM,I_IN_AB_SM  
             END IF
C?           WRITE(6,*) ' I_IN_C_SM, I_IN_A_SM = ',
C?   &                    I_IN_C_SM, I_IN_A_SM 
             IF(NTEST.GE.10000) THEN
               WRITE(6,'(A,4I4)') 
     &         'I_IN_CA_SM, I_IN_CB_SM, I_IN_AA_SM, I_IN_AB_SM = ',
     &          I_IN_CA_SM, I_IN_CB_SM, I_IN_AA_SM, I_IN_AB_SM
               WRITE(6,'(A,4I4)') 
     &         'I_EX_CA_SM, I_EX_CB_SM, I_EX_AA_SM, I_EX_AB_SM = ',
     &          I_EX_CA_SM, I_EX_CB_SM, I_EX_AA_SM, I_EX_AB_SM
               WRITE(6,'(A,4I4)') 
     &         'I_EX_CA_SM, I_EX_CB_SM, I_EX_AA_SM, I_EX_AB_SM = ',
     &          I_EX_CA_SM, I_EX_CB_SM, I_EX_AA_SM, I_EX_AB_SM
             END IF
*.  Symmetry of complete CA,CB,AA and AB strings
*.  Symmetry of complete CA,CB,AA and AB strings
             I_CA_SM = MULTD2H(I_EX_CA_SM,I_IN_CA_SM)
             I_CB_SM = MULTD2H(I_EX_CB_SM,I_IN_CB_SM)
             I_AA_SM = MULTD2H(I_EX_AA_SM,I_IN_AA_SM)
             I_AB_SM = MULTD2H(I_EX_AB_SM,I_IN_AB_SM)
* Sign to bring ECA ECB EAA EAB ICA ICB IAA IAB into
*            ECA ICA ECB ICB EAA IAB EAB IAB
             NPERMG = NEL_EX_AB*
     &       (NEL_IN_CA(JJ_INTP)+NEL_IN_CB(JJ_INTP)+NEL_IN_AA(JJ_INTP))
     &      +         NEL_EX_AA*
     &       (NEL_IN_CA(JJ_INTP)+NEL_IN_CB(JJ_INTP))
     &      +         NEL_EX_CB*
     &       NEL_IN_CA(JJ_INTP)
             IF(MOD(NPERMG,2).EQ.0) THEN
              ISIGNG = 1
             ELSE
              ISIGNG = -1
             END IF
*
             IF(NTEST.GE.10000) THEN
               WRITE(6,'(A,4I4)') 
     &         'I_CA_SM, I_CB_SM, I_AA_SM, I_AB_SM = ',
     &         I_CA_SM, I_CB_SM, I_AA_SM, I_AB_SM
             END IF
*. And internal strings
             DO I_IN_AB = 1, NST_IN_AB(I_IN_AB_SM,JJ_INTP)
             DO I_IN_AA = 1, NST_IN_AA(I_IN_AA_SM,JJ_INTP)
             DO I_IN_CB = 1, NST_IN_CB(I_IN_CB_SM,JJ_INTP)
             DO I_IN_CA = 1, NST_IN_CA(I_IN_CA_SM,JJ_INTP)
               IF(NTEST.GE.10000) THEN
                WRITE(6,*) ' I_IN_AB, I_IN_AA, I_IN_CB, I_IN_CA = ',
     &                       I_IN_AB, I_IN_AA, I_IN_CB, I_IN_CA
               END IF
              I_EI_OP = I_EI_OP + 1
*. Merge strings to obtain complete CA,CB,AA,AB strings
*. CA
              IOFF_IN_CA = 1+(I_IN_CA-1)*NEL_IN_CA(JJ_INTP)
              IOFF_EX_CA = 1+(I_EX_CA-1)*NEL_EX_CA
              CALL PROD_TWO_STRINGS(I_EI_CA,
     &        IST_EX_CA(IOFF_EX_CA,I_EX_CA_SM),
     &        IST_IN_CA(IOFF_IN_CA,I_IN_CA_SM,JJ_INTP),
     &        NEL_EX_CA,NEL_IN_CA(JJ_INTP),IS_CA)
*. CB
              IOFF_IN_CB = 1+(I_IN_CB-1)*NEL_IN_CB(JJ_INTP)
              IOFF_EX_CB = 1+(I_EX_CB-1)*NEL_EX_CB
              CALL PROD_TWO_STRINGS(I_EI_CB,
     &        IST_EX_CB(IOFF_EX_CB,I_EX_CB_SM),
     &        IST_IN_CB(IOFF_IN_CB,I_IN_CB_SM,JJ_INTP),
     &        NEL_EX_CB,NEL_IN_CB(JJ_INTP),IS_CB)
*. AA
              IOFF_IN_AA = 1+(I_IN_AA-1)*NEL_IN_AA(JJ_INTP)
              IOFF_EX_AA = 1+(I_EX_AA-1)*NEL_EX_AA
              CALL PROD_TWO_STRINGS(I_EI_AA,
     &        IST_EX_AA(IOFF_EX_AA,I_EX_AA_SM),
     &        IST_IN_AA(IOFF_IN_AA,I_IN_AA_SM,JJ_INTP),
     &        NEL_EX_AA,NEL_IN_AA(JJ_INTP),IS_AA)
*. AB
              IOFF_IN_AB = 1+(I_IN_AB-1)*NEL_IN_AB(JJ_INTP)
              IOFF_EX_AB = 1+(I_EX_AB-1)*NEL_EX_AB
              CALL PROD_TWO_STRINGS(I_EI_AB,
     &        IST_EX_AB(IOFF_EX_AB,I_EX_AB_SM),
     &        IST_IN_AB(IOFF_IN_AB,I_IN_AB_SM,JJ_INTP),
     &        NEL_EX_AB,NEL_IN_AB(JJ_INTP),IS_AB)
*. And Adresses of combined strings
C ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
*. CA
              I_CA_ADR =  ISTRNM(I_EI_CA,NTOOB,NEL_CA(JJ_INTP),
     &                    IZ_CA(1,JJ_INTP),IREO_CA(1,JJ_INTP),1)
*. CB
              I_CB_ADR =  ISTRNM(I_EI_CB,NTOOB,NEL_CB(JJ_INTP),
     &                    IZ_CB(1,JJ_INTP),IREO_CB(1,JJ_INTP),1)
*. AA
              I_AA_ADR =  ISTRNM(I_EI_AA,NTOOB,NEL_AA(JJ_INTP),
     &                    IZ_AA(1,JJ_INTP),IREO_AA(1,JJ_INTP),1)
*. AB
              I_AB_ADR =  ISTRNM(I_EI_AB,NTOOB,NEL_AB(JJ_INTP),
     &                    IZ_AB(1,JJ_INTP),IREO_AB(1,JJ_INTP),1)
*
              IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' CA, CB, AA, AB strings : '
               CALL IWRTMA(I_EI_CA,1,NEL_CA(JJ_INTP),1,NEL_CA(JJ_INTP))
               CALL IWRTMA(I_EI_CB,1,NEL_CB(JJ_INTP),1,NEL_CB(JJ_INTP))
               CALL IWRTMA(I_EI_AA,1,NEL_AA(JJ_INTP),1,NEL_AA(JJ_INTP))
               CALL IWRTMA(I_EI_AB,1,NEL_AB(JJ_INTP),1,NEL_AB(JJ_INTP))
               WRITE(6,'(A,4I5)') ' Adresses of CA, CB, AA, AB',
     &         I_CA_ADR, I_CB_ADR, I_AA_ADR, I_AB_ADR
               WRITE(6,*) ' Offset of symmetry-blocks in ST',
     &         IB_ST_SSS(I_CA_SM,I_CB_SM,I_AA_SM,JJ_INTP)
              END IF
* CA, CB, AA, AB adress
              IST_ADR = 
     &      (I_AB_ADR-1)*NST_AA(I_AA_SM,JJ_INTP)*NST_CB(I_CB_SM,JJ_INTP)
     &                  *NST_CA(I_CA_SM,JJ_INTP)
     &     +(I_AA_ADR-1)                        *NST_CB(I_CB_SM,JJ_INTP)
     &                  *NST_CA(I_CA_SM,JJ_INTP)
     &     +(I_CB_ADR-1)*NST_CA(I_CA_SM,JJ_INTP)
     &     + I_CA_ADR + IB_ST_TP(JJ_INTP) - 1 
     &     + IB_ST_SSS(I_CA_SM,I_CB_SM,I_AA_SM,JJ_INTP)-1
*
*. Signs
*. The sign is composed of two parts
* 1: Sign to bring ECA ECB EAA EAB ICA ICB IAA IAB into
*            ECA ICA ECB ICB EAA IAB EAB IAB
* 2. Sign to bring each of ECA ICA, ECB ICB, EAA IAB, EAB IAB
* into standard order
              IST_SIGN = IS_CA*IS_CB*IS_AA*IS_AB*ISIGNG
              IF(NTEST.GE.10000)  
     &        WRITE(6,*) ' I_EI_OP, IST_ADR, ISIGNG; IST_SIGN',
     &        I_EI_OP,  IST_ADR, ISIGNG, IST_SIGN
*. And now : what we all have been waiting for a few hundred lines..
              IREO_EI_ST(I_EI_OP) = IST_SIGN*IST_ADR
*
             END DO
             END DO
             END DO
             END DO
*.           ^ End of loop over internal CA,CB,AA,AB strings of given sym
            END DO
            END DO
            END DO
*.          ^ End of loop over symmetry of internal CA, CB, AA, AB strings
           END DO
*.         ^ End of loop over internal types
          END DO
          END DO
          END DO
          END DO
*         ^ End of loop over external CA, CB,AA,AB strings of given sym
        END DO
        END DO
        END DO
*.      ^ End of loop over symmetry of external CA, CB, AA, AB strings
       END DO
*.     ^ End of loop over symmetry of external operators
      END DO
*.    ^ End of loop over external types
      N_EI_OP = I_EI_OP
*
*. Sum check
      I_DO_SUMCHECK = 1
      IF(I_DO_SUMCHECK.EQ.1.AND.N_EI_OP.LT.30000) THEN
*. Check that sum SUM_I IREO_EI_ST(I) = N_EI_OP*(N_EI_OP+1)/2
        ISUM = 0
        DO I = 1, N_EI_OP
          ISUM = ISUM + ABS(IREO_EI_ST(I))
        END DO
        IF(ISUM.EQ.N_EI_OP*(N_EI_OP+1)/2) THEN
          WRITE(6,*) ' Sumcheck passed '
        ELSE
          WRITE(6,*) ' Sumcheck failed in REO_EI..'
          WRITE(6,'(A,2I7)') 
     &    ' Observed and expected sum', ISUM,N_EI_OP*(N_EI_OP+1)/2
*. Determine elements that were not obtained exactly one time- 
*, N-squared algorithm in use, as I do not want to allocate another
*. array..
          DO ITARGET = 1, N_EI_OP
           NFOUND = 0
           DO I = 1,N_EI_OP
             IF(IREO_EI_ST(I).EQ.ITARGET) NFOUND = NFOUND + 1
           END DO
           IF(NFOUND.NE.1) THEN
             WRITE(6,*) 
     &       ' Element ', ITARGET ,' Found ', NFOUND ,' times'
           END IF
          END DO
*. Print reorder array before stop
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' The reorder array EI => standard order '
            CALL IWRTMA(IREO_EI_ST,1,N_EI_OP,1,N_EI_OP,1)
          END IF
C         STOP ' sumcheck failed '
        END IF
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The reorder array EI => standard order '
        CALL IWRTMA(IREO_EI_ST,1,N_EI_OP,1,N_EI_OP,1)
      END IF
*
      RETURN
      END 
      SUBROUTINE PROD_TWO_STRINGS(ISTR_AB,ISTR_A,ISTR_B,NEL_A,NEL_B,IS)
*
* Two strings are given in standard ascending order
* Obtain product of string in ascending order, and sign 
* required for conversion
*
*. Jeppe Olsen, October 2006
*. Input
      INTEGER ISTR_A(NEL_A),ISTR_B(NEL_B)
*. Output
      INTEGER ISTR_AB(NEL_A+NEL_B)
*
      IEL_A = 1
      IEL_B = 1
      NPERM = 0
      DO IEL_AB = 1, NEL_A+NEL_B
        IF(IEL_B.GT.NEL_B) THEN
*. No more B-electrons so next electron is from A
           ISTR_AB(IEL_AB) = ISTR_A(IEL_A)
           IEL_A = IEL_A + 1
        ELSE IF (IEL_A.GT.NEL_A) THEN
*. No more A-electrons so next electron is from B
           ISTR_AB(IEL_AB) = ISTR_B(IEL_B)
           IEL_B = IEL_B + 1
        ELSE
*. There are both A and B-electrons left so compare
          IF(ISTR_A(IEL_A).LE.ISTR_B(IEL_B)) THEN
*. Next electron is from IEL_A
            ISTR_AB(IEL_AB) = ISTR_A(IEL_A)
            IEL_A = IEL_A + 1
          ELSE 
*. Next electron is from IEL_B
            ISTR_AB(IEL_AB) = ISTR_B(IEL_B)
            IEL_B = IEL_B + 1
*. Number of permutations over A-electrons to get it in place
            NPERM = NPERM + NEL_A-IEL_A+1
          END IF
        END IF
      END DO
*. Well, I am not sure that the above count of permutations is correct
*. so here is another try-starting with the smallest b, and checking how
*. many a's it should be moved passed
      NPERM = 0
      DO IEL_B = 1, NEL_B
        IORB_B = ISTR_B(IEL_B)
*. Find largest index in A that is smaller than IORB_B
        I_SMALL = 0
        DO IEL_A = 1, NEL_A
         IF(ISTR_A(IEL_A).LT.IORB_B) I_SMALL = IEL_A
        END DO
        NPERM = NPERM + NEL_A-I_SMALL
      END DO
*
      IF(MOD(NPERM,2).EQ.0) THEN
        IS = 1
      ELSE
        IS = -1
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        NEL_AB = NEL_A + NEL_B
        WRITE(6,*) 'Merging two strings : AB, A, B '
        WRITE(6,*) 'Number of operators in A,B', NEL_A, NEL_B
        CALL IWRTMA(ISTR_AB,1,NEL_AB,1,NEL_AB)
        CALL IWRTMA(ISTR_A,1,NEL_A,1,NEL_A)
        CALL IWRTMA(ISTR_B,1,NEL_B,1,NEL_B)
      END IF
*
      RETURN
      END 
      SUBROUTINE EITP_TO_SPOXTP_REO(I_EI_TP_REO,ISPOBEX_TP,NSPOBEX_TP,
     &           IB_EXTP, IB_INTP,NGAS,N_EXTP,N_INTP,
     &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX)
*
* Obtain reorder array of spinorbital excitation types from EI
* to standard order
*
*. Jeppe Olsen
*
*. October 2006
*
      INCLUDE 'wrkspc.inc'
*. Input : CAAB for all spinorbitalexcitations
      INTEGER ISPOBEX_TP(4,NGAS,*)
      INTEGER N_IN_FOR_EX(N_EXTP),IB_IN_FOR_EX(N_EXTP), I_IN_FOR_EX(*)
*. Output
      INTEGER I_EI_TP_REO(NSPOBEX_TP)
*. Local scratch
      INTEGER I_EI_MERGE(4*MXPNGAS)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from EITP_TO_SPOXTP_REO '
        WRITE(6,*) ' ============================='
      END IF
      IZERO = 0
      CALL ISETVC(I_EI_MERGE,IZERO,4*NGAS)
      I_EI_TP = 0
*
C?    WRITE(6,*) ' NSPOBEX_TP, N_EXTP = ', NSPOBEX_TP, N_EXTP
      DO J_EXTP = 1, N_EXTP
C?      WRITE(6,*) ' J_EXTP = ', J_EXTP
        L_INTP = N_IN_FOR_EX(J_EXTP)
        IB_JEX  = IB_IN_FOR_EX(J_EXTP)
        DO JJ_INTP = 1, L_INTP
          I_EI_TP = I_EI_TP + 1
          J_INTP = I_IN_FOR_EX(IB_JEX-1+JJ_INTP)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' J_EXTP, J_INTP, I_EI_TP = ',
     &                   J_EXTP, J_INTP, I_EI_TP
          END IF
          I_IN_ADR = IB_INTP - 1 + J_INTP
          I_EX_ADR = IB_EXTP - 1 + J_EXTP
C IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
*. Merge Internal and external types
          IONE = 1
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' J_EXTP, J_INTP: '
            CALL WRT_SPOX_TP(ISPOBEX_TP(1,1,I_EX_ADR),1)
            CALL WRT_SPOX_TP(ISPOBEX_TP(1,1,I_IN_ADR),1)
          END IF
*
          CALL IVCSUM(I_EI_MERGE,
     &         ISPOBEX_TP(1,1,I_IN_ADR),ISPOBEX_TP(1,1,I_EX_ADR),
     &         1,1,4*NGAS)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' I_EI_MERGE: '
            CALL  WRT_SPOX_TP(I_EI_MERGE,1)
          END IF
*. Find this type in SPOBEX
          IFOUND = 0
          DO JSPOBEX_TP = 1, NSPOBEX_TP
C                COMPARE_TWO_INTARRAYS(IA,IB,NAB,IDENT)
            CALL COMPARE_TWO_INTARRAYS(
     &           I_EI_MERGE,ISPOBEX_TP(1,1,JSPOBEX_TP),4*NGAS,IDENT)
            IF(IDENT.NE.0) IFOUND = JSPOBEX_TP
          END DO
          I_EI_TP_REO(I_EI_TP) = IFOUND
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Reorder array for types: EI => standard order '
        CALL IWRTMA(I_EI_TP_REO,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      RETURN
      END
C     CALL GET_INTERNAL_STATES(N_EXTOP_TP,N_INTOP_TP,
C    &     WORK(KLSOBEX),WORK(KL_N_INT_FOR_EXT),WORK(KL_IB_INT_FOR_EXT),
C    &     WORK(KL_I_INT_FOR_EXT),WORK(KL_NDIM_IN_SE),
C    &     WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
C    &     WORK(KL_X1_INT_EI_FOR_SE), WORK(KL_X2_INT_EI_FOR_SE),
C    &     WORK(KL_SG_INT_EI_FOR_SE)
C    &     WORK(KL_IBX1_INT_EI_FOR_SE), WORK(KL_IBX2_INT_EI_FOR_SE),
C    &     WORK(KL_IBSG_INT_EI_FOR_SE)
C    &     WORK(KL_X2L_INT_EI_FOR_SE),
C    &     I_IN_TP,I_INT_OFF,I_EXT_OFF)
      SUBROUTINE GET_INTERNAL_STATES(N_EXTP,N_INTP,ISPOBEX,
     &           N_IN_FOR_EX,IB_IN_FOR_EX,I_IN_FOR_EX,
     &           N_IEL_FOR_SE,N_ORTN_FOR_SE,N_INT_FOR_SE,
     &           X1_INT_EI_FOR_SE,X2_INT_EI_FOR_SE,SG_INT_EI_FOR_SE,
     &           S_INT_EI_FOR_SE,
     &           IBX1_INT_EI_FOR_SE,IBX2_INT_EI_FOR_SE,
     &           IBSG_INT_EI_FOR_SE, IBS_INT_EI_FOR_SE,
     &           X2L_INT_EI_FOR_SE,
     &           I_INT_TP,IB_INTP,IB_EXTP)
*
*. Obtain the orthonormal internal states for each EI class.
*. If I_INT_TP = 1, then the internal states are the 
*  the nonorthogonal states of the metric
*. IF I_INT_TP = 2, then the internal states are 
*. obtained by diagonalizing the internal Hamiltonian in the 
*. space of orthonormal internal states
*. If I_INT_TP = 3, then the Jacobian in the internal space 
*. is constructed and diagonalized and the left and right hand 
*. side eigenvectors are obtained.
*. The transformation matrix between actual orthonormal and elementary internal states
*
* Currently a CASSCF state is assumed
*
*. October 28, 2007 in Crete, being at conference at Stines 28 years birthday..
*
*. Continued March 17, 2008, on the train back from Warsaw 
*
* Finalized(!) in Pisa Oct 6 (2008 (no comments, please...))
*. Nope, modified march 12, 2009- and silence is still appreciated
*
*. Metric over internal states added Oct. 2009
*
*. Improving, Oct. 2012
*
* Last Modification; Oct. 18, 2012; Jeppe Olsen, Madrid; Batching to disc, etc
*.(Debugged in the air from Madrid to Amsterdam....)
*.
*
*. Form of internal Hamiltonian is defined by  I_INT_HAM from CRUN
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc' 
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'crun.inc'
*. Input
      INTEGER N_IN_FOR_EX(N_EXTP),IB_IN_FOR_EX(N_EXTP),I_IN_FOR_EX(*)
*. Number of elementary internal operators of given sym
*. for given external type
      INTEGER N_IEL_FOR_SE(NSMOB,N_EXTP)
      INTEGER ISPOBEX(NGAS,4,*)
*. Output
*. Number of orthonormal internal states per symmetry and external type
      INTEGER N_ORTN_FOR_SE(NSMOB,N_EXTP)
*. Number of included internal states per symmetry and external types
      INTEGER N_INT_FOR_SE(NSMOB,N_EXTP)
*. Transformation matrix X1(EI,INT) Diagonalizing metric
*. various internal states
      DIMENSION X1_INT_EI_FOR_SE(*)
*.    Nonvanishing eigenvalues sigma of metric
      DIMENSION SG_INT_EI_FOR_SE(*)
*. Overlap of internal states
      DIMENSION S_INT_EI_FOR_SE(*)
*. Tranformation matrix diagonalizing zero order Hamiltonian in
*  orthonormal basis
      DIMENSION X2_INT_EI_FOR_SE(*)
*. lefthand side transformation for X2 if I_IN_TP = 3
      DIMENSION X2L_INT_EI_FOR_SE(*)
*. and offsets for X1
      INTEGER IBX1_INT_EI_FOR_SE(NSMOB,N_EXTP)
*. and offsets for X2
      INTEGER IBX2_INT_EI_FOR_SE(NSMOB,N_EXTP)
*. and offsets for Sigma
      INTEGER IBSG_INT_EI_FOR_SE(NSMOB,N_EXTP)
*. and offsets for S
      INTEGER IBS_INT_EI_FOR_SE(NSMOB,N_EXTP)
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ---------------------------- '
        WRITE(6,*) ' GET_INTERNAL_STATES Speaking '
        WRITE(6,*) ' ---------------------------- '
        WRITE(6,*) 
        WRITE(6,*) ' I_INT_TP, I_INT_HAM = ', I_INT_TP,I_INT_HAM
      END IF
*
      I12_SAVE = I12
      PSSIGN_SAVE = PSSIGN
      IDC_SAVE = IDC
*
      PSSIGN = 0.0D0
      IDC = 1
*. Incore or out-of core construction of matrices
      I_IN_OR_OUT = 2
      IF(I_IN_OR_OUT.EQ.1) THEN
        WRITE(6,*) ' In core construction of matrices'
        LUSCR_INT = -1
        LUSCR_INT2 = -1
      ELSE
        WRITE(6,*) ' Out of core construction of matrices'
C            FILEMAN_MINI(IFILE,ITASK)
        CALL FILEMAN_MINI(LUSCR_INT,'ASSIGN')
        CALL FILEMAN_MINI(LUSCR_INT2,'ASSIGN')
        CALL REWINO(LUSCR_INT)
        CALL REWINO(LUSCR_INT2)
      END IF 
*
      IB_X1 = 1
      IB_X2 = 1
      IB_SG = 1
      IB_S = 1
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',1,'GT_INS')
      CALL LUCIAQENTER('GT_INS')
*. Scratch space  : Two matrices S(INT1,INT2), where INT1, INT2 runs over 
*                   internal states with given symmetry 
*                   belonging to a given type of external state
*
*. Obtain occupations of alpha- and beta-strings of reference space - 
*. is currently assumed to be a single space
C GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BEC
      CALL MEMMAN(KL_REF_AL,NGAS,'ADDL  ',2,'REF_AL')
      CALL MEMMAN(KL_REF_BE,NGAS,'ADDL  ',2,'REF_BE')
      CALL MEMMAN(KL_IREF_AL,NGAS,'ADDL  ',2,'IEF_AL')
      CALL MEMMAN(KL_IREF_BE,NGAS,'ADDL  ',2,'IEF_BE')
*. NOTE : reference space is assumed to be space 1
      IREFSPC = 1
      CALL GET_REF_ALBE_OCC(IREFSPC,WORK(KL_REF_AL),WORK(KL_REF_BE))
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Occupation of alpha- and beta-strings in reference'
        CALL IWRTMA(WORK(KL_REF_AL),1,NGAS,1,NGAS)
        CALL IWRTMA(WORK(KL_REF_BE),1,NGAS,1,NGAS)
      END IF
*Two matrices over internal states for given external type
      NDIM_II = IMXMN(1,N_IEL_FOR_SE,NSMOB*N_EXTP)
      IF(NTEST.GE.10) WRITE(6,*) 
     &' Largest number of internal state for given sym and E-type', 
     & NDIM_II
      CALL MEMMAN(KL_MATII1,NDIM_II**2,'ADDL  ',2,'MATII1')
      CALL MEMMAN(KL_MATII2,NDIM_II**2,'ADDL  ',2,'MATII2')
*. and two vectors over internal states
      CALL MEMMAN(KL_VECII1,NDIM_II,'ADDL  ',2,'VECII1')
      CALL MEMMAN(KL_VECII2,NDIM_II,'ADDL  ',2,'VECII2')
      
*. In the following we are going to use/misuse the standard
*. spinorbital types stored as the first elements in the spinorbitalexcitation arrays.
*. Save the information usually stored there
*. Arrays
      CALL MEMMAN(KLSOBEX_SAVE,(NSPOBEX_TP+1)*NGAS*4,'ADDL ',1,'SPOXSV')
      CALL ICOPVE(ISPOBEX,WORK(KLSOBEX_SAVE),(NSPOBEX_TP+1)*NGAS*4)
      NSPOBEX_TP_SAVE = NSPOBEX_TP
*. Loop over the various EI spaces
      DO I_EXTP = 1, N_EXTP
       IF(NTEST.GE.10) THEN
         WRITE(6,*) 'Output for external type = ', I_EXTP
         WRITE(6,*) '==============================='
       END IF
       L_INTP = N_IN_FOR_EX(I_EXTP)
       IOFF = IB_IN_FOR_EX(I_EXTP)
C?     WRITE(6,*) ' L_INTP, IOFF = ', L_INTP, IOFF
*. Prepare the arrays defining calculations for this space
       NSPOBEX_TP = L_INTP
       DO II_INTP = 1, L_INTP
          I_INTP = I_IN_FOR_EX(IOFF-1+II_INTP)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' II_INTP, I_INTP ', II_INTP, I_INTP
            WRITE(6,*) ' IB_INTP,NGAS =', IB_INTP,NGAS
          END IF
          CALL ICOPVE(ISPOBEX(1,1,IB_INTP-1+I_INTP),
     &                ISPOBEX(1,1,II_INTP),4*NGAS)
       END DO
*
*----------------------------------
*. Obtain metric in internal space
* ---------------------------------
*
*
* Type of C-coefficients will be type of reference
       ICATP = 1
       ICBTP = 2
*
       N_ALEL_C = NELFTP(ICATP)
       N_BEEL_C = NELFTP(ICBTP)
*. Find how action of internal operator changes occupation, 
*. as all operators makes identical changes of occupations, 
*. it is sufficient to look on the first operator
COLD   I_INTP1 = I_IN_FOR_EX(IB_INTP-1+1)
       IALDEL = IELSUM(ISPOBEX(1,1,1),NGAS)-IELSUM(ISPOBEX(1,3,1),NGAS)
       IBEDEL = IELSUM(ISPOBEX(1,2,1),NGAS)-IELSUM(ISPOBEX(1,4,1),NGAS)
C?     WRITE(6,*) ' IALDEL, IBEDEL = ', IALDEL, IBEDEL
*. Occupation of internal operator times reference space
C     CAAB_T_ABSTR(ICAAB,IAOC_IN,IBOC_IN,IAOC_OUT,IBOC_OUT,NGAS )
       CALL CAAB_T_ABSTR(ISPOBEX(1,1,1),WORK(KL_REF_AL),WORK(KL_REF_BE),
     &                   WORK(KL_IREF_AL),WORK(KL_IREF_BE),NGAS)
*. We now have occupation of resulting strings, find string numbers
C FIND_SPGRP_FROM_OCC(IOCC,ISPGRP_NUM,ITP)
*. Types of alpha and beta strings 
       N_ALEL_S = N_ALEL_C + IALDEL
       N_BEEL_S = N_BEEL_C + IBEDEL
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Number of alpha-electrons in IOP*|ref> ',N_ALEL_S
         WRITE(6,*) ' Number of beta- electrons in IOP*|ref> ',N_BEEL_S
       END IF
*. Find supergroup types with these number of electrons
       ISATP = 0
       ISBTP = 0
       DO ISPGP_TP = 1, NSTTP
         IF(NELFTP(ISPGP_TP).EQ.N_ALEL_S) ISATP = ISPGP_TP
         IF(NELFTP(ISPGP_TP).EQ.N_BEEL_S) ISBTP = ISPGP_TP
       END DO
       IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' a-,b-Types of internal operator times ref',
     &             ISATP,ISBTP
       END IF
*. Below is not neccesary
       CALL FIND_SPGRP_FROM_OCC(WORK(KL_IREF_AL),IAL_SPGP_S,ISATP)
       CALL FIND_SPGRP_FROM_OCC(WORK(KL_IREF_BE),IBE_SPGP_S,ISBTP)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' a-b-supergroups of internal operator times ref',
     &   ISATP,ISBTP
       END IF
*
       IF(ISATP.EQ.0.OR.ISBTP.EQ.0) THEN
         WRITE(6,*) 
     &   ' String type with modified electroncount does not exist'
         WRITE(6,*) ' NAEL, NBEL, ISATP, ISBTP = ',
     &                N_ALEL_S, N_BEEL_S, ISATP, ISBTP
         STOP       
     &   ' String type with modified electroncount does not exist'
       END IF
*
*. Occupations of internal operator times reference space
*
* Generate a space with this occupation. This space is 
* stored as the last CI space and combination space. If the numbers 
* of these equals MXPICI I am trouble 
*
       IF(NCISPC.EQ.MXPICI.OR.NCMBSPC.EQ.MXPICI) THEN
         WRITE(6,*) ' Not space for an extra CI space '
         WRITE(6,*) ' Increase parameter MXPICI'
         WRITE(6,*) ' NCISPC, MXPICI ', NCISPC, MXPICI
         WRITE(6,*) ' NCMBSPC, MXPICI ', NCMBSPC, MXPICI
         STOP ' Increase parameter MXPICI'
       END IF
       IMODREF_CISPC = NCISPC + 1
       IMODREF_CMBSPC = NCMBSPC + 1
       LCMBSPC(IMODREF_CMBSPC) = 1
       ICMBSPC(1,IMODREF_CMBSPC) = IMODREF_CISPC
*
       IALTP_FOR_GAS(IMODREF_CMBSPC) = ISATP
       IBETP_FOR_GAS(IMODREF_CMBSPC) = ISBTP
*. Occupation constraints for modified reference : electrons are 
*. added or removed from the active space. At the moment a single 
* active space is assumed.
*. Copy reference space occupation 
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Min and max for reference space '
         CALL IWRTMA(IGSOCCX(1,1,IREFSPC),1,NGAS,1,NGAS)
         CALL IWRTMA(IGSOCCX(1,2,IREFSPC),1,NGAS,1,NGAS)
       END IF
*
       CALL ICOPVE(IGSOCCX(1,1,IREFSPC),IGSOCCX(1,1,IMODREF_CISPC),
     &             NGAS)
       CALL ICOPVE(IGSOCCX(1,2,IREFSPC),IGSOCCX(1,2,IMODREF_CISPC),
     &             NGAS)
*. The active space
       IACT_GAS = 0
       DO IGAS = 1, NGAS
         IF(IHPVGAS(IGAS).EQ.3) IACT_GAS = IGAS
       END DO
       IF(NTEST.GE.1000) WRITE(6,*) ' The active GASpace ', IACT_GAS
*. Change number of electrons in active space
       DO IGAS = 1, NGAS
         IF(IGAS.GE.IACT_GAS) THEN
           IGSOCCX(IGAS,1,IMODREF_CISPC) =
     &     IGSOCCX(IGAS,1,IMODREF_CISPC) + IALDEL + IBEDEL
           IGSOCCX(IGAS,2,IMODREF_CISPC) =
     &     IGSOCCX(IGAS,2,IMODREF_CISPC) + IALDEL + IBEDEL
         END IF
       END DO
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Min and max for modified reference space '
         CALL IWRTMA(IGSOCCX(1,1,IMODREF_CISPC),1,NGAS,1,NGAS)
         CALL IWRTMA(IGSOCCX(1,2,IMODREF_CISPC),1,NGAS,1,NGAS)
       END IF
*. Loop over symmetry of internal operator
       DO IOPSM = 1, NSMOB
        IF(NTEST.GE.10) WRITE(6,'(A,2I5)') 
     &  ' Info for I_EXTP, IOPSM = ', I_EXTP, IOPSM
        IBX1_INT_EI_FOR_SE(IOPSM,I_EXTP) = IB_X1
        IBX2_INT_EI_FOR_SE(IOPSM,I_EXTP) = IB_X2
        IBSG_INT_EI_FOR_SE(IOPSM,I_EXTP) = IB_SG
        IBS_INT_EI_FOR_SE(IOPSM,I_EXTP)  = IB_S
*.(IB_X1, IB_X2, IB_SG will be updated at end of loop over IOPSM)
*. Symmetry of operator times reference state
        IOPREFSM = MULTD2H(IOPSM,IREFSM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOPSM, IREFSM, IOPREFSM = ', 
     &                 IOPSM, IREFSM, IOPREFSM
        END IF
*. Number of determinants in modified internal space with symmetry ISYM
        LDET = LEN_CISPC(IMODREF_CMBSPC,IOPREFSM,NTEST)
*. Number of internal operators of this sym
        LINT = N_IEL_FOR_SE(IOPSM,I_EXTP)
*. Copy to N_INT_FOR_SE
        N_INT_FOR_SE(IOPSM,I_EXTP) = LINT
        
        IF(NTEST.GE.10)
     &  WRITE(6,*) ' IOPSM, LDET, LINT = ', IOPSM,LDET,LINT
        IF(LDET.EQ.0.OR.LINT.EQ.0) THEN
*. No novanishing states of this symmetry, so
          N_ORTN_FOR_SE(IOPSM,I_EXTP) = 0
        ELSE
*. Not trivial zero, so look further ( for some hundred of lines..)
*. Space for expansion Int |ref> in SD for all Internal states
        IDUM = 0
        CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'INTMAT')
        IF(I_IN_OR_OUT.EQ.1) THEN
          LDIM = LINT*LDET
        ELSE
          LDIM = MAX(LINT,LDET)
        END IF
        CALL MEMMAN(KL_INTMAT,LDIM,'ADDL  ',2,'INTMAT')
        CALL MEMMAN(KL_INTMAT2,LDIM,'ADDL  ',2,'INTMT2')
        CALL MEMMAN(KL_INTOP,LINT**2,'ADDL  ',2,'INTOP ')
        CALL MEMMAN(KL_INTOP2,LINT**2,'ADDL  ',2,'INTOP2')
        CALL MEMMAN(KL_INTOPV,LINT,'ADDL  ',2,'INTOP ')
*
        LWORK = 4*LINT
        CALL MEMMAN(KL_WORK,LWORK,'ADDL  ',2,'INTOP ')
        CALL MEMMAN(KL_IWORK,LWORK,'ADDL  ',2,'INTOP ')
*
        ICSPC = IREFSPC
        ISSPC = IMODREF_CISPC
        ICSM = IREFSM
        ISSM = IOPREFSM
        IF(I_IN_OR_OUT.EQ.2) CALL REWINO(LUSCR_INT)
        DO INTOP = 1, LINT
          ZERO = 0.0D0
          CALL SETVEC(WORK(KL_INTOP),ZERO,LINT)
          WORK(KL_INTOP-1+INTOP) = 1.0D0
          ICSPC = IREFSPC
          ISSPC = IMODREF_CISPC
          ICSM = IREFSM
          ISSM = IOPREFSM
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &         WORK(KL_INTOP),1)
          CALL REWINO(LUHC)
          IF(NTEST.GE.10000) THEN
            WRITE(6,*) ' Result of SIGDEN_CC on LUHC '
            CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
          END IF
          CALL REWINO(LUHC)
          IF(I_IN_OR_OUT.EQ.1) THEN
            CALL FRMDSCN(WORK(KL_INTMAT+(INTOP-1)*LDET),-1,-1,LUHC)
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' INTOP * |ref> for INTOP = ', INTOP
              CALL WRTMAT(WORK(KL_INTMAT+(INTOP-1)*LDET),1,LDET,1,LDET)
            END IF
          ELSE
            CALL COPVCD(LUHC,LUSCR_INT,WORK(KL_INTMAT),0,-1)
          END IF
C?        STOP ' Enforced stop after first call to SIGDEN_CC'
        END DO ! Loop over internal states
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' The matrix X(IDET,IELOP) '
          IF(I_IN_OR_OUT.EQ.1) THEN
            CALL WRTMAT(WORK(KL_INTMAT),LDET,LINT,LDET,LINT)
          ELSE
            CALL REWINO(LUSCR_INT)
            DO INTOP = 1, LINT
              CALL WRTVCD(WORK(KL_INTMAT),LUSCR_INT,0,-1)
            END DO
          END IF
        END IF !NTEST is sufficiently large
*. Set up the overlap matrix S(I,J) = <ref!op(i)+ op(j)!ref>
        IF(I_IN_OR_OUT.EQ.1) THEN
*. In house construction 
          IJ = 0
          DO I = 1, LINT
            DO J = 1, I
              IJ = IJ + 1
              WORK(KL_MATII1-1+I*(I-1)/2+J) =
     &        INPROD(WORK(KL_INTMAT+(I-1)*LDET),
     &               WORK(KL_INTMAT+(J-1)*LDET),LDET)
            END DO
          END DO
        ELSE 
*. Disc based construction
*. Determine amount of memory that can be used for storing expansion
*. of external states
C MEMMAN(KBASE,KADD,TASK,IR,IDENT)
          KBASE = 0
          CALL MEMMAN(KBASE,KFREE,'FREE  ',IDUM,'CHKFRE')
          KMAX = 100000000
          KDISC = MIN(KMAX,KFREE,LDET*LINT)
          WRITE(6,*) 
     &    ' Amount of memory to be used for batches of int. states ',
     &    KDISC
*. Number of states per batch
          NSTA_BAT = KDISC/LDET
          IF(NSTA_BAT.EQ.0) THEN
            WRITE(6,*) 
     &      ' Problem: Batch cannot contain a single internal state'
            WRITE(6,*) ' KDISC, LDET = ', KDISC, LDET
            STOP
     &      ' Problem: Batch cannot contain a single internal state'
          END IF
          NBAT = LINT/NSTA_BAT
          IF(NBAT*NSTA_BAT.LT.LINT) NBAT = NBAT + 1
          CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'BATINT')
          CALL MEMMAN(KLCBAT,KDISC,'ADDL  ',2,'CBAT  ')
          DO JBAT = 1, NBAT
            J_INI = (JBAT-1)*NSTA_BAT + 1
            J_END = MIN(JBAT*NSTA_BAT, LINT)
            LBAT = J_END + 1 - J_INI
*. Read J-states in
            CALL SKPVCD(LUSCR_INT,J_INI-1,WORK(KL_INTMAT),1,-1)
            DO JSTA = 1, LBAT
             CALL FRMDSCN(WORK(KLCBAT+(JSTA-1)*LDET),-1,-1,LUSCR_INT)
            END DO
*. Loop over I-states and generate overlap
            CALL REWINO(LUSCR_INT)
            DO I = 1, LINT
             CALL FRMDSCN(WORK(KL_INTMAT),-1,-1,LUSCR_INT)
             DO J = J_INI, I
               IJ = I*(I-1)/2 + J
               WORK(KL_MATII1-1+I*(I-1)/2+J) =
     &              INPROD(WORK(KL_INTMAT),
     &              WORK(KLCBAT+(J-J_INI)*LDET),LDET)
             END DO
            END DO ! End of loop over pair of states
          END DO ! End of loop over batches of states
          CALL MEMMAN(IDUM,IDUN,'FLUSM ',IDUM,'BATINT')
        END IF ! in house or disc version
*
        IF(NTEST.GE.100) THEN
         WRITE(6,*) ' The overlap  <ref!op(i)+ op(j)!ref>'
         CALL PRSYM(WORK(KL_MATII1),LINT)
        END IF
        CALL
     &  COPVEC(WORK(KL_MATII1),S_INT_EI_FOR_SE(IB_S),LINT*(LINT+1)/2)
*
* -------------------
*. Diagonalize metric
* -------------------
*
*. Expand to complete form
C TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
        CALL TRIPAK(WORK(KL_MATII2),WORK(KL_MATII1),2,LINT,LINT)
*. Obtain nonsingular orthogonal internal states
C E CHK_S_FOR_SING(S,NDIM,NSING,X,SCR,SCR2)
        CALL CHK_S_FOR_SING2(WORK(KL_MATII2),LINT,NSING,
     &       WORK(KL_MATII2),WORK(KL_VECII1),WORK(KL_VECII2),
     &       THRES_SINGU)
        NNONSING = LINT - NSING
        IF(NTEST.GE.10)
     &  WRITE(6,*) ' Number of singular and nonsingular states ', 
     &  NSING,NNONSING
        N_ORTN_FOR_SE(IOPSM,I_EXTP) = NNONSING
*. The nonsingular basis are the last NNONSING vectors in WORK(KLMATII2)
*. copy these to first X1
        DO IORTN = 1, NNONSING
          CALL COPVEC(WORK(KL_MATII2+(NSING+IORTN-1)*LINT),
     &                X1_INT_EI_FOR_SE(IB_X1+(IORTN-1)*LINT),LINT)
          SG_INT_EI_FOR_SE(IB_SG-1+IORTN) =
     &    WORK(KL_VECII1-1+NSING+IORTN)
        END DO
*
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' internal states diagonalizing metric'
         CALL WRTMAT(X1_INT_EI_FOR_SE(IB_X1),
     &               LINT,NNONSING,LINT,NNONSING)
        END IF
        IF(NTEST.GE.20) THEN
         WRITE(6,*) ' Eigenvalues of metric'
         CALL WRTMAT(SG_INT_EI_FOR_SE(IB_SG),1,NNONSING,1,NNONSING)
        END IF
        IF(NTEST.GE.20) THEN
          WRITE(6,*) 
     &    ' Info for construction of matrices over orth. states'
        END IF
*
        IF(I_INT_TP.GE.2) THEN
*. Obtain - if required - zero-order Hamiltonian in internal 
*. space and diagonalize this
C. H0 Intop |ref> for all orthonormal states
          IF(I_IN_OR_OUT.EQ.2) CALL REWINO(LUSCR_INT2)
          IF(I_IN_OR_OUT.EQ.2) CALL REWINO(LUSCR_INT)
          DO IORTN = 1, NNONSING
            IF(NTEST.GE.1000) WRITE(6,*) ' IORTN = ', IORTN
*. Scale to obtain orthonormal state
            FACTOR = 1.0D0/SQRT(SG_INT_EI_FOR_SE(IB_SG-1+IORTN))
            CALL COPVEC(X1_INT_EI_FOR_SE(IB_X1+(IORTN-1)*LINT),
     &      WORK(KL_VECII1),LINT)
            CALL SCALVE(WORK(KL_VECII1),FACTOR,LINT)
*
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
            ICSM = IREFSM
            ISSM = IOPREFSM
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUHC)
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KL_VECII1),1)
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' After sigden '
              WRITE(6,*) ' Resulting vector : '
              CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
            END IF
            IF(I_IN_OR_OUT.EQ.1) THEN
             CALL REWINO(LUHC)
             CALL FRMDSCN(WORK(KL_INTMAT+(IORTN-1)*LDET),-1,-1,LUHC)
             IF(NTEST.GE.10000) THEN
             WRITE(6,*) ' Intop (IORTN) |ref>'
              CALL WRTMAT(WORK(KL_INTMAT+(IORTN-1)*LDET),1,LDET,1,LDET)
             END IF
            ELSE
              CALL REWINO(LUHC)
              CALL COPVCD(LUHC,LUSCR_INT,WORK(KL_INTMAT),0,-1)
            END IF
*. zero-order Hamiltonian times Int(iortn)|ref>
            ICSPC = IMODREF_CMBSPC
            ISSPC = IMODREF_CMBSPC
            ICSM = IOPREFSM
            ISSM = IOPREFSM
            IF(I_INT_HAM.EQ.1) THEN
C?           WRITE(6,*) 
C?   &       ' One-body H0 used for internal zero-order states'
             I12 = 1
             CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
            ELSE
C?           WRITE(6,*) 
C?   &       ' Two-body H used for internal zero-order states'
            END IF
*
            CALL REWINO(LUHC)
            CALL REWINO(LUSC51)
            CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC51,0,0)
            CALL MEMCHK2('AFTMV7')
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' IMODREF_CMBSPC =', IMODREF_CMBSPC
              WRITE(6,*) ' After mv7 '
              WRITE(6,*) ' Resulting vector : '
              CALL WRTVCD(WORK(KVEC1P),LUSC51,1,-1)
            END IF
            IF(I_INT_HAM.EQ.1) 
     &      CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
            CALL REWINO(LUSC51)
            IF(I_IN_OR_OUT.EQ.1) THEN
              CALL FRMDSCN(WORK(KL_INTMAT2+(IORTN-1)*LDET),-1,-1,LUSC51)
            ELSE
              CALL COPVCD(LUSC51,LUSCR_INT2,WORK(KL_INTMAT),0,-1)
            END IF
          END DO
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Oint!0> and H0 Oint!0> expansions:'
            IF(I_IN_OR_OUT.EQ.1) THEN
               CALL WRTMAT(WORK(KL_INTMAT),LDET,NNONSING,
     &                     LDET,NNONSING)
               WRITE(6,*)
               CALL WRTMAT(WORK(KL_INTMAT2),LDET,NNONSING,
     &                     LDET,NNONSING)
               WRITE(6,*)
            ELSE
              CALL REWINO(LUSCR_INT)
              DO I = 1, NNONSING 
                CALL WRTVCD(WORK(KVEC1P),LUSCR_INT,0,-1)
              END DO
              WRITE(6,*)
              CALL REWINO(LUSCR_INT2)
              DO I = 1, NNONSING 
                CALL WRTVCD(WORK(KVEC1P),LUSCR_INT2,0,-1)
              END DO
            END IF
          END IF

*. <ref|intop+ H+ intop|ref>: for test complete matrix and metric are 
*. calculated
          IF(I_IN_OR_OUT.EQ.1) THEN
           DO IORTN = 1, NNONSING
            DO JORTN = 1, NNONSING
              IJ = (JORTN-1)*NNONSING + IORTN
*. H0
              WORK(KL_INTOP + IJ - 1) =
     &        INPROD(WORK(KL_INTMAT2+(JORTN-1)*LDET),
     &               WORK(KL_INTMAT+(IORTN-1)*LDET),LDET)
*. S
              WORK(KL_INTOP2+ IJ - 1) =
     &        INPROD(WORK(KL_INTMAT+(JORTN-1)*LDET),
     &               WORK(KL_INTMAT+(IORTN-1)*LDET),LDET)
            END DO
           END DO
          ELSE
*. Disk version
*. Determine amount of memory that can be used for storing expansion
*. of external states
C MEMMAN(KBASE,KADD,TASK,IR,IDENT)
            KBASE = 0
            CALL MEMMAN(KBASE,KFREE,'FREE  ',IDUM,'CHKFRE')
            KMAX = 100000000
            KDISC = MIN(KMAX,KFREE,LDET*LINT)
            WRITE(6,*) 
     &      ' Amount of memory to be used for batches of int. states ',
     &      KDISC
*. Number of states per batch
            NSTA_BAT = KDISC/LDET
            IF(NSTA_BAT.EQ.0) THEN
              WRITE(6,*) 
     &        ' Problem: Batch cannot contain a single internal state'
              WRITE(6,*) ' KDISC, LDET = ', KDISC, LDET
              STOP
     &      ' Problem: Batch cannot contain a single internal state'
            END IF
            NBAT = NNONSING/NSTA_BAT
            IF(NBAT*NSTA_BAT.LT.NNONSING) NBAT = NBAT + 1
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' NSTA_BAT, NBAT = ',
     &                     NSTA_BAT, NBAT 
            END IF
            CALL MEMMAN(IDUM,IDUN,'MARK  ',IDUM,'BATINT')
            CALL MEMMAN(KLCBAT,KDISC,'ADDL  ',2,'CBAT  ')
            DO JBAT = 1, NBAT
              J_INI = (JBAT-1)*NSTA_BAT + 1
              J_END = MIN(JBAT*NSTA_BAT, NNONSING)
              LBAT = J_END + 1 - J_INI
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' JBAT, J_INI, J_END = ',
     &                       JBAT, J_INI, J_END 
              END IF
*. Read states in
              CALL SKPVCD(LUSCR_INT,J_INI-1,WORK(KL_INTMAT),1,-1)
              IF(NTEST.GE.1000) WRITE(6,*) ' SKPVCD passed '
              DO JSTA = 1, LBAT
               CALL FRMDSCN(WORK(KLCBAT+(JSTA-1)*LDET),-1,-1,LUSCR_INT)
               IF(NTEST.GE.1000) 
     &         WRITE(6,*) ' FRMDSCN passed for JSTA =', JSTA
              END DO
* <I!HJ>
              CALL REWINO(LUSCR_INT2)
              DO I = 1, NNONSING
               CALL FRMDSCN(WORK(KL_INTMAT),-1,-1,LUSCR_INT2)
               IF(NTEST.GE.1000) WRITE(6,*)
     &         ' FRMDSCN, LUSCR_INT2 passed for I = ', I
               DO J = J_INI, J_END
                 IJ = (J-1)*NNONSING + I
                 WORK(KL_INTOP-1+IJ) =
     &           INPROD(WORK(KL_INTMAT),WORK(KLCBAT+(J-J_INI)*LDET),
     &            LDET)
               END DO
              END DO
*.<I!J>
              CALL REWINO(LUSCR_INT)
              DO I = 1, NNONSING
               CALL FRMDSCN(WORK(KL_INTMAT),-1,-1,LUSCR_INT)
               DO J = J_INI, J_END
                 IJ = (J-1)*NNONSING+ I
                 WORK(KL_INTOP2-1+IJ) =
     &           INPROD(WORK(KL_INTMAT),WORK(KLCBAT+(J-J_INI)*LDET),
     &           LDET)
               END DO
              END DO ! End of loop over pair of states
            END DO ! End of loop over batches of states
            CALL MEMMAN(IDUM,IDUN,'FLUSM ',IDUM,'BATINT')
          END IF ! disk version
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Metric over orthonormal states '
            CALL WRTMAT(WORK(KL_INTOP2),NNONSING,NNONSING,
     &           NNONSING,NNONSING)
            WRITE(6,*) ' H0 over orthonormal states '
            CALL WRTMAT(WORK(KL_INTOP),NNONSING,NNONSING,
     &           NNONSING,NNONSING)
          END IF
          IF(I_INT_TP.EQ.2.AND.NNONSING.NE.0) THEN
*. Diagonalize zero-Hamiltonian 
C  DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
            CALL DIAG_SYMMAT_EISPACK(WORK(KL_INTOP),WORK(KL_VECII1),
     &           WORK(KL_VECII2),NNONSING,IRETURN)
*. In output eigenvectors are in WORK(KL_INTOP) and eigenvalues are in
*. WORK(KL_VECII1)
            CALL COPVEC(WORK(KL_INTOP),X2_INT_EI_FOR_SE(IB_X2),
     &      NNONSING*NNONSING)
            IF(NTEST.GE.20) THEN
              WRITE(6,*) ' Energies of internal states '
              CALL WRTMAT(WORK(KL_VECII1),1,NNONSING,1,NNONSING)
            ELSE IF(NTEST.GE.10.AND.NNONSING.GE.1) THEN
              WRITE(6,*) ' Energy of lowest internal state '
              CALL WRTMAT(WORK(KL_VECII1),1,1,1,1)
            END IF
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) 
     &        ' Transformation from zero-order to orthonormal states'
              CALL WRTMAT(X2_INT_EI_FOR_SE(IB_X2),
     &        NNONSING,NNONSING,NNONSING,NNONSING)
            END IF
          ELSE IF(I_INT_TP.EQ.3) THEN
            IF(I_IN_OR_OUT.EQ.2) THEN
              WRITE(6,*) ' Disc version for I_INT_TP = 3 not programmed'
              STOP ' Disc version for I_INT_TP = 3 not programmed'
            END IF
*. <0!O+I [H0,OJ]|0> should be obtained and diagonalized
*. On input WORK(KL_INTOP) contains <0!O+I H0 OJ|0>.
*. Obtain <0!O+I  OJ H0 |0>
*.  H0 |0> on LUHC
            ICSPC = IREFSPC
            ISSPC = IREFSPC
            ICSM = IREFSM
            ISSM = IREFSM
            IF(I_INT_HAM.EQ.1) THEN
C?           WRITE(6,*) 
C?   &       ' One-body H0 used for internal zero-order states'
             I12 = 1
             CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
            ELSE
C?           WRITE(6,*) 
C?   &       ' Two-body H used for internal zero-order states'
            END IF
            CALL REWINO(LUC)
            CALL REWINO(LUHC)
            CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,0,0)
            IF(I_INT_HAM.EQ.1)
     &      CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
*. Generate OI H0 |ref>
            DO IORTN = 1, NNONSING
              IF(NTEST.GE.1000) WRITE(6,*) ' IORTN = ', IORTN
*. Scale to obtain orthonormal state
              FACTOR = 1.0D0/SQRT(SG_INT_EI_FOR_SE(IB_SG-1+IORTN))
              CALL COPVEC(X1_INT_EI_FOR_SE(IB_X1+(IORTN-1)*LINT),
     &        WORK(KL_VECII1),LINT)
              CALL SCALVE(WORK(KL_VECII1),FACTOR,LINT)
              ICSPC = IREFSPC
              ISSPC = IMODREF_CMBSPC
C                  SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
              CALL REWINO(LUSC51)
              CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC51,
     &             WORK(KL_VECII1),1)
              IF(NTEST.GE.10000) THEN
                WRITE(6,*) ' After sigden '
                WRITE(6,*) ' Resulting vector : '
                CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
              END IF
              CALL FRMDSCN(WORK(KL_INTMAT2+(IORTN-1)*LDET),-1,-1,LUSC51)
            END DO
*. Generate OI |ref>
            DO JORTN = 1, NNONSING
              IF(NTEST.GE.10000) WRITE(6,*) ' JORTN = ', JORTN
*. Scale to obtain orthonormal state
              FACTOR = 1.0D0/SQRT(SG_INT_EI_FOR_SE(IB_SG-1+IORTN))
              CALL COPVEC(X1_INT_EI_FOR_SE(IB_X1+(JORTN-1)*LINT),
     &        WORK(KL_VECII1),LINT)
              CALL SCALVE(WORK(KL_VECII1),FACTOR,LINT)
              ICSPC = IREFSPC
              ISSPC = IMODREF_CMBSPC
C                  SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
              CALL REWINO(LUHC)
              CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &             WORK(KL_VECII1),1)
              CALL REWINO(LUSC51)
              CALL FRMDSCN(WORK(KL_INTMAT+(JORTN-1)*LDET),-1,-1,LUSC51)
            END DO
* Obtain <ref!o+i [H0,o j]|ref> =  <ref!o+i H0 o j|ref> - 
*                                  <ref!o+i o j H0|ref>
            DO IORTN = 1, NNONSING
            DO JORTN = 1, NNONSING
              IJ = (JORTN-1)*NNONSING + IORTN
*. H0
              WORK(KL_INTOP + IJ - 1) =
     &        WORK(KL_INTOP + IJ - 1) -
     &        INPROD(WORK(KL_INTMAT2+(JORTN-1)*LDET),
     &             WORK(KL_INTMAT+(IORTN-1)*LDET),LDET)
            END DO
            END DO
*
            IF(NTEST.GE.100) THEN
              WRITE(6,*) 
     &        ' <ref!o+i [H0,o j]|ref> over orthonormal states'
              CALL WRTMAT(WORK(KL_INTOP),NNONSING,NNONSING,
     &             NNONSING,NNONSING)
            END IF
*. Diagonalize:
*. Obtain eigenvalues in WORK(KL_VECII1),lefteigenvectors in
*. X2L,righteigenvectors in X2.., work(KL_MATII1),
*- ONLY PROGRAMMED FOR REAL EIGENVALUES...
            CALL GET_LR_EIGVEC_GENMAT(WORK(KL_INTOP),NNONSING,
     &           X2_INT_EI_FOR_SE(IB_X2),X2L_INT_EI_FOR_SE(IB_X1),
     &           WORK(KL_VECII1),
     &           WORK(KL_VECII2),WORK(KL_WORK),WORK(KL_IWORK),
     &           LWORK,WORK(KL_MATII1))
C     GET_LR_EIGVEC_GENMAT(A,NDIM,VCR,VCL,VLR,
C    &           VLI,WORK,IWORK,LWORK,SCRMAT)
          END IF
*.        ^ End if I_INT_TP = 3
        END IF
*      ^ End if I_INT_TP.GE.2
        CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'INTMAT')
        IB_X1 = IB_X1 + LINT*NNONSING
        IB_X2 = IB_X2 + NNONSING*NNONSING
        IB_SG = IB_SG + NNONSING
        IB_S  = IB_S  + LINT*(LINT+1)/2
       END IF
*      ^ End if there was a nonvanishing number of EI states
       END DO
*      ^ End of loop over symmetries 
      END DO
*     ^ End of loop over external states
*- Clean-up time: Move Spin-orbital excitation  types back
      NSPOBEX_TP = NSPOBEX_TP_SAVE 
      CALL ICOPVE(WORK(KLSOBEX_SAVE),ISPOBEX,(NSPOBEX_TP+1)*NGAS*4)
      I12 = I12_SAVE
      PSSIGN = PSSIGN_SAVE
      IDC = IDC_SAVE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',1,'GT_INS')
*
      IF(I_IN_OR_OUT.EQ.2) THEN
C            FILEMAN_MINI(IFILE,ITASK)
        CALL FILEMAN_MINI(LUSCR_INT,'FREE  ')
        CALL FILEMAN_MINI(LUSCR_INT2,'FREE  ')
      END IF
*
      CALL LUCIAQEXIT('GT_INS')
      RETURN
      END 
      FUNCTION IMXMN(MAX_OR_MIN,IVEC,NELMNT)
*
*. Find largest (MAX_OR_MIN = 1) or smallest (MAX_OR_MIN = 2) element 
*. in integer vector IVEC
*
      INCLUDE 'implicit.inc'
      INTEGER IVEC(NELMNT)
*
      IVAL = IVEC(1)
      IF(MAX_OR_MIN.EQ.1) THEN
        DO I = 2, NELMNT
          IVAL = MAX(IVAL,IVEC(I))
        END DO
      ELSE
        DO I = 2, NELMNT
          IVAL = MIN(IVAL,IVEC(I))
        END DO
      END IF
*
      IMXMN = IVAL
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        IF(MAX_OR_MIN.EQ.1) THEN
          WRITE(6,*) ' Largest element found by IMXMN ', IMXMN
        ELSE
          WRITE(6,*) ' Smallest element found by IMXMN ', IMXMN
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE CAAB_T_ABSTR(ICAAB,IAOC_IN,IBOC_IN,
     &                       IAOC_OUT,IBOC_OUT,NGAS )
*
* A CAAB operator ICAAB and occupation of alpha and betastrings
* IAOC_IN, IBOC_IN are given
*
*. Find occupation of resulting strings
*
* Jeppe Olsen, March 17, 2007, trying once again to get momentum to MRCC
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICAAB(NGAS,4), IAOC_IN(NGAS),IBOC_IN(NGAS)
*. Output
      INTEGER IAOC_OUT(NGAS),IBOC_OUT(NGAS)
*
      DO IGAS = 1, NGAS
        IAOC_OUT(IGAS) = IAOC_IN(IGAS) + ICAAB(IGAS,1)-ICAAB(IGAS,3)
        IBOC_OUT(IGAS) = IBOC_IN(IGAS) + ICAAB(IGAS,2)-ICAAB(IGAS,4)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CAAB_T_ABSTR '
        WRITE(6,*) ' Input CAAB type '
        CALL WRT_SPOX_TP(ICAAB,1)
        WRITE(6,*) ' Input alpha- and beta-types'
        CALL IWRTMA(IAOC_IN,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC_IN,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      FUNCTION ISQELSUM(IVEC,NVEC,ISYM)
*
* ISYM = 0:  SUM_I IVEC(I)*IVEC(I)
* ISYM = 1:  SUM_I IVEC(I)*(IVEC(I)+1)/2
*
*. Jeppe Olsen, March 2007
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IVEC(NVEC)
*
      ISUM = 0
      IF(ISYM.EQ.0) THEN
        DO I = 1, NVEC
         ISUM = ISUM + IVEC(I)**2
        END DO
      ELSE
        DO I = 1, NVEC
         ISUM = ISUM + IVEC(I)*(IVEC(I)+1)/2
        END DO
      END IF
      ISQELSUM = ISUM
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' ISQELSUM : sum of squared elements = ', ISUM
      END IF
*
      RETURN
      END
      FUNCTION LEN_CISPC(JCMBSPC,ISYM,IPRNT)
*
* Number of dets and combinations for given sym and combination space
*
* Jeppe Olsen, obtained from LCISPC
*
* ===================
*.Input common blocks
* ===================
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cicisp.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'LNCISP')
      CALL LUCIAQENTER('LCISP')
*. Obtain types 
      ICISPC1 = ICMBSPC(1,JCMBSPC)
*. Type of alpha- and beta strings
      IF(ICISPC1.LE.NCISPC) THEN
        IATP = 1
        IBTP = 2
      ELSE
        IATP = IALTP_FOR_GAS(ICISPC1)
        IBTP = IBETP_FOR_GAS(ICISPC1)
      END IF
*
      NOCTPA =  NOCTYP(IATP)
      NOCTPB =  NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*.Local memory
      CALL MEMMAN(KLBLTP,NSMST,'ADDL  ',2,'KLBLTP')
      IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
     &CALL MEMMAN(KLCVST,NSMST,'ADDL  ',2,'KLCVST')
      CALL MEMMAN(KLIOIO,NOCTPA*NOCTPB,   'ADDL  ',2,'KLIOIO')
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
      IF(IDC.EQ.3.OR.IDC.EQ.4)
     &CALL SIGVST(WORK(KLCVST),NSMST)
*
*. Note : size of max blocks not recalculated
*. Array defining symmetry combinations of internal strings
      CALL SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*. allowed combination of types
      CALL IAIBCM(JCMBSPC,WORK(KLIOIO))
      CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,WORK(KLBLTP),WORK(KLCVST))
      CALL NGASDT(IGSOCCX(1,1,1),IGSOCCX(1,2,1),NGAS,ISYM,
     &   NSMST,NOCTPA,NOCTPB,WORK(KNSTSO(IATP)),
     &   WORK(KNSTSO(IBTP)),
     &   ISPGPFTP(1,IBSPGPFTP(IATP)),
     &   ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &   NCOMB,XNCOMB,MXS,MXSOO,WORK(KLBLTP),NTTSBL,
     &   LCOL,WORK(KLIOIO),MXSOO_AS,XMXSOO,XMXSOO_AS)
*
      LEN_CISPC = NCOMB
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CMB space and symmetry ', JCMBSPC,ISYM
       WRITE(6,*) ' Number of determinants ', LEN_CISPC
      END IF
*
      CALL LUCIAQEXIT('LCISP')
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'LNCISP')
*
      RETURN
      END
      SUBROUTINE FIND_TYPSTR_WITH_TOTOCC(NEL,ITYP)
*
* Find type with NEL electrons. The program finds 
* the last set type with given number of electrons
*
*. Jeppe Olsen, Billund airport, March 22, 2006
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'gasstr.inc'
#include "global.fh"
*
      ITYP = 0
      DO ISTTP = 1, NSTTP
        IF(NELFTP(ISTTP).EQ.NEL) ITYP = ISTTP
      END DO
*
      IF(ISTTP.EQ.0) THEN
        WRITE(6,*) 
     &  ' Error, stringtype with given number of elecs not found'
        WRITE(6,*) ' Required number of electrons = ', NEL
        STOP
     &  ' Error, stringtype with given number of elecs not found'
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Number of electrons and stringtype ',
     &  NEL, ITYP
      END IF
*
      RETURN
      END
      SUBROUTINE MATRIX_TIMES_SPARSEVEC
     &(A,VECELM,IVECIND,NDIM,NVEC,VECOUT)
*
*. VECOUT(I) = sum_(k=1,NVEC) A(I,(IVECIND(K))*VECELM(K)
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NDIM,*)
      DIMENSION VECELM(NVEC),IVECIND(NVEC)
*. Output
      DIMENSION VECOUT(NDIM)
*
      ZERO = 0.0D0
      CALL SETVEC(VECOUT,ZERO,NDIM)
*
      DO K = 1, NVEC
       KCOL = IVECIND(K)
       COEF = VECELM(K)
       DO I = 1, NDIM
        VECOUT(I) = VECOUT(I) + COEF*A(I,KCOL)
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GET_LR_EIGVEC_GENMAT(A,NDIM,VCR,VCL,VLR,
     &           VLI,WORK,IWORK,LWORK,SCRMAT)
*
*. Obtain eigenvalues and left and right eigenvectors of real 
*. general matrix - outer routine for DGEEV
*
* Has not been generalized to complex eigenvalues 
* and eigenvectors 
*
*. Jeppe Olsen, September 07
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NDIM,NDIM)
*. Output
      DIMENSION VCR(NDIM,NDIM),VCL(NDIM,NDIM) 
      DIMENSION VLR(NDIM),VLI(NDIM)
*. Scratch: WORK: Atleast 4*NDIM, preferably longer, IWORK
      DIMENSION WORK(LWORK), IWORK(LWORK)
*. And a scratch matrix - used for constructing overlap of eigenvectors
      DIMENSION SCRMAT(NDIM*NDIM)
*
      REAL*8 INPROD
      INFO = 0
      CALL DGEEV('V','V',NDIM,A,NDIM,VLR,VLI,VCL,NDIM,VCR,NDIM,
     &            WORK,LWORK,INFO)
      IF(INFO.NE.0) THEN
        WRITE(6,*) ' Error in GET_LT_EIGVEC_GENMAT'
        WRITE(6,*) ' INFO=',INFO
        STOP       ' Error in GET_LT_EIGVEC_GENMAT'
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from DGEEV '
        WRITE(6,*) ' List of eigenvalues, real and imaginary parts'
        CALL WRT_2VEC(VLR,VLI,NDIM)
      END IF
*. test if there are imaginary parts of eigenvalues
      THRES = 1.0D-10
      XMAX_IMAG = 0.0D0
      N_IMAG = 0
      DO IVAL = 1, NDIM
        IF(ABS(VLI(IVAL)).GT.THRES) N_IMAG = N_IMAG + 1
        XMAX_IMAG = MAX(XMAX_IMAG,ABS(VLI(IVAL)))
      END DO
*
      IF(N_IMAG.NE.0) THEN
        WRITE(6,*) 
     &  ' GET_LR_EIGVEC_GENMAT Complex eigenvalues encountered '
        WRITE(6,*) ' Number of imaginary eigenvalues ', N_IMAG
        WRITE(6,*) ' Largest imaginary component ', XMAX_IMAG
        STOP
     &  ' GET_LR_EIGVEC_GENMAT Complex eigenvalues encountered '
      END IF
      K1 = 1
      K2 = K1 + NDIM
      K3 = K2 + NDIM
*
      DO I = 1, NDIM
        IWORK(K1-1+I) = 0
      END DO
*. Make sure that <L(I)!R(J)> = delta(I,J) also holds 
*. degenerate eigenvectors
*
*. Loop over eigenvalues and collect degenerate sets
*
      DO I = 1, NDIM
*. Ensure that this eigenvalue has not been studied before
        IF(IWORK(K1-1+I).EQ.0) THEN
*. Check  if eigenvalue I is degenerate with other
          NDEG = 0
          DO J = 1, NDIM
            IF(IWORK(K1-1+J).EQ.0.AND.ABS(VLR(I)-VLR(J)).LE.THRES) THEN
             NDEG = NDEG + 1 
             IWORK(K2-1+NDEG) = J
            END IF
          END DO
*. The NDEG eigenvalues IWORK(K2-1+J),J=1,NDEG are degenerate, set
*. up overlap matrix in this space
          IF(NDEG.EQ.1) THEN
*. No problem- just scale L(I) so <L(I)!R(I)> = 1  
             IWORK(K1-1+I) = 1
             XLR = INPROD(VCL(1,I),VCR(1,I),NDIM)
             SCALE = 1.0D0/SCALE
             CALL SCALVE(VCL(1,I),SCALE,NDIM)
          ELSE
*. Mark that these vectors have been accessed
             DO K = 1, NDEG
               IWORK(K1-1+IWORK(K2-1+K)) = 1
             END DO
*. construct overlap matrix of degenerate vectors
             DO II = 1, NDEG
             DO J = 1, NDEG
               SCRMAT((J-1)*NDEG + II) =
     &         INPROD(VCL(1,II),VCR(1,J),NDIM)
             END DO
             END DO
             IF(NTEST.GE.100) THEN
              WRITE(6,*)  
     &        ' Overlap matrix <L(I)!R(J)> for degenerate eigvectors'
              WRITE(6,*) ' Eigenvectors: ', (IWORK(K2-1+II),II=1,NDEG)
              CALL WRTMAT(SCRMAT,NDEG,NDEG,NDEG)
             END IF
*. Find inverse of overlap matrix- use original matrix as scratch space,
*. inverse is returned in SCRMAT
             CALL INVMAT(SCRMAT,A,NDIM,NDIM,ISING)
C                 INVMAT(A,B,MATDIM,NDIM,ISING)
             IF(ISING.NE.0) THEN
              WRITE(6,*) ' Problems with matrix inversion'
              WRITE(6,*) ' I was programmed by an optimistic person'
              WRITE(6,*) ' so I continue'
             END IF
*.Transform the left eigenvectors, L'(i) = sum_k L_k S^-1 _ik, save in A
             DO II = 1, NDEG 
               DO K = 1, NDEG
                 WORK(K3-1+K) = SCRMAT((K-1)*NDEG+II)
               END DO
C     MATRIX_TIMES_SPARSEVEC(A,VECELM,IVECIND,NDIM,NVEC,VECOUT)
               CALL MATRIX_TIMES_SPARSEVEC(VCL,WORK(K3),IWORK(K2),NDIM,
     &              NDEG,A(1,II))
             END DO
*. And copy back to the matrix of lefteigenvectors
             DO  II = 1, NDEG
               ICOL = IWORK(K2-1+II)
               CALL COPVEC(A(1,II),VCL(1,ICOL),NDIM)
             END DO
          END IF
*         ^ End if degenerate space contained more than one element
        END IF
*       ^ End if I had not been used before
      END DO
*     ^ End of loop over I
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues: real and imaginary parts'
        CALL WRT_2VEC(VLR,VLI,NDIM)
*            WRT_2VEC(VEC1,VEC2,NDIM)
        WRITE(6,*) ' Space of bioorthonormal L and R eigenvectors'
        WRITE(6,*) ' (Real eigenvalues assumed...)'
        CALL WRTMAT(VCL,NDIM,NDIM,NDIM,NDIM)
        CALL WRTMAT(VCR,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_BLOCK_OF_HS_IN_INTERNAL_SPACE(
     &           IEXTP,IINTSM,I_HS,HSBLCK,I_INT_TP,I_ONLY_DIA)
*
*. Obtain block of H0(I_HS=1) or S(I_HS=2) over internal states
*. for external type IEXTP and internal symmetry IINTSM.
*  IF_I_ONLY_DIA = 1, then only the diagonal is calculated
*
* Jeppe Olsen, March 11, 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'ctcc.inc'
*. Output
      DIMENSION HSBLCK(*)
*
      CALL GET_BLOCK_OF_HS_IN_INTERNAL_SPACE_SLAVE(
     &     IEXTP,IINTSM,I_HS,HSBLCK,
     &     N_INTOP_TP,WORK(KLSOBEX),N_EXTOP_TP,
     &     WORK(KL_N_INT_FOR_EXT),WORK(KL_IB_INT_FOR_EXT),
     &     WORK(KL_I_INT_FOR_EXT),
     &     WORK(KL_N_INT_FOR_SE),WORK(KL_N_ORTN_FOR_SE),
     &     I_IN_TP,I_INT_OFF,I_EXT_OFF,I_ONLY_DIA)
*
      RETURN
      END 
      SUBROUTINE GET_BLOCK_OF_HS_IN_INTERNAL_SPACE_SLAVE(
     %           IEXTP,IINTSM,I_HS,HSBLCK,
     &           N_INTP,ISPOBEX,N_EXTP,
     &           N_IN_FOR_EX,IB_IN_FOR_EX,
     &           I_IN_FOR_EX,
     &           N_INT_FOR_SE,N_ORTN_FOR_SE,
     &           I_IN_TP,IB_INTP,IB_EXTP,I_ONLY_DIA)
*
*. Obtain block of H (I_HS=1) or S (I_HS=2) for external 
*  type IEXTP and Internal symmetry IINTSM
*
* If I_ONLY_DIA = 1, then only the diagonal is constructed and stored in 
* HSBLCK
*
* X_INT_EI_FOR_SE contains righthandside zero-order states
* whereas XL_INT_EI_FOR_SE contains lefthand side zero-order states
*
*. Jeppe Olsen, October 6, 2008
*
*. Note current version assumes that expansion of all internal 
*. states over determinants may be kept in memory..
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc' 
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'crun.inc'
*. Input: Complete matrices, i.e. matrices for all external types 
*.        and internal symmetries are used.
*. Number of internal types for given external type, base for internal 
*- types for given external types and the actual internal types for
*  given external type
      INTEGER N_IN_FOR_EX(*),IB_IN_FOR_EX(*),I_IN_FOR_EX(*)
*. Number of elementary internal operators of given sym
*. for given external type
      INTEGER ISPOBEX(NGAS,4,*)
*. Number of orthonormal internal states per symmetry and external type
      INTEGER N_ORTN_FOR_SE(NSMOB,*)
*. Number of internal states per symmetry and external types
      INTEGER N_INT_FOR_SE(NSMOB,*)
*. Output: complete matrix or diagonal
      DIMENSION HSBLCK(*)
*
      NTEST = 0
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' ---------------------------------- '
        WRITE(6,*) ' GET_BLOCK_OF_HS_IN_INTERNAL_SPACE '
        WRITE(6,*) ' ---------------------------------- '
*
        IF(I_HS.EQ.1) THEN
          WRITE(6,*) ' H-block  will be constructed'
        ELSE
          WRITE(6,*) ' S-block wil be constructed'
        END IF
        WRITE(6,*) ' External type of block: ', IEXTP
        WRITE(6,*) ' Symmetry of internal block: ', IINTSM
        IF(I_ONLY_DIA.EQ.1) 
     &  WRITE(6,*) ' Only diagonal terms calculated'
        WRITE(6,*) ' I_IN_TP, I_INT_HAM = ', I_IN_TP, I_INT_HAM
C?      WRITE(6,*) ' WORK(KINT1) = ', WORK(KINT1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',1,'GT_HS ')
*
      ICSM_SAVE = ICSM
      ISSM_SAVE = ISSM
      PSSIGN_SAVE = PSSIGN
      IDC_SAVE = IDC
      PSSIGN = 0.0D0
      IDC = 1
*
* ---------------------------------------------------------------
*. Offsets and general info for given external type and internal
*  symmetry
* ---------------------------------------------------------------
*
*. NOTE : reference space is assumed to be space 1
      IREFSPC = 1
*. Number of elementary internal states
      LINT = N_INT_FOR_SE(IINTSM,IEXTP)
*. Number  of orthonormal internal states
      LINT_ORTN = N_ORTN_FOR_SE(IINTSM,IEXTP)
C?    WRITE(6,*) ' LINT,LINT_ORTN, ', LINT,LINT_ORTN
*. Number of internal types for given external type
      L_INTP = N_IN_FOR_EX(IEXTP)
*. Offset in I_IN_FOR_EX for given external type
      IOFF = IB_IN_FOR_EX(IEXTP)
C?    WRITE(6,*) ' L_INTP, IOFF = ', L_INTP, IOFF
*. The active space
      IACT_GAS = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) IACT_GAS = IGAS
      END DO
C?    WRITE(6,*) ' The active GASpace ', IACT_GAS
*. Symmetry of internal operator times reference state
      IINTREFSM = MULTD2H(IINTSM,IREFSM)
C?    WRITE(6,*) ' IINTSM, IREFSM, IINTREFSM = ', 
C?   &             IINTSM, IREFSM, IINTREFSM
*. Two vectors of length = number of elementary internal operators
      CALL MEMMAN(KLVEC1,LINT,'ADDL  ',2,'LVEC1 ')
      CALL MEMMAN(KLVEC2,LINT,'ADDL  ',2,'LVEC1 ')
*. Obtain occupations of alpha- and beta-strings of reference space - 
*. is currently assumed to be a single space
C GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BEC
      CALL MEMMAN(KL_REF_AL,NGAS,'ADDL  ',2,'REF_AL')
      CALL MEMMAN(KL_REF_BE,NGAS,'ADDL  ',2,'REF_BE')
      CALL MEMMAN(KL_IREF_AL,NGAS,'ADDL  ',2,'IEF_AL')
      CALL MEMMAN(KL_IREF_BE,NGAS,'ADDL  ',2,'IEF_BE')
      CALL GET_REF_ALBE_OCC(IREFSPC,WORK(KL_REF_AL),WORK(KL_REF_BE))
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of alpha- and beta-strings in reference'
        CALL IWRTMA(WORK(KL_REF_AL),1,NGAS,1,NGAS)
        CALL IWRTMA(WORK(KL_REF_BE),1,NGAS,1,NGAS)
      END IF
      
*. In the following we are going to use/misuse the standard
*. spinorbital types stored as the first elements in the spinorbitalexcitation arrays.
*. Save the information usually stored there
*. Arrays
      CALL MEMMAN(KLSOBEX_SAVE,(NSPOBEX_TP+1)*NGAS*4,'ADDL ',1,'SPOXSV')
      CALL ICOPVE(ISPOBEX,WORK(KLSOBEX_SAVE),(NSPOBEX_TP+1)*NGAS*4)
      NSPOBEX_TP_SAVE = NSPOBEX_TP
*. Prepare the arrays defining calculations for this space
      NSPOBEX_TP = L_INTP
      DO II_INTP = 1, L_INTP
        I_INTP = I_IN_FOR_EX(IOFF-1+II_INTP)
C?      WRITE(6,*) ' II_INTP, I_INTP ', II_INTP, I_INTP
C?      WRITE(6,*) ' IB_INTP = ', IB_INTP
        CALL ICOPVE(ISPOBEX(1,1,IB_INTP-1+I_INTP),
     &              ISPOBEX(1,1,II_INTP),4*NGAS)
      END DO
* Type of C-coefficients will be type of reference
      ICATP = 1
      ICBTP = 2
*
      N_ALEL_C = NELFTP(ICATP)
      N_BEEL_C = NELFTP(ICBTP)
*. Find how action of internal operator changes occupation, 
*. as all operators makes identical changes of occupations, 
*. it is sufficient to look on the first operator
      IALDEL = IELSUM(ISPOBEX(1,1,1),NGAS)-IELSUM(ISPOBEX(1,3,1),NGAS)
      IBEDEL = IELSUM(ISPOBEX(1,2,1),NGAS)-IELSUM(ISPOBEX(1,4,1),NGAS)
C?    WRITE(6,*) ' IALDEL, IBEDEL = ', IALDEL, IBEDEL
*. Occupation of internal operator times reference space
      CALL CAAB_T_ABSTR(ISPOBEX(1,1,1),WORK(KL_REF_AL),WORK(KL_REF_BE),
     &                  WORK(KL_IREF_AL),WORK(KL_IREF_BE),NGAS)
*. We now have occupation of resulting strings, find corresponding
*. supergroups
      N_ALEL_S = N_ALEL_C + IALDEL
      N_BEEL_S = N_BEEL_C + IBEDEL
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of alpha-electrons in IOP*|ref> ',N_ALEL_S
        WRITE(6,*) ' Number of beta- electrons in IOP*|ref> ',N_BEEL_S
      END IF
*. Find supergroup types with these number of electrons
      ISATP = 0
      ISBTP = 0
      DO ISPGP_TP = 1, NSTTP
        IF(NELFTP(ISPGP_TP).EQ.N_ALEL_S) ISATP = ISPGP_TP
        IF(NELFTP(ISPGP_TP).EQ.N_BEEL_S) ISBTP = ISPGP_TP
      END DO
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' a-,b-supergroups of internal operator times ref(I)',
     &ISATP,ISBTP
*. Below is not neccesary
      CALL FIND_SPGRP_FROM_OCC(WORK(KL_IREF_AL),IAL_SPGP_S,ISATP)
      CALL FIND_SPGRP_FROM_OCC(WORK(KL_IREF_BE),IBE_SPGP_S,ISBTP)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' a-,b-supergroups of internal operator times ref(II)',
     &  ISATP,ISBTP
      END IF
*. End of not neccesary
*
* -------------------------------------------------------------------
* Generate the space with this occupation. This space is 
* stored as the last CI space and combination space. If the numbers 
* of these equals MXPICI I am trouble 
* -------------------------------------------------------------------
*
      IF(NCISPC.EQ.MXPICI.OR.NCMBSPC.EQ.MXPICI) THEN
        WRITE(6,*) ' Not space for an extra CI space '
        WRITE(6,*) ' Increase parameter MXPICI'
        WRITE(6,*) ' NCISPC, MXPICI ', NCISPC, MXPICI
        WRITE(6,*) ' NCMBSPC, MXPICI ', NCMBSPC, MXPICI
        STOP ' Increase parameter MXPICI'
      END IF
      IMODREF_CISPC = NCISPC + 1
      IMODREF_CMBSPC = NCMBSPC + 1
C?    WRITE(6,*) ' IMODREF_CISPC, IMODREF_CMBSPC =',
C?   &             IMODREF_CISPC, IMODREF_CMBSPC
      LCMBSPC(IMODREF_CMBSPC) = 1
      ICMBSPC(1,IMODREF_CMBSPC) = IMODREF_CISPC
*
      IALTP_FOR_GAS(IMODREF_CMBSPC) = ISATP
      IBETP_FOR_GAS(IMODREF_CMBSPC) = ISBTP
*. Occupation constraints for modified reference: electrons are 
*. added or removed from the active space. At the moment a single 
* active space is assumed.
*. Copy reference space occupation
      CALL ICOPVE(IGSOCCX(1,1,IREFSPC),IGSOCCX(1,1,IMODREF_CISPC),
     &            NGAS)
      CALL ICOPVE(IGSOCCX(1,2,IREFSPC),IGSOCCX(1,2,IMODREF_CISPC),
     &            NGAS)
*. Change number of electrons in active space
      DO IGAS = 1, NGAS
        IF(IGAS.GE.IACT_GAS) THEN
          IGSOCCX(IGAS,1,IMODREF_CISPC) =
     &    IGSOCCX(IGAS,1,IMODREF_CISPC) + IALDEL + IBEDEL
          IGSOCCX(IGAS,2,IMODREF_CISPC) =
     &    IGSOCCX(IGAS,2,IMODREF_CISPC) + IALDEL + IBEDEL
        END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Min and max for modified reference space '
        CALL IWRTMA(IGSOCCX(1,1,IMODREF_CISPC),1,NGAS,1,NGAS)
        CALL IWRTMA(IGSOCCX(1,2,IMODREF_CISPC),1,NGAS,1,NGAS)
      END IF
*. Number of determinants in modified reference space
      LDET = LEN_CISPC(IMODREF_CMBSPC,IINTREFSM,NTEST)
C     WRITE(6,*) ' LDET = ', LDET
*. Space for expansion Int |ref> in SD for all Internal states
      IDUM = 0
      IF(I_ONLY_DIA.EQ.0) THEN
        NTERMS = LINT_ORTN
      ELSE
        NTERMS = LDET
      END IF
      LDIM = LINT_ORTN*NTERMS
      CALL MEMMAN(KL_INTSDMAT, LDIM,'ADDL  ',2,'INTMAT')
      CALL MEMMAN(KL_INTSDMAT2,LDIM,'ADDL  ',2,'INTMAT')
*
      ICSPC = IREFSPC
      ISSPC = IMODREF_CISPC
      ICSM = IREFSM
      ISSM = IINTREFSM
      IF(I_HS.EQ.2) THEN
       IF(I_ONLY_DIA.EQ.0) THEN
*. Obtain complete block
*. O+(I,right)|ref> in KL_INTSDMAT
        DO INTOP = 1, LINT_ORTN
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
*. obtain internal operator INTOP
C     GET_ZERO_ORDER_STATE(IEXTP,INTSM,IORTN,X,ILR)
          CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),2)
          ICSM = IREFSM
          ISSM = IINTREFSM
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &         WORK(KLVEC1),1)
          CALL REWINO(LUHC)
          CALL FRMDSCN(WORK(KL_INTSDMAT+(INTOP-1)*LDET),-1,-1,LUHC)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' O+(I,right) |ref> in SD-basis for I = ', INTOP
            CALL WRTMAT(WORK(KL_INTSDMAT+(INTOP-1)*LDET),1,LDET,1,LDET)
          END IF
        END DO
*. O+(J.left)|ref> in KL_INTSDMAT2
        DO INTOP = 1, LINT_ORTN
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
*. obtain internal operator INTOP
C              GET_ZERO_ORDER_STATE(IEXTP,INTSM,IORTN,X,ILR)
          ICSM = IREFSM
          ISSM = IINTREFSM
          CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),1)
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &         WORK(KLVEC1),1)
          CALL REWINO(LUHC)
          CALL FRMDSCN(WORK(KL_INTSDMAT2+(INTOP-1)*LDET),-1,-1,LUHC)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' O+(I,left) |ref> in SD-basis for I = ', INTOP
            CALL WRTMAT(WORK(KL_INTSDMAT2+(INTOP-1)*LDET),1,LDET,1,LDET)
          END IF
        END DO
*. S(I,J) = <ref!op(i,left) op(j,right)!ref>
        DO I = 1, LINT_ORTN
          DO J = 1, LINT_ORTN
            HSBLCK((J-1)*LINT_ORTN + I) =
     &      INPROD(WORK(KL_INTSDMAT2+(I-1)*LDET),
     &             WORK(KL_INTSDMAT+(J-1)*LDET),LDET)
          END DO
        END DO
*
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' The overlap  <ref!op(i)+ op(j)!ref>'
         CALL PRSYM(HSBLCK,LINT)
        END IF
       ELSE IF(I_ONLY_DIA.EQ.1) THEN
*. Obtain only diagonal block
        DO INTOP = 1, LINT_ORTN
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
          CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),2)
          ICSM = IREFSM
          ISSM = IINTREFSM
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &         WORK(KLVEC1),1)
          CALL REWINO(LUHC)
          CALL FRMDSCN(WORK(KL_INTSDMAT),-1,-1,LUHC)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' O+(I,right) |ref> in SD-basis for I = ', INTOP
            CALL WRTMAT(WORK(KL_INTSDMAT),1,LDET,1,LDET)
          END IF
*. O+(I.left)|ref> in KL_INTSDMAT2
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
          CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),1)
          ICSM = IREFSM
          ISSM = IINTREFSM
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &         WORK(KLVEC1),1)
          CALL REWINO(LUHC)
          CALL FRMDSCN(WORK(KL_INTSDMAT2),-1,-1,LUHC)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' O+(I,left) |ref> in SD-basis for I = ', INTOP
            CALL WRTMAT(WORK(KL_INTSDMAT2),1,LDET,1,LDET)
          END IF
*. S(I,I) = <ref!op(i,left) op(i,right)!ref>
          HSBLCK(INTOP) = 
     &    INPROD(WORK(KL_INTSDMAT2),WORK(KL_INTSDMAT),LDET)
        END DO
*
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' The overlapdiagonal  <ref!op(i)+ op(i)!ref>'
         CALL WRTMAT(HSBLCK,LINT_ORTN)
        END IF
       END IF
*      ^ End if I_ONLY_DIA = 1
      ELSE IF(I_HS.EQ.1) THEN
       IF(I_ONLY_DIA.EQ.0) THEN
        IF(I_IN_TP.GE.2) THEN
*. Obtain zero-order Hamiltonian in internal space 
C. H0 Intop(right) |ref> for all orthonormal states
          DO INTOP = 1, LINT_ORTN
            CALL REWINO(LUHC)
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),2)
            ICSM = IREFSM
            ISSM = IINTREFSM
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KLVEC1),1)
*. zero-order Hamiltonian times Int(intop,right)|ref> in INTSDMAT
            ICSPC = IMODREF_CMBSPC
            ISSPC = IMODREF_CMBSPC
            ISSM = IINTREFSM
            ICSM = IINTREFSM
            IF(I_INT_HAM.EQ.1) THEN
C?           WRITE(6,*) 
C?   &       ' One-body H0 used for internal zero-order states'
             I12 = 1
             CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
            ELSE
C?           WRITE(6,*) 
C?   &       ' Two-body H used for internal zero-order states'
            END IF
            CALL REWINO(LUHC)
            CALL REWINO(LUSC51)
            ICSM = IINTREFSM
            ISSM = IINTREFSM
            CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC51,0,0)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' After mv7 '
              WRITE(6,*) ' Resulting vector : '
              CALL WRTVCD(WORK(KVEC1P),LUSC51,1,-1)
            END IF
            IF(I_INT_HAM.EQ.1) 
     &      CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
*
            CALL REWINO(LUSC51)
            CALL FRMDSCN(WORK(KL_INTSDMAT+(INTOP-1)*LDET),-1,-1,
     &      LUSC51)
          END DO
C. Intop(left) |ref> for all orthonormal states IN INTSDMAT2
          DO INTOP = 1, LINT_ORTN
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
            ICSM  = IREFSM
            ISSM  = IINTREFSM
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUC)
            CALL REWINO(LUHC)
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),1)
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KLVEC1),1)
            CALL REWINO(LUHC)
            CALL FRMDSCN(WORK(KL_INTSDMAT2+(INTOP-1)*LDET),-1,-1,LUHC)
          END DO
*. <ref|intop(left) H0 intop(right)|ref> 
          DO IORTN = 1, LINT_ORTN
            DO JORTN = 1, LINT_ORTN
              IJ = (JORTN-1)*LINT_ORTN + IORTN
*. H0
              HSBLCK(IJ) = 
     &        INPROD(WORK(KL_INTSDMAT+(JORTN-1)*LDET),
     &               WORK(KL_INTSDMAT2+(IORTN-1)*LDET),LDET)
            END DO
          END DO
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' H0 over orthonormal states '
            CALL WRTMAT(HSBLCK,LINT_ORTN,LINT_ORTN,LINT_ORTN,LINT_ORTN)
          END IF
        END IF
        IF(I_IN_TP.EQ.3) THEN
*. <0!O+I [H0,OJ]|0> should be obtained 
*. On input HSBLCK contains <0!O+I H0 OJ|0>.
*. Obtain <0!O+I  OJ H0 |0>
*.  H0 |0> on LUHC
          ICSPC = IREFSPC
          ISSPC = IREFSPC
          ICSM = IREFSM
          ISSM = IREFSM
          IF(I_INT_HAM.EQ.1) THEN
C?          WRITE(6,*) 
C?   &      ' One-body H0 used for internal zero-order states'
            I12 = 1
            CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
           ELSE
C?          WRITE(6,*) 
C?   &      ' Two-body H used for internal zero-order states'
          END IF
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
          CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,0,0)
          IF(I_INT_HAM.EQ.1) CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
*. Generate OI H0 |ref> and save in KL_INTSDMAT
          DO INTOP= 1, LINT_ORTN
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUSC51)
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),2)
            ICSM = IINTREFSM
            ISSM = IREFSM
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC51,
     &           WORK(KLVEC1),1)
            CALL FRMDSCN(WORK(KL_INTSDMAT+(INTOP-1)*LDET),-1,-1,LUSC51)
          END DO
*. Generate OI |ref> in KL_INTSDMAT2
          DO INTOP = 1, LINT_ORTN
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUHC)
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),1)
            ICSM = IREFSM
            ISSM = IINTREFSM
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KLVEC1),1)
            CALL REWINO(LUSC51)
            CALL FRMDSCN(WORK(KL_INTSDMAT2+(INTOP-1)*LDET),-1,-1,
     &      LUSC51)
          END DO
* Obtain <ref!o+i [H0,o j]|ref> =  <ref!o+i H0 o j|ref> - 
*                                  <ref!o+i o j H0|ref>
          DO IORTN = 1, LINT_ORTN
            DO JORTN = 1, LINT_ORTN
              IJ = (JORTN-1)*LINT_ORTN + IORTN
*. H0
              HSBLCK(IJ) = HSBLCK(IJ) -
     &        INPROD(WORK(KL_INTSDMAT2+(JORTN-1)*LDET),
     &             WORK(KL_INTSDMAT+(IORTN-1)*LDET),LDET)
            END DO
          END DO
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) 
     &      ' <ref!o+i [H0,o j]|ref> over orthonormal states'
            CALL WRTMAT(HSBLCK,LINT_ORTN,LINT_ORTN,LINT_ORTN,LINT_ORTN)
          END IF
*
        END IF
*.        ^ End if I_IN_TP = 3
       ELSE IF(I_ONLY_DIA.EQ.1) THEN
        IF(I_IN_TP.GE.2) THEN
*. Obtain diagonal of zero-order Hamiltonian in internal space 
          DO INTOP = 1, LINT_ORTN
C. H0 Intop(right) |ref> 
            CALL REWINO(LUHC)
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
            ICSM = IREFSM
            ISSM = IINTREFSM
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),2)
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KLVEC1),1)
*. zero-order Hamiltonian times Int(intop,right)|ref> in INTSDMAT
            IF(I_INT_HAM.EQ.1) THEN
C?           WRITE(6,*) 
C?   &       ' One-body H0 used for internal zero-order states'
             I12 = 1
             CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
            ELSE
C?           WRITE(6,*) 
C?   &       ' Two-body H used for internal zero-order states'
            END IF
            ICSPC = IMODREF_CMBSPC
            ISSPC = IMODREF_CMBSPC
            CALL REWINO(LUHC)
            CALL REWINO(LUSC51)
            ICSM = IINTREFSM
            ISSM = IINTREFSM
            CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC51,0,0)
            IF(I_INT_HAM.EQ.1) 
     &      CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' After mv7 '
              WRITE(6,*) ' Resulting vector : '
              CALL WRTVCD(WORK(KVEC1P),LUSC51,1,-1)
            END IF
            CALL REWINO(LUSC51)
            CALL FRMDSCN(WORK(KL_INTSDMAT),-1,-1,
     &      LUSC51)
C. Intop(left) |ref> INTSDMAT2
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
            ICSM = IREFSM
            ISSM = IINTREFSM
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),1)
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUC)
            CALL REWINO(LUHC)
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KLVEC1),1)
            CALL REWINO(LUHC)
            CALL FRMDSCN(WORK(KL_INTSDMAT2),-1,-1,LUHC)
*. <ref|intop(left) H0 intop(right)|ref> 
            HSBLCK(INTOP) = 
     &      INPROD(WORK(KL_INTSDMAT),WORK(KL_INTSDMAT2),LDET)
          END DO
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Diagonal of H0 over orthonormal states '
            CALL WRTMAT(HSBLCK,LINT_ORTN,1,LINT_ORTN,1)
          END IF
        END IF
        IF(I_IN_TP.EQ.3) THEN
*. Not prepared for I_INT_HAM = 2 !!!
*. <0!O+I [H0,OI]|0> should be obtained 
*. On input HSBLCK contains <0!O+I H0 OI|0>.
*. Obtain <0!O+I  OI H0 |0>
*.  H0 |0> on LUHC
          ICSPC = IREFSPC
          ISSPC = IREFSPC
          ICSM = IREFSM
          ISSM = IREFSM
          I12 = 1
          CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
          CALL REWINO(LUC)
          CALL REWINO(LUHC)
          CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,0,0)
          CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
*. Generate OI H0 |ref> and save in KL_INTSDMAT
          DO INTOP= 1, LINT_ORTN
C?          WRITE(6,*) ' INTOP= ', INTOP
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
            ICSM  = IREFSM
            ISSM  = IINTREFSM
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUSC51)
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),2)
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC51,
     &           WORK(KLVEC1),1)
            CALL FRMDSCN(WORK(KL_INTSDMAT),-1,-1,LUSC51)
*. Generate OI |ref> in KL_INTSDMAT2
            ICSPC = IREFSPC
            ISSPC = IMODREF_CMBSPC
C                SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
            CALL REWINO(LUHC)
            CALL GET_ZERO_ORDER_STATE(IEXTP,IINTSM,INTOP,WORK(KLVEC1),1)
            CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &           WORK(KLVEC1),1)
            CALL REWINO(LUSC51)
            CALL FRMDSCN(WORK(KL_INTSDMAT2),-1,-1,
     &      LUSC51)
* Obtain <ref!o+i [H0,o j]|ref> =  <ref!o+i H0 o j|ref> - 
            HSBLCK(INTOP) = HSBLCK(INTOP) -
     &      INPROD(WORK(KL_INTSDMAT2),WORK(KL_INTSDMAT),LDET)
          END DO
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) 
     &      ' <ref!o+i [H0,o i]|ref> over orthonormal states'
            CALL WRTMAT(HSBLCK,1,LINT_ORTN,1,LINT_ORTN)
          END IF
*
        END IF
*.        ^ End if I_IN_TP = 3
       END IF
*.     ^ End of I_ONLY_DIA switch
      END IF
*     ^
*     End if I_HS = 1
*. Clean up: restore spinorbitaltypes and symmetries
      NSPOBEX_TP = NSPOBEX_TP_SAVE 
      CALL ICOPVE(WORK(KLSOBEX_SAVE),ISPOBEX,(NSPOBEX_TP+1)*NGAS*4)
*
      ICSM = ICSM_SAVE
      ISSM = ISSM_SAVE
      PSSIGN = PSSIGN_SAVE
      IDC = IDC_SAVE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',1,'GT_HS ')
*
      RETURN
      END 
      FUNCTION E1_FOR_STRING(HDIAG,ISTRING,NEL)
*
* Obtain one-electron contribution to energy from given string
*
*. Jeppe Olsen, March 2009
*
      INCLUDE 'implicit.inc'
*
*. Input:
      DIMENSION  ISTRING(NEL)
      DIMENSION HDIAG(*)
*
      E1 = 0.0D0
      DO IEL = 1, NEL
       E1 = E1 + HDIAG(ISTRING(IEL))
      END DO
*
      E1_FOR_STRING = E1
*
      RETURN
      END
      FUNCTION GET_E0REF_EXT(FI,IPHGASX)
*
*. Obtain external part of zero-order energy.
*. State is assumed to contain a double occupied hole part
*. and an unoccupied particle part. P/H is flagged by IPHGASX)
*
*. Jeppe Olsen, March 11, 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      DIMENSION FI(*)
      DIMENSION IPHGASX(*)
*
      NTEST = 100
      WRITE(6,*) ' I am not working'
      STOP ' GET_E0REF_EXT not working'
*
      E0REF_EXT = 0.0D0
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          IB = 1
        ELSE
          IB = IB + NOBPT(IGAS)
        END IF
        IF(IPHGASX(IGAS).EQ.2) THEN
          DO I = 1, NOBPT(IGAS)
*. Absolute number of orbital in ST ordering
           I_ABS = IREOTS(IB-1+I)
           E0REF_EXT = E0REF_EXT + 2.0D0*FI((I_ABS+1)*I_ABS/2)
          END DO
        END IF
      END DO
*
      GET_E0REF_EXT = E0REF_EXT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' External part of zero-order energy ', E0REF_EXT
      END IF
*
      RETURN
      END
      FUNCTION N_ZERO_ORDER_STATES(NORTN_FOR_SE,NDIM_EX_ST,N_EXTP,
     &         ITOTSM)
*
*. Number of zero-order states of given symmetry
*
*. Jeppe Olsen, March 2009
*
      INCLUDE 'wrkspc.inc'
*. General input
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
*. Specific input
      INTEGER NORTN_FOR_SE(NSMOB,N_EXTP),NDIM_EX_ST(NSMOB,N_EXTP)
*
      NTEST = 1000
*. 
      N = 0
      DO I_EXTP = 1, N_EXTP
       DO I_EXSM = 1, NSMOB
        I_INSM =  MULTD2H(I_EXSM,ITOTSM)
        N = N + NORTN_FOR_SE(I_INSM,I_EXTP)*NDIM_EX_ST(I_EXSM,I_EXTP)
C?      WRITE(6,*) ' I_EXTP,  I_EXSM, I_INSM = ',I_EXTP,I_EXSM,I_INSM
C?      WRITE(6,*) ' NORTN_FOR_SE, NDIM_EX_ST = ',
C?   &  NORTN_FOR_SE(I_INSM,I_EXTP),NDIM_EX_ST(I_EXSM,I_EXTP)
       END DO
      END DO
*
      N_ZERO_ORDER_STATES =  N
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of zero-order states = ',
     &  N_ZERO_ORDER_STATES
      END IF
*
      RETURN
      END
      SUBROUTINE GET_ZERO_ORDER_STATE(IEXTP,INTSM,IORTN,X,ILR)
*
* Obtain internal state IORTN of symmetri INTSM and corresponding to
* external type IEXTP
*
*. Master code
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cei.inc'
*.Output
      DIMENSION X(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ', IDUM,'GT_ZOST')
      CALL MEMMAN(KLSCR,N_INT_MAX,'ADDL  ',2,'SCRINT')
*
      IF(ILR.EQ.2) THEN
        CALL GET_ZERO_ORDER_STATE_SLAVE(IEXTP,INTSM,IORTN,X,
     &       WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &       WORK(KL_X1_INT_EI_FOR_SE),WORK(KL_X2_INT_EI_FOR_SE),
     &       WORK(KL_SG_INT_EI_FOR_SE),
     &       WORK(KL_IBX1_INT_EI_FOR_SE),WORK(KL_IBX2_INT_EI_FOR_SE),
     &       WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KLSCR),
     &       NSMOB)
      ELSE
        CALL GET_ZERO_ORDER_STATE_SLAVE(IEXTP,INTSM,IORTN,X,
     &       WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &       WORK(KL_X1_INT_EI_FOR_SE),WORK(KL_X2L_INT_EI_FOR_SE),
     &       WORK(KL_SG_INT_EI_FOR_SE),
     &       WORK(KL_IBX1_INT_EI_FOR_SE),WORK(KL_IBX2_INT_EI_FOR_SE),
     &       WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KLSCR),
     &       NSMOB)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM', IDUM,'GT_ZOST')
*
      RETURN
      END
      SUBROUTINE GET_ZERO_ORDER_STATE_SLAVE(IEXTP,INTSM,IORTN,X,
     &       N_ORTN_FOR_SE,N_INT_FOR_SE,
     &       X1,X2,SG,IBX1,IBX2,IBSG,XSCR,NSMOB)
*
*. Obtain zero-order state IORTN of symmetry INTSM corresponding
*  to external type IEXTP
*
*. Jeppe Olsen, March 12, 2009
*
*. General Input
      include 'implicit.inc'
*
      INTEGER N_ORTN_FOR_SE(NSMOB,*),N_INT_FOR_SE(NSMOB,*)
      DIMENSION X1(*),X2(*),SG(*)
      INTEGER IBX1(NSMOB,*),IBX2(NSMOB,*),IBSG(NSMOB,*)
*. Scratch
      DIMENSION XSCR(*)
*. Output
      DIMENSION X(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_ZERO_ORDER_STATE_SLAVE in action '
        WRITE(6,*) ' IEXTP,INTSM = ', IEXTP,INTSM
      END IF
*
* Zero-order state  X(IORTN)(J) =(X1 SG(-1/2( X2)) K,IORTN
*
      IB_X1 = IBX1(INTSM,IEXTP)
      IB_X2 = IBX2(INTSM,IEXTP)
      IB_SG = IBSG(INTSM,IEXTP)
*  
      LINT  = N_INT_FOR_SE(INTSM,IEXTP)
      LORTN = N_ORTN_FOR_SE(INTSM,IEXTP)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' LINT, LORTN = ',LINT,LORTN
        WRITE(6,*) ' IB_X1, IB_X2, IB_SG = ', IB_X1,IB_X2,IB_SG
        WRITE(6,*) ' Blocks of Sigma, X1, X2:'
        CALL WRTMAT(SG(IB_SG),1,LORTN,1,LORTN)
        CALL WRTMAT(X1(IB_X1),LINT,LORTN,LINT,LORTN)
        CALL WRTMAT(X2(IB_X2),LORTN,LORTN,LORTN,LORTN)
      END IF
* SG(-1/2) X(2) K, IORTN
      DO K = 1, LORTN
        XSCR(K) = 
     &  1.0D0/SQRT(SG(IB_SG-1+K))*X2(IB_X2-1+(IORTN-1)*LORTN + K)
      END DO
* (X1 SG(-1/2( X2)) K,IORTN
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(X1(IB_X1),XSCR,X,LINT,LORTN,0)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Internal orthonormal zero-order state '
        WRITE(6,*) ' IEXTP, INTSM, IORTN =', IEXTP,INTSM,IORTN
        WRITE(6,*) ' Expansion in elementary operators '
        CALL WRTMAT(X,LINT,1,LINT,1)
      END IF
*
      RETURN
      END
      SUBROUTINE TRANS_EI_ORTN(T_EI,T_ORTN,ITSYM,IEO,ILR,ICOCON)
*
* Transform between elementary ei and orthonormal form of 
* vector
* 
* IEO = 1 => elementary ei to orthonormal order
* IEO = 2 => orthonormal to elementary ei order
*
* ICOCON = 1 => Covariant transformation (transform as eigenvectors)
* ICOCON = 2 => Contravariant transformation (transform to ensure
*               invariance of scalar)
*
* Zero-order state is explicitly transferred as last elements in
* T_EI and T_ORTN, respectively (if I_INCLUDE_UNI=1)
*
* Jeppe Olsen, July 2009, finished on the train to Dusseldorf, aug14
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'lucinp.inc'
*
*. Explicit input and output
*
      DIMENSION T_EI(*), T_ORTN(*)
*
      IF(IEO.EQ.2.AND.ICOCON.EQ.1) THEN
        WRITE(6,*) 
     &  ' Covariant backtransformation to elementary basis not defined'
        STOP
     &  ' Covariant backtransformation to elementary basis not defined'
      END IF
*
      CALL TRANS_EI_ORTN_SLAVE(T_EI,T_ORTN,ITSYM,IEO,ILR,N_EXTOP_TP,
     &     WORK(KL_NDIM_EX_ST),WORK(KL_NDIM_IN_SE),
     &     WORK(KL_N_ORTN_FOR_SE),NSMOB,ICOCON,I_INCLUDE_UNI)
*
      RETURN
      END 
      SUBROUTINE TRANS_EI_ORTN_SLAVE(T_EI,T_ORTN,ITSYM,IEO,ILR,N_EXTP,
     &     NDIM_EX_ST,NDIM_IN_ST,NDIM_ORT_ST,NSMOB,ICOCON,
     &     I_INCLUDE_UNI)
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'multd2h.inc'
      DIMENSION NDIM_EX_ST(NSMOB,N_EXTP),NDIM_IN_ST(NSMOB,N_EXTP),
     &          NDIM_ORT_ST(NSMOB,N_EXTP)
*. Specific input and output
      DIMENSION T_EI(*), T_ORTN(*)
*
      NTEST = 00
*
      IOFF_EI=1
      IOFF_ORT=1
      DO I_EXTP = 1, N_EXTP
       DO I_EXSM = 1, NSMOB
        I_INSM = MULTD2H(I_EXSM,ITSYM)
*
        N_EX = NDIM_EX_ST(I_EXSM,I_EXTP)
        N_IN = NDIM_IN_ST(I_INSM,I_EXTP)
        N_ORT = NDIM_ORT_ST(I_INSM,I_EXTP)
*
        CALL TRANS_EI_ORTN_BL(T_EI(IOFF_EI),T_ORTN(IOFF_ORT),
     &       I_EXTP,I_INSM,N_EX,IEO,ILR,ICOCON)
        IOFF_EI = IOFF_EI + N_EX*N_IN
        IOFF_ORT = IOFF_ORT + N_EX*N_ORT
       END DO
      END DO
*
      IF(I_INCLUDE_UNI.EQ.1) THEN
        IF(IEO.EQ.1) THEN
           T_ORTN(IOFF_ORT) = T_EI(IOFF_EI)
        ELSE
           T_EI(IOFF_EI) = T_ORTN(IOFF_ORT)
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
*
       IF(IEO.EQ.1) THEN
         WRITE(6,*) ' Transformation:  T(I,E) => T(ORT,E)'
       ELSE
         WRITE(6,*) ' Transformation:  T(ORT,E) => T(I,E)'
       END IF
*
       IF(ICOCON.EQ.1) THEN
         WRITE(6,*) ' Covariant transformation '
       ELSE 
         WRITE(6,*) ' Contravariant transformation'
       END IF
*
       WRITE(6,*) ' Vector T(I,E):'
         CALL PRINT_T_EI(T_EI,1,ITSYM)
       WRITE(6,*) ' Vector T(Ort,E):'
         CALL PRINT_T_EI(T_ORTN,2,ITSYM)
*
COLD   IF(I_INCLUDE_UNI.EQ.1)  THEN
COLD     WRITE(6,*) ' And the coefficient for unitop:'
COLD     IF(IEO.EQ.1) THEN
COLD       WRITE(6,*) T_EI(IOFF_EI)
COLD     ELSE 
COLD       WRITE(6,*) T_ORTN(IOFF_ORT)
COLD     END IF
COLD   END IF
      END IF
*
      RETURN
      END
      SUBROUTINE TRANS_EI_ORTN_BL(T_EI,T_ORTN,IEXTP,INTSM,NVEC,IEO,ILR,
     &                            ICOCON)
*
*. A block of coefficients given by external type IEXTP and 
*. internal symmetry INTSM
*. Transform between orthonormal form and elementary ei form of 
*. operators. 
* IEO = 1 => elementary ei to orthonormal order
* IEO = 2 => orthonormal to elementaty ei order
*
* ICOCON = 1 => Covariant transformation (transform as eigenvectors)
* ICOCON = 2 => Contravariant transformation (transform to ensure
*               invariance)
*
*  The elements are supposed to be in EI order
*
*. Jeppe Olsen, March 12, 2009, on the train to Koeln,
*               July 29, 2009, ICOCON added
*
* Last modification; Oct. 18, 2012; Jeppe Olsen, reduced allocation
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cei.inc'
*. Specific input and output
      DIMENSION T_EI(*), T_ORTN(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'TREIOR')
      CALL MEMMAN(KLSCR,N_INT_MAX*NVEC,'ADDL  ',2,'SCRINT')
*
      IF(ILR.EQ.2) THEN
        CALL  TRANS_EI_ORTN_EL_SLAVE
     &       (T_EI,T_ORTN,IEXTP,INTSM,NVEC,IEO,
     &       WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &       WORK(KL_X1_INT_EI_FOR_SE),WORK(KL_X2_INT_EI_FOR_SE),
     &       WORK(KL_SG_INT_EI_FOR_SE),
     &       WORK(KL_IBX1_INT_EI_FOR_SE),WORK(KL_IBX2_INT_EI_FOR_SE),
     &       WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KLSCR),
     &       NSMOB,ICOCON)
      ELSE
        CALL  TRANS_EI_ORTN_EL_SLAVE
     &       (T_EI,T_ORTN,IEXTP,INTSM,NVEC,IEO,
     &       WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &       WORK(KL_X1_INT_EI_FOR_SE),WORK(KL_X2L_INT_EI_FOR_SE),
     &       WORK(KL_SG_INT_EI_FOR_SE),
     &       WORK(KL_IBX1_INT_EI_FOR_SE),WORK(KL_IBX2_INT_EI_FOR_SE),
     &       WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KLSCR),
     &       NSMOB,ICOCON)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM,'TREIOR')
*
      RETURN
      END
      SUBROUTINE TRANS_EI_ORTN_EL_SLAVE(
     &       T_EI,T_ORTN,IEXTP,INTSM,NVEC,IEO,
     &       N_ORTN_FOR_SE,N_INT_FOR_SE,
     &       X1,X2,SG,IBX1,IBX2,IBSG,XSCR,NSMOB,ICOCON)
*
* Transform between elementary and orthonormal form of 
* a block of coefficients given by external type IEXTP and 
* internal symmetry INTSM
*
*
* IEO = 1 => elementary to orthonormal order
* IEO = 2 => orthonormal to elementaty order
*
* ICOCON = 1 => Covariant transformation (transform as eigenvectors)
* ICOCON = 2 => Contravariant transformation (transform to ensure
*               invariance)
*
*. Jeppe Olsen, March 12, 2009
*
*. Contravariant transformation:
*
* T_ORTN = X(2)T sigma(1/2) X(1)T T_EI
* T_EI   = X(1) Sigma(-1/2) X(2)  T_ORTN
*
*. Covariant transformation:
*
* V_ORTN = X(2)T sigma(-1/2)X1(T) V_EI
* Transformation from V_ORTN to V_EI is not defined or needed (I hope)





*. Generel Input
      include 'implicit.inc'
*
      INTEGER N_ORTN_FOR_SE(NSMOB,*),N_INT_FOR_SE(NSMOB,*)
      DIMENSION X1(*),X2(*),SG(*)
      INTEGER IBX1(NSMOB,*),IBX2(NSMOB,*),IBSG(NSMOB,*)
*. Scratch- should hold a block of coefficients
      DIMENSION XSCR(*)
*. Input and Output (Depending on IEO)
      DIMENSION T_EI(*), T_ORTN(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' TRANS_EI_ORTN_EL_SLAVE in action'
        WRITE(6,*) ' IEXTP,INTSM,IEO,ICOCON =', IEXTP,INTSM,IEO,ICOCON
      END IF
*
      IF(IEO.EQ.2.AND.ICOCON.EQ.1) THEN
        WRITE(6,*) 
     &  ' Covariant backtransformation to elementary basis not defined'
        STOP
     &  ' Covariant backtransformation to elementary basis not defined'
      END IF
*
      IB_X1 = IBX1(INTSM,IEXTP)
      IB_X2 = IBX2(INTSM,IEXTP)
      IB_SG = IBSG(INTSM,IEXTP)
*  
      LINT  = N_INT_FOR_SE(INTSM,IEXTP)
      LORTN = N_ORTN_FOR_SE(INTSM,IEXTP)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' LINT, LORTN = ',LINT,LORTN
        WRITE(6,*) ' IB_X1, IB_X2, IB_SG = ', IB_X1,IB_X2,IB_SG
      END IF
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
*
      IF(IEO.EQ.1) THEN
*. Elementary => orthonormal transformation
* T_ORTN = X(2)T sigma(1/2) X(1)T T_EI
*. X(1)T T_EI
       CALL MATML7(XSCR,X1(IB_X1),T_EI,LORTN,NVEC,LINT,LORTN,LINT,NVEC,
     &             FACTORC, FACTORAB,1)
* Sigma(+/- 1/2) (X(1)T T_EI) - done here columnwise to reduce number of
* sqrts - probable set up an array of sigma(-1/2) and proceed rowwise
       DO IROW = 1, LORTN
        IF(ICOCON.EQ.1) THEN
          FACTOR = 1.0D0/SQRT(SG(IB_SG-1+IROW))
        ELSE
          FACTOR = SQRT(SG(IB_SG-1+IROW))
        END IF
        DO ICOL = 1, NVEC
          XSCR((ICOL-1)*LORTN+IROW) = FACTOR*XSCR((ICOL-1)*LORTN+IROW)
        END DO
       END DO
* X(2)T (sigma(+/- 1/2) X(1)T T_EI) = X(2)T XSCR
       CALL MATML7(T_ORTN,X2(IB_X2),XSCR,
     &             LORTN,NVEC,LORTN,LORTN,LORTN,NVEC,
     &             FACTORC,FACTORAB,1)
      ELSE IF(IEO.EQ.2) THEN
*. Orthonormal to elementary order
* T_EI   = X(1) Sigma(-1/2) X(2)  T_ORTN
*
*. X(2) T_ORTN in XSCR
       CALL MATML7(XSCR,X2(IB_X2),T_ORTN,
     &             LORTN,NVEC,LORTN,LORTN,LORTN,NVEC,
     &             FACTORC,FACTORAB,0)
* Sigma(-1/2) X(2) T_ORTN = Sigma(-1/2) XSCR
       DO IROW = 1, LORTN
         FACTOR = 1.0D0/SQRT(SG(IB_SG-1+IROW))
         DO ICOL = 1, NVEC
           XSCR((ICOL-1)*LORTN+IROW) = FACTOR*XSCR((ICOL-1)*LORTN+IROW)
         END DO
       END DO
* X(1) Sigma(-1/2) X(2)  T_ORTN = X(1) XSCR in T_EI
       CALL MATML7(T_EI,X1(IB_X1),XSCR,
     &      LINT,NVEC,LINT,LORTN,LORTN,NVEC,
     &      FACTORC, FACTORAB)
      END IF
*     ^End of IEO switch
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Transformation between orthonormal and internal expansion'
        WRITE(6,*) ' IEXTP, INTSM =', IEXTP,INTSM
        IF(IEO.EQ.1) THEN
         WRITE(6,*) ' Elementary => orthonormal transformation'
        ELSE
         WRITE(6,*) ' Orthonormal=> elementary transformation'
        END IF 
        IF(ICOCON.EQ.1) THEN
          WRITE(6,*) ' Convariant transformation '
        ELSE
          WRITE(6,*) ' Contravariant transformation '
        END IF
        WRITE(6,*) ' Coefficients in elementary basis'
        CALL WRTMAT(T_EI,LINT,NVEC,LINT,NVEC)
        WRITE(6,*) ' Coefficients in orthonormal basis'
        CALL WRTMAT(T_ORTN,LORTN,NVEC,LORTN,NVEC)
      END IF
*
      RETURN
      END
      FUNCTION LARGEST_BLOCK_IN_MAT(NBLK,LR,LC)
*
* A matrix is given with NBLK blocks with row dim LR anc column dim LC
* Find largest block
*
* Jeppe Olsen, March 12, 2009
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LR(*), LC(*)
*
      LMAX = 0  
      DO I = 1, NBLK
        LMAX = MAX(LMAX,LR(I)*LC(I))
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Largest block = ', LMAX
      END IF
*
      LARGEST_BLOCK_IN_MAT = LMAX
*
      RETURN
      END
      SUBROUTINE COMPARE_TO_UNI(A,NDIM)
*
* A matrix A is given. Find largest deviation from unit matrix
* Jeppe Olsen, Aug. 2009, Red Roof Hotel, Washington
*
      INCLUDE 'implicit.inc'
      DIMENSION A(NDIM,NDIM)
*
      IOFF_R = 0
      IOFF_C = 0
      I_DIAG = 0
      XMAX_OFF = 0.0D0
      XMAX_DIAG = 0.0D0
      DO I = 1, NDIM
        DO J = 1,NDIM
          IF(I.NE.J) THEN
           IF(ABS(A(I,J)).GT.XMAX_OFF) THEN
             XMAX_OFF = ABS(A(I,J))
             IOFF_R = I
             IOFF_C = J
           END IF
          ELSE
           IF(ABS(A(I,I)-1.0D0).GT.XMAX_DIAG) THEN
             XMAX_DIAG = ABS(A(I,I)-1)
             I_DIAG = I
           END IF
          END IF
*.        ^ End of diagonal/off diagonal switch
        END DO
      END DO
*
      WRITE(6,*) ' Comparison of matrix to unit matrix '
      WRITE(6,*) ' Largest offdiagonal element : value, row, column =',
     &            XMAX_OFF,IOFF_R, IOFF_C
      WRITE(6,*) 
     &' Largest deviation of unit element from 1: value and row ',
     &  XMAX_DIAG,I_DIAG
*
      RETURN
      END 
      SUBROUTINE NORM_T_EI(T,IEO,ITSYM,XNORM_EI,IPRT)
*
* Norm of the various blocks of a  T(I,E) vector given in elementary(IEO=1) or 
* orthonormal(IEO=2) form
*
*. Jeppe Olsen, Nov. 12, 2009
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'lucinp.inc'
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION XNORM_EI(*)
*
      IF(IEO.EQ.1) THEN
        CALL NORM_T_EI_SLAVE(T,ITSYM,N_EXTOP_TP,
     &       WORK(KL_NDIM_EX_ST),WORK(KL_NDIM_IN_SE),NSMOB,XNORM_EI)
      ELSE
        CALL NORM_T_EI_SLAVE(T,ITSYM,N_EXTOP_TP,
     &       WORK(KL_NDIM_EX_ST),WORK(KL_N_ORTN_FOR_SE),NSMOB,XNORM_EI)
      END IF
*
      IF(IPRT.NE.0) THEN
        WRITE(6,*) ' Norm of T-EI vector for various E-types'
        CALL WRTMAT(XNORM_EI,1,N_EXTOP_TP,1,N_EXTOP_TP)
      END IF
*
      RETURN
      END
      SUBROUTINE NORM_T_EI_SLAVE(T,ITSYM,N_EXTP,
     &           NDIM_EX_ST,NDIM_IN_ST,NSMOB,XNORM_EI)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*
      REAL*8
     &INPROD
*. Specific input
      DIMENSION T(*)
      DIMENSION NDIM_EX_ST(NSMOB,N_EXTP),NDIM_IN_ST(NSMOB,N_EXTP)
*. Output
      DIMENSION XNORM_EI(*)
*
      IOFF = 1
      DO I_EXTP = 1, N_EXTP
       X = 0.0D0
       DO I_EXSM = 1, NSMOB
        I_INSM = MULTD2H(I_EXSM,ITSYM)
        N_EX = NDIM_EX_ST(I_EXSM,I_EXTP)
        N_IN = NDIM_IN_ST(I_INSM,I_EXTP)
C?      WRITE(6,*) ' N_EX, N_IN = ', N_EX, N_IN
        X = X +INPROD(T(IOFF),T(IOFF),N_EX*N_IN)
        IOFF = IOFF + N_EX*N_IN
       END DO
       XNORM_EI(I_EXTP) = SQRT(X)
      END DO
*
      RETURN
      END
      SUBROUTINE TEST_E
*
* Test calc of E
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      REAL*8 INPRDD
*
      DIMENSION XJ1(10000),XJ2(10000)
*
      WRITE(6,*) ' First element of INT1 = ', WORK(KINT1)
      WRITE(6,*) ' IPHGAS: '
      CALL IWRTMA(IPHGAS,1,NGAS,1,NGAS)
      CALL  MV7(XJ1,XJ2,LUC,LUHC,0,0)
*
      EREF = INPRDD(XJ1,XJ2,LUC,LUHC,1,-1)
      WRITE(6,*) ' EREF, ECORE calc in ....', EREF+ECORE, ECORE
*
      RETURN
      END

      
      
c $Id$
