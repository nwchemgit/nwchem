      SUBROUTINE DIIS_ACC(LUVEC_IN,NVEC,NPVEC,SCR,BMAT,LUVEC_OUT,
     &                    VEC1,VEC2,LUSCR,IACTCOR,IRESET)
*
* DIIS acceleration
*
* A set of corrections, e_k,  are residing  on LUVEC_IN, 
* combine these to an obtained improved approximation and 
* write improved eigenvector on LUVEC_OUT as first vectors
*
* IF IACTCOR.EQ.1 Then overwrite the last correction vector 
*                 with the actual difference between current 
*                 and previous approximations
*                 The previous approximation is assumed to be 
*                 the sum of the first NVEC-1 vectors
*
* IRESET is advice given by this routine to control routine.
* IRESET = 0 suggest no reset
* IRESET = 1 suggest reset
*
* The routine minimizes the vector
*
* sum_k w_k e_k with the constraint that sum_k w_k = 1
*
* Jeppe Olsen, Feb. 1999
*
      include 'implicit.inc'
*. Scratch for holding two blocks of vectors
      DIMENSION VEC1(*),VEC2(*)
*. Scratch Atleast length 3*(NVEC+1)**2 + NVEC+1
      DIMENSION SCR(*)
*. Input/Output
      DIMENSION BMAT(*)
*
      REAL*8 INPRDD
*. Form of files
      LBLK = -1
*
      IRESET = 0
      NVECP1 = NVEC + 1
      KLINV  = 1 
      KLSCR1 = KLINV  + NVECP1**2
      KLSCR2 = KLSCR1 + NVECP1**2
      KLW    = KLSCR2 + NVECP1**2
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to DIIS_ACC '
        WRITE(6,*) ' Number of vectors in subspace', NVEC
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vectors on file LUVEC_IN' 
        CALL REWINO(LUVEC_IN)
        DO IVEC = 1, NVEC
          CALL WRTVCD(VEC1,LUVEC_IN,0,LBLK)
        END DO
      END IF
*. Currently no reuse of B matrix so
      IF(NPVEC.GT.0) THEN
        WRITE(6,*) ' DIIS_ACC : ERROR NPVEC .ne. 0 '
        STOP       ' DIIS_ACC : ERROR NPVEC .ne. 0 '
      END IF
*
*
*. Augment matrix B = <e_i!e_j>
*. It is assumed that the B-matrix previously has been constructed
*. for NPVEC vectors and is stored in BMAT in packed form. 
      DO IVEC = NPVEC+1, NVEC
*. Copy vector IVEC to LUSCR
        CALL SKPVCD(LUVEC_IN,IVEC-1,VEC1,1,LBLK)
        CALL REWINO(LUSCR)
        CALL COPVCD(LUVEC_IN,LUSCR,VEC1,0,LBLK)
        CALL REWINO(LUVEC_IN)
        DO JVEC = 1, IVEC
          CALL REWINO(LUSCR)
          X = INPRDD(VEC1,VEC2,LUSCR,LUVEC_IN,0,LBLK)
          BMAT(IVEC*(IVEC-1)/2+JVEC) = X
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN 
         WRITE(6,*) ' Updated B = EE  matrix '
         CALL PRSYM(BMAT,NVEC)
      END IF
*
*. We will start by allowing individual weights for
*. all NVEC vectors. If the system is to loose -
*. as indicated by the first vectors getting a high 
*. weight, we will first combine the first two vectors 
*. Into a single vector, then the first three vectors in 
*. a single vector etc.
*
        I_DO_WCHECK = 1
        IF(I_DO_WCHECK.EQ.1) THEN
          NTRUNC_MX = 2
        ELSE
          NTRUNC_MX = 1
        END IF
        DO NTRUNC = 1, NTRUNC_MX
          NVEC_ACT   = NVEC-NTRUNC
          NVEC_ACTP1 = NVEC_ACT +1
*. Set up matrix A
*
*       ( B11  B12  B13 ... B1N  -1)
*       ( B21  B22  B23 ... B1N  -1)
*       .
*       .
*       ( BN1  BN2  BN3 ... BNN  -1)
*       (  -1  -1   -1  ... -1    0)
*
      CALL DIIS_AMAT(BMAT,SCR(KLSCR1),NVEC)
      CALL COPVEC(SCR(KLSCR1),BMAT,NVECP1*(NVECP1+1)/2)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) NVECP1
        WRITE(6,*) ' DIIS-Matrix before truncation:'
        CALL PRSYM(BMAT,NVECP1)
      END IF
*. Modifications connected with truncation of basis set
* B(1,1) = SUM(I=1,NTRUNC,J=1,NTRUNC) B(I,J)
      B11 = 0.0D0
      IJ = 0
      DO I = 1, NTRUNC
        DO J =1 , I-1   
         IJ = IJ + 1
         B11 = B11 + 2.0D0*BMAT(IJ)
        END DO
        IJ = IJ + 1
        B11 = B11 + BMAT(IJ)
      END DO
      BMAT(1) = B11
*. B(NTRUNC,1) = 0
      IF(NTRUNC.GT.1) THEN
        BMAT(NTRUNC*(NTRUNC-1)/2+1) = 0.0D0 
      END IF
*. B(I,1),I .gt. NTRUNC
      DO I = NTRUNC+1, NVEC
       BI1 = 0.0D0
*. B(I,1) = Sum(J=1,NTRUNC) B(I,J)
       DO J = 1, NTRUNC
         IJ = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
         BI1 = BI1 + BMAT(IJ)
       END DO
       BMAT(I*(I-1)/2+1) = BI1
      END DO
*
C     Zero all off diagonal elemnents for truncated vectors
      DO J = 2, NTRUNC
       DO I = J+1,NVECP1
         IJ = I*(I-1)/2+J
         BMAT(IJ) = 0.0D0
       END DO
      END DO
      IF(NTEST.GE.10) THEN
        WRITE(6,*) NVECP1
        WRITE(6,*) ' DIIS-Matrix after truncation:'
        CALL PRSYM(BMAT,NVECP1)
      END IF
C?    WRITE(6,*) ' BMAT after truncation '
C?    CALL PRSYM(BMAT,NVECP1)
*
      CALL TRIPAK(SCR(KLINV),BMAT,2,NVECP1,NVECP1)
*. Invert 
      CALL INVERT_BY_DIAG2(SCR(KLINV),SCR(KLSCR1),SCR(KLSCR2),
     &                    SCR(KLW),NVECP1)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Inverted DIIS-Matrix:'
        CALL WRTMAT(SCR(KLINV),NVECP1,NVECP1,NVECP1,NVECP1)
      END IF
*. Multiply inverse matrix and vector (0,0, .. ,-1)
c      CALL COPVEC(SCR(KLINV+NVEC*(NVEC+1)),SCR(KLW),NVEC)
      CALL COPVEC(SCR(KLINV+NVECP1*(NVECP1-1)),SCR(KLW),NVEC)
*. Copy the weight from the first vector to vectors 2-NTRUNC 
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Weight vector (1)'
       CALL WRTMAT(SCR(KLW),1,NVEC,1,NVEC) 
      END IF
      DO I = 2, NTRUNC
        SCR(KLW-1+I) = SCR(KLW-1+1)
      END DO
*
      ONEM = -1.0D0
      CALL SCALVE(SCR(KLW),ONEM,NVEC)
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Weight vector '
       CALL WRTMAT(SCR(KLW),1,NVEC,1,NVEC) 
      END IF
*. Is this weight vector okay ?
*. Current criteria sum_i=1,NVEC-1 Abs(W_I) .le. ABS(W_NVEC)
      WPREV = 0.0D0
      DO I = 1, NVEC-1
       WPREV = WPREV + ABS(SCR(KLW-1+I))
      END DO
      FACTOR = 1.0D0
      WNVEC = ABS(SCR(KLW-1+NVEC))
      IF(WPREV.LE.FACTOR*WNVEC.OR.NVEC.LE.3.OR.I_DO_WCHECK.EQ.0)THEN
        I_AM_OKAY = 1
        IF(NTRUNC.EQ.1) IRESET = 0
      ELSE
        I_AM_OKAY = 0
        IF(NTRUNC.EQ.1) IRESET = 1
      END IF
      IF(NTEST.GE.1) WRITE(6,'(A,I3,2E13.6)') 
     &' DIIS : nvec, w_prev, w_nvec',NVEC,WPREV,WNVEC
      IF(I_AM_OKAY.EQ.1) GOTO 1001
      END DO
*     ^ End of loop over truncation level
 1001 CONTINUE
*. Form new vector 
*. First we have only the corrections stored on LUVEC_IN, 
*. weights are given as som of sum of solution vectors, so 
*. change weights
*. sum w_k x_k = sum (n+1-k)w_k e_k
*
      DO I = 1, NVEC
        SUM = 0
        DO K = I, NVEC
          SUM = SUM + SCR(KLW-1+K)
        END DO
        SCR(KLINV-1+I) = SUM
      END DO
c Alternative: add weighted correction vector to old one
c ----> 
c      DO I = 1, NVEC-1
c        SCR(KLINV-1+I) = SCR(KLW-1+I)+1d0
c      END DO
c      SCR(KLINV+NVEC-1)=SCR(KLW+NVEC-1)
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Weight vector used:'
       CALL WRTMAT(SCR(KLINV),1,NVEC,1,NVEC)         
      END IF

      IF(IACTCOR.EQ.0) THEN
       IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Weight vector used: (2a)'
         CALL WRTMAT(SCR(KLINV),1,NVEC,1,NVEC)         
       END IF
        CALL MVCSMD(LUVEC_IN,SCR(KLINV),LUVEC_OUT,LUSCR,
     &              VEC1,VEC2,NVEC,1,LBLK)
      ELSE 
*. Obtain correction vector
       DO IVEC = 1, NVEC - 1
         SCR(KLINV-1+IVEC) = SCR(KLINV-1+IVEC) -1.0D0
       END DO 
       IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Weight vector used: (2b)'
         CALL WRTMAT(SCR(KLINV),1,NVEC,1,NVEC)         
       END IF
       CALL MVCSMD(LUVEC_IN,SCR(KLINV),LUVEC_OUT,LUSCR,
     &              VEC1,VEC2,NVEC,1,LBLK)
*. And save on LUVEC_IN
       CALL SKPVCD(LUVEC_IN,NVEC-1,VEC1,1,LBLK)
       CALL REWINO(LUVEC_OUT)
       CALL COPVCD(LUVEC_OUT,LUVEC_IN,VEC1,0,LBLK)
*. And obtain new total vector on LUVEC_OUT
       ONE = 1.0D0
       CALL SETVEC(SCR(KLINV),ONE,NVEC)
       CALL MVCSMD(LUVEC_IN,SCR(KLINV),LUVEC_OUT,LUSCR,
     &              VEC1,VEC2,NVEC,1,LBLK)
      END IF
*     ^ End of IACTCOR switch 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' New approximation to solution vector'
        CALL WRTVCD(VEC1,LUVEC_OUT,1,LBLK)
      END IF
*
      RETURN
      END 
      SUBROUTINE DIIS_AMAT(EE,A,NVEC)
*
* Construct B matrix of DIIS method from error matrix
*
* Input :
* 
* NVEC : Number of vectors in DIIS subspace
* EE(I,J) =  <e_i!e_j>
*
* Output
*
*   A =  
*
*       ( E11  E12  E13 ... E1N  -1)
*       ( E21  E22  E23 ... E1N  -1)
*       .
*       .
*       ( EN1  EN2  EN3 ... ENN  -1)
*       (  -1  -1   -1  ... -1    0)
*
* Input and output matrices are assumed in lower diag packed form
*
* Jeppe Olsen, Feb28, 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION EE(*)
*. Output
      DIMENSION A(*)
*
      LEN = NVEC*(NVEC+1)/2
      CALL COPVEC(EE,A,LEN)
      ONEM = -1.0D0
      CALL SETVEC(A(LEN+1),ONEM,NVEC)
      NVECP1 = NVEC+1
      A(NVECP1*(NVECP1+1)/2) = 0.0D0
*
      NTEST = 00
      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' A matrix from DIIS_AMAT'
        CALL PRSYM(A,NVECP1)
      END IF
*
      RETURN
      END
      SUBROUTINE ADDVEC_DIIS_SUBSPC(MAXVEC,NVECIN,NVECADD,
     &           LUVEC,LUADD,LUSCR1,LUSCR2,VEC1,VEC2,NVECOUT,NVECDEL,
     &           IRESET)
*
* Add NVECADD vectors to DIIS subspace. If the number of vectors in subspace 
* before addition is MAXVEC, remove also the oldest correction
*
* LUVEC : Initial approximation to solution vector and 
*         NVECIN - 1 correction correction
*
* IF IRESET = 1, then the LUVEC file is reset : current approximation
* is constructed on LUVEC and the NVECADD vectors are added
* Jeppe Olsen, Feb. 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
*. Local scratch : Assuming max dim of subspace is 100
      PARAMETER (MXP_DIM_SUBSPACE = 100)
      DIMENSION SCR(MXP_DIM_SUBSPACE)
*
      IF(MAXVEC.GT.MXP_DIM_SUBSPACE) THEN
        WRITE(6,*) ' Potential problem for ADDVEC_DIIS_SUBSPC '
        WRITE(6,*) ' Max dim of subspace larger than program limit'
        WRITE(6,*) ' MAXVEC, MXP_DIM_SUBSPACE = ', 
     &               MAXVEC, MXP_DIM_SUBSPACE
        WRITE(6,*) ' Increase  MXP_DIM_SUBSPACE or decrease MAXVEC'
        STOP       ' Increase  MXP_DIM_SUBSPACE or decrease MAXVEC'
      END IF
*
      LBLK = -1
*
      IF(IRESET.EQ.0) THEN
        NVECOUT = MIN(MAXVEC,NVECIN+NVECADD)
        NVECDEL = MAX(0,NVECIN+NVECADD-MAXVEC)
      ELSE 
        NVECOUT = NVECADD + 1
        NVECDEL = NVECIN - 1 
      END IF
C?    WRITE(6,*) ' NVECIN, NVECADD,NVECOUT, NVECDEL',
C?   &             NVECIN, NVECADD,NVECOUT, NVECDEL
      IF(NVECDEL.EQ.0) THEN
*. Just add vectors
       CALL SKPVCD(LUVEC,NVECIN,VEC1,1,LBLK)
       CALL REWINO(LUADD)
       DO IVEC = 1, NVECADD
         CALL COPVCD(LUADD,LUVEC,VEC1,0,LBLK)
       END DO
      ELSE IF( NVECDEL.GT.0 ) THEN
*. Reset and add 
*. New first vector (initial approximation) 
*  should now by sum of first NVECDEL+1 vectors 
*. on LUVEC  
       NVECDELP1 = NVECDEL + 1
       CALL REWINO(LUVEC)
       ONE = 1.0D0
       CALL SETVEC(SCR,ONE,NVECDELP1) 
C           MVCSMD(LUIN, FAC,LUOUT, LUSCR, VEC1,VEC2,NVEC,IREW,LBL
       CALL MVCSMD(LUVEC,SCR,LUSCR1,LUSCR2,VEC1,VEC2,NVECDELP1,1,LBLK)
*. Copy to LUSCR1 the vectors that should be kept after new vector 
       CALL SKPVCD(LUVEC,NVECDELP1,VEC1,1,LBLK)
       DO IVEC = NVECDELP1+1,NVECIN
         CALL COPVCD(LUVEC,LUSCR1,VEC1,0,LBLK)
       END DO
       CALL REWINO(LUSCR1)
       CALL REWINO(LUVEC)
       DO IVEC = 1, NVECIN-NVECDEL
         CALL COPVCD(LUSCR1,LUVEC,VEC1,0,LBLK)
       END DO
*. And new vectors from LUADD
       CALL REWINO(LUADD)
       DO IVEC = 1, NVECADD
         CALL COPVCD(LUADD,LUVEC,VEC1,0,LBLK)
       END DO
      END IF
*     ^ End of switch : NVECIN
      RETURN
      END
      SUBROUTINE INVERT_BY_DIAG2(A,B,SCR,VEC,NDIM)!,IMOD,NSKIP)
*
* Invert symmetric  - hopefully nonsingular - matrix A 
* by diagonalization
*
* Jeppe Olsen, Oct 97 to check INVMAT
*              March 00 : Scale initial matrix to obtain unit diagonal
*
* Modification: Eliminate automatically vectors (starting at pos. 1)
*               until A becomes nonsingular (IMOD=1); the number
*               of skipped vectors is returned on NSKIP
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output matrix
      DIMENSION A(*)       
*. Scratch matrices and vector
      DIMENSION B(*),SCR(*),VEC(*)
*
      LOGICAL INVERT_IT
*
      NTEST = 00
cc      NSKIP = 0
cc      INVERT_IT=.TRUE.
cc
cc      DO WHILE (INVERT_IT)
cc
cc        IF (NSKIP.GE.NDIM) STOP ' Logical error!'

*. Reform a to symmetric packed form
      CALL TRIPAK(A,SCR,1,NDIM,NDIM)
cc        DO I = NSKIP+1,NDIM
cc         DO J = NSKIP+1, I-1
cc            A(I-NSKIP,J-NSKIP) = 
*. Extract diagonal
      CALL COPDIA(SCR,VEC,NDIM,1)
*
*.scale 
*
      DO I = 1, NDIM
*. Scaling vector
        IF(VEC(I).EQ.0.0D0) THEN
          VEC(I) = 1.0D0
        ELSE
          VEC(I) = 1.0D0/SQRT(ABS(VEC(I)))
        END IF
      END DO
*. Scale matrix
      IJ = 0
      DO I = 1, NDIM
        DO J = 1, I
          IJ = IJ + 1
          SCR(IJ) = SCR(IJ)*VEC(I)*VEC(J)
        END DO
      END DO
C     DO I = 1, NDIM
C       VEC(I) = 1.0D0/VEC(I)
C     END DO
*. Diagonalize
      CALL EIGENL(SCR,B,NDIM,0,1)
*. Scale eigenvectors
      DO IVEC = 1, NDIM 
        IOFF = 1 + (IVEC-1)*NDIM
        CALL VVTOV(B(IOFF),VEC,B(IOFF),NDIM)
      END DO
*. 
      CALL COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of scaled matrix : '
        CALL WRTMAT(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements 
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode activated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END 
      SUBROUTINE APRJAC_TV(NVEC,LUIN,LUOUT,LUVEC,LUJVEC,
     &           LUJDIA,VEC1,VEC2,SCR,N_CC_AMP,LUSCR,LUSCR2,
     &           MAXVEC)
*
* An approximate Jacobian is given in the form of 
* a diagonal approximation, JACDIA, and  
* NVEC vectors(in LUVEC)  and Jacobian times these vectors(LUJVEC).
*
* Find Inverse approximate Jacobian times the vector in luin to
* obtain vector in luout.
* Largest number vectors in subspave of approximate in Jacobian is MAXVEC
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
*. Local vectors
      DIMENSION VEC1(*),VEC2(*)
*. Scratch space : Should atleast be length :
      DIMENSION SCR(*)
*.
      REAL*8 INPRDD
*
      NTEST = 1000
      LBLK = -1
*. 
*. Obtain overlap of input vectors
*
      KLS = 1
      KLFREE = KLS + NVEC**2
*
      KLS2 = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLEIGVEC = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLU = KLFREE        
      KLFREE = KLFREE + NVEC**2
*
      KLSCR = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLVEC = KLFREE
      KLFREE = KLFREE + NVEC 
*
      KLVEC2 = KLFREE
      KLFREE = KLFREE + NVEC 
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'APRJA_TV reporting to work '
        WRITE(6,*) '==========================='
        WRITE(6,*) ' LUVEC, LUJVEC ',LUVEC,LUJVEC
        WRITE(6,*) ' NVEC, MAXVEC = ', NVEC,MAXVEC
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vectors on file LUVEC '
        CALL REWINO(LUVEC)
        DO IVEC = 1, NVEC
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
          CALL WRTVCD(VEC1,LUVEC,0,LBLK)
        END DO
        WRITE(6,*) ' Vectors in file LUJVEC '
        CALL REWINO(LUJVEC)
        DO IVEC = 1, NVEC
          CALL WRTVCD(VEC1,LUJVEC,0,LBLK)
        END DO
      END IF
*
      CALL REWINO(LUVEC)
      DO JVEC = 1, NVEC
        CALL REWINO(LUSCR)
        CALL SKPVCD(LUVEC,JVEC-1,VEC1,1,LBLK)
        CALL COPVCD(LUVEC,LUSCR,VEC1,0,LBLK)
C  INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
        SJJ = INPRDD(VEC1,VEC2,LUSCR,LUSCR,1,LBLK)
        JJ = JVEC*(JVEC+1)/2
        SCR(KLS-1+JJ) = SJJ
        DO IVEC = JVEC+1,NVEC
          IJ = IVEC*(IVEC-1)/2+JVEC
          CALL REWINO(LUSCR)
          SIJ = INPRDD(VEC1,VEC2,LUSCR,LUVEC,0,LBLK)
          SCR(KLS-1+IJ) = SIJ
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix '
        CALL PRSYM(SCR(KLS),NVEC)
      END IF
*
*
*. Check for linear dependencies by diagonalizing
*. Eliminate vectors ( start with the oldest) to obtain 
*. a non-singular basis
*
*
*. Allowed ratio between largest and smallest eigenvalue
      XMXRATIO = 1.0D10
*. for the case where the are no vectors in subspace 
      NVEC_EFF = 0
*
      DO NELI = 0, NVEC-1
       I1VEC = NELI + 1
       NVEC_EFF = NVEC - NELI 
       LEN = NVEC_EFF*(NVEC_EFF+1)/2
*. Obtain overlap of last NVEC_EFF vectors
       DO I = NELI+1,NVEC
         DO J = NELI+1,I
          IJ = (I-NELI)*(I-NELI-1)/2 + J - NELI
          SCR(KLS2-1+IJ) = SCR(KLS-1+I*(I-1)/2  + J) 
         END DO
       END DO
*. Diagonalize 
C      CALL EIGEN(A,R,N,MV,MFKR)
       CALL EIGENL(SCR(KLS2),SCR(KLEIGVEC),NVEC_EFF,0,1)
*. Extract eigenvalues
       CALL COPDIA(SCR(KLS2),SCR(KLS2),NVEC_EFF,1)
*. check eigenvalue ratio and explicit singularty
       I_AM_FINE = 1
       IF(SCR(KLS2).LE.0) THEN
         I_AM_FINE = 0
       ELSE      
         RATIO  = SCR(KLS2-1+NVEC_EFF)/SCR(KLS2-1+1) 
         IF(RATIO.GT.XMXRATIO) I_AM_FINE = 0
       END IF
       IF(I_AM_FINE.EQ.1) GOTO 1001
      END DO
 1001 CONTINUE
      WRITE(6,*) ' Eigenvalues of nonsingular Metric  : '
      CALL WRTMAT(SCR(KLS2),1,NVEC_EFF,1,NVEC_EFF)
*     ^ End of loop over number of vectors to be truncated
*. We now know that the last NVEC_EFF vectors spans an     
*. orthonormal basis
*
      NVEC_SUB = MIN(MAXVEC,NVEC_EFF)
      NVEC_SKIP = NVEC - NVEC_SUB
      I1VEC = NVEC-NVEC_SUB + 1
      IF(NVEC_SUB.EQ.0) GOTO 1002
*
      WRITE(6,*) ' APRJAC_TV : nvec_sub, nvec_skip ',
     &             NVEC_SUB, NVEC_SKIP
      CALL COP_SYMMAT(SCR(KLS),SCR(KLS2),I1VEC,NVEC_SUB)
*. Metric resides now in SCR(KLS2) in packed form
*
*. Obtain inverse of metrix 
*
*. 
      CALL TRIPAK(SCR(KLS),SCR(KLS2),2,NVEC_SUB,NVEC_SUB)
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
      CALL INVERT_BY_DIAG(SCR(KLS),SCR(KLEIGVEC),SCR(KLSCR),
     &                    SCR(KLVEC),NVEC_SUB)
*. And pack it
      CALL TRIPAK(SCR(KLEIGVEC),SCR(KLS),1,NVEC_SUB,NVEC_SUB)
      LEN = NVEC_SUB*(NVEC_SUB+1)/2
      CALL COPVEC(SCR(KLEIGVEC),SCR(KLS),LEN)
*. Inverse matrix resides now in SCR(KLS) in packed form
  
*
* We will use vectors I1VEC to NVEC to span the 
* subspace used for 
*
* The projection to the subspace is   
*
* P = sum_{ij} x_i s_{ij}-1 x^t_j
*
* The approximate Jacobian is thus (J0 is diagonal part of Jacobian)
*
* J_apr = J0 (1-P) + Jex P
*       = J0 + (Jex-J0) P
*       = J0 + sum_ij (s_i - s0_i) S_{ij}-1 x^t_j (s_i = Jex x_i,s0_i = J0 x_i)
*       = J0  + R S^{-1} T^t
*
* 
* where R is the collection of NVEC_SUB column vectors (s_i - s0_i)
* and T is the collection of NVEC_SUB column vectors x_i.
*
* the Sherman-Morrison formula gives the inverse of a rank-n updated matrix
* 
* A' = A + R X T^t 
*
* as
*
*   A'^{-1} = A^{-1} -  A^{-1}R U^{-1} T^t A^{-1}
*
* U = X^{-1} + T^t A^{-1} R
*
* 
* In our case the U matrix become
* U_ij =  S - T^t_i A^{-1} R_j = 
*          x_i^t J_0^{-1} J x_j 
*
      CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO J = 1,NVEC_SUB
*. Diagonal times Jx_j
        ZERO = 0.0D0
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
        CALL REWINO(LUJDIA)
        CALL REWINO(LUSCR)
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUJVEC,LUSCR,ZERO,0,1,LBLK)
        CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
        DO I = 1, NVEC_SUB
*. x_i^t J0^{-1}Jx_j 
          CALL REWINO(LUSCR)
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
          UIJ = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
          IJEXP = (J-1)*NVEC_SUB + I
          SCR(KLU-1+IJEXP) = UIJ
CERROR. and remember the overlap term
CERROR          IJSYM = I*(I-1)/2 + J
CERROR          SCR(KLU-1+ IJEXP) = SCR(KLU-1+IJEXP)
CERROR     &                      + SCR(KLS-1+IJSYM) - SCR(KLS2-1+IJSYM)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' U matrix '
        CALL WRTMAT(SCR(KLU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*. Inverse U matrix and keep in KLU
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
      CALL INVERT_BY_DIAG(SCR(KLU),SCR(KLEIGVEC),SCR(KLSCR),
     &                    SCR(KLVEC),NVEC_SUB)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverted U matrix '
        CALL WRTMAT(SCR(KLU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*
* And then J0{-1} R U{-1} T^t J0^{-1} VEC
*
*. 1 : J0^{-1} VEC
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,ZERO,1,1,LBLK)
*. 2 : T^t J0^{-1} VEC 
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO I = 1, NVEC_SUB
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
        CALL REWINO(LUSCR)
        SCR(KLVEC-1+I) = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*. 3 :   VEC' = U{-1} T^t J0^{-1} VEC (Vec' is a vector in the subspace)
      CALL MATVCB(SCR(KLU),SCR(KLVEC),SCR(KLVEC2),NVEC_SUB,
     &            NVEC_SUB,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' U(-1) T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC2),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*  4 : J0^{-1} R U{-1} T^t J0^{-1} VEC
*      = J0^{-1} ( sum_j Jx_j vec_j)  - sum_j x_j vec_j
*
*. sum_j Jx_j vec_j on LUSCR
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  LUSCR before call to MVCSMD'
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
      CALL REWINO(LUSCR)
      CALL REWINO(LUSCR2)
      CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUJVEC,SCR(KLVEC2),LUSCR,LUSCR2,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*.  J0^{-1}  sum_j Jx_j vec_j on LUSCR2
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUSCR,LUSCR2,ZERO,1,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' J0(-1) (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR2,1,LBLK)
      END IF
*.  sum_j x_j vec_j on LUSCR (LUOUT is used as scratch)
      CALL REWINO(LUSCR)
      CALL REWINO(LUOUT)
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUVEC,SCR(KLVEC2),LUSCR,LUOUT,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' xU(-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*. and the synthesis in LUOUT
C  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK
      ONE = 1.0D0
      ONEM = -1.0D0
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR2,LUSCR,LUOUT,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Total low rank correction'
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
*. And the remainders 
*. J0^{-1} Vecin on LUSCR
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,ZERO,1,1,LBLK)
*. (J0^{-1} -  J0^{-1} R U{-1} T^t J0^{-1} ) VECIN
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR,LUOUT,LUSCR2,1,LBLK)
*. Well the result ended up on LUSCR2, place it properly on LUOUT
C      COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      CALL COPVCD(LUSCR2,LUOUT,VEC1,1,LBLK)
*
 1002 CONTINUE
      IF(NVEC_SUB.EQ.0) THEN
        ZERO = 0.0D0
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUOUT,ZERO,1,1,LBLK)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector '
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUIN,1,LBLK)
        WRITE(6,*) ' Inverted approximate Jacobian times Vector '
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE APRJAC_TV2(NVEC,NNEWVEC,LUIN,LUOUT,LUVEC,
     &           LUJVEC,LUJDIA,VEC1,VEC2,
     &           SMAT,UMAT,SCR,
     &           N_CC_AMP,LUSCR,LUSCR2,MAXVEC,MAXVECMAX,
     &           XMXSTP,XMXSTP_APRJ)
*
* An approximate Jacobian is given in the form of 
* a diagonal approximation, JACDIA, and  
* NVEC vectors(in LUVEC)  and Jacobian times these vectors(LUJVEC).
*
* Find Inverse approximate Jacobian times the vector in luin to
* obtain vector in luout.
* Largest number vectors in subspave of approximate in Jacobian is MAXVEC
*
* Jeppe Olsen, March 2000
*
* Modified version: reduced I/O to ord(3*NVEC)
*  Andreas, Jan 2004
*
      INCLUDE 'implicit.inc'
*. Local vectors
      DIMENSION VEC1(*),VEC2(*)
*. Scratch space : Should atleast be length :
      DIMENSION UMAT(*),SMAT(*),SCR(*)
*. Local mini-scratch-space:
      DIMENSION IPIV(NVEC)
*. A SAVE variable to remember skipped vectors:
      INTEGER, SAVE :: NVEC_SKIP=0
      REAL*8, SAVE :: DAMP_OLD=-999.999d0, STEPNORM_LAST=1000d0
*.
      LOGICAL FACT_U, TOO_LARGE_STEP
      REAL*8 INPRDD
*
      NTEST = 150
      LBLK = -1
*. 
*. Obtain overlap of input vectors
*
      KLFREE = 1
*
      KLSU = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLEIGVEC = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLSCR = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLVEC = KLFREE
      KLFREE = KLFREE + NVEC 
*
      KLVEC2 = KLFREE
      KLFREE = KLFREE + NVEC 
*
* dirty, but we need to rewrite this subroutine soon
      MAXVECMAX=100
*
      IF(NNEWVEC.LE.0) THEN
        WRITE(6,'(x,a,i2,a)')
     &   'ARPRJA_TV: The input for NNEWVEC',NNEWVEC, 
     &   'is definitely silly! It is actually hard to continue....'
        STOP 'NNEWVEC is silly'
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'APRJA_TV2 reporting to work '
        WRITE(6,*) '==========================='
        WRITE(6,*) ' LUIN,  LUOUT         ', LUIN, LUOUT
        WRITE(6,*) ' LUVEC, LUJVEC, LUJDIA ',LUVEC,LUJVEC,LUJDIA
        WRITE(6,*) ' LUSCR, LUSCR2         ',LUSCR,LUSCR2
        WRITE(6,*) ' NVEC, NNEWVEC, MAXVEC = ', NVEC,NNEWVEC,MAXVEC 
        WRITE(6,*) ' XMXSTP, XMXSTP_APRJ = ',XMXSTP,XMXSTP_APRJ
      END IF
      IF(NTEST.GE.150) THEN
        XNRM = SQRT(INPRDD(VEC1,VEC2,LUIN,LUIN,1,LBLK))
        WRITE(6,*) 'Contents of LUIN: norm = ',XNRM
        XNRM = SQRT(INPRDD(VEC1,VEC2,LUJDIA,LUJDIA,1,LBLK))
        WRITE(6,*) 'Contents of LUJDIA: norm = ',XNRM 
        WRITE(6,*) 'Contents of LUVEC:'
        CALL REWINO(LUVEC)
        DO II = 1, NVEC
          XNRM = SQRT(INPRDD(VEC1,VEC2,LUVEC,LUVEC,0,LBLK))
          WRITE(6,*) II,'. ',XNRM
        END DO  
        WRITE(6,*) 'Contents of LUJVEC:'
        CALL REWINO(LUJVEC)
        DO II = 1, NVEC
          XNRM = SQRT(INPRDD(VEC1,VEC2,LUJVEC,LUJVEC,0,LBLK))
          WRITE(6,*) II,'. ',XNRM
        END DO  
      END IF
*
* Initialize NVEC_SKIP
      IF (NVEC.LE.1) NVEC_SKIP=0
      IF (NVEC.LE.1) DAMP_OLD=-999.999d0
      IF (NVEC.LE.1) STEPNORM_LAST=1000d0
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vectors on file LUVEC '
        CALL REWINO(LUVEC)
        DO IVEC = 1, NVEC
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
          CALL WRTVCD(VEC1,LUVEC,0,LBLK)
        END DO
        WRITE(6,*) ' Vectors in file LUJVEC '
        CALL REWINO(LUJVEC)
        DO IVEC = 1, NVEC
          CALL WRTVCD(VEC1,LUJVEC,0,LBLK)
        END DO
      END IF
    
      ! we first probe the simple perturbation correction:
      TOO_LARGE_STEP = .TRUE.
      DAMP = 0d0
      DO WHILE(TOO_LARGE_STEP)
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,DAMP,1,1,LBLK)
        CALL REWINO(LUSCR)
        XNORM = SQRT(INPRDD(VEC1,VEC2,LUSCR,LUSCR,1,LBLK))
        WRITE(6,*) 'NORM,DAMP : ',XNORM, DAMP
        IF (XNORM.GT.XMXSTP) THEN
          DAMP = DAMP+0.1d0
        ELSE
          TOO_LARGE_STEP=.FALSE.              
        END IF
      END DO

c      IF (DAMP.GT.0d0 .OR. XNORM.GT.XMAXSTEP2) THEN
c if the step is too large:
      IF (XNORM.GT.XMXSTP_APRJ) THEN
c flush subspace expansion
         NVEC_SKIP=NVEC+1
      END IF

      ! normally only the new vectors' contribution to U is needed:
      NDOVEC = NNEWVEC 
      NVEC_SUB = MAX(0,NVEC-NVEC_SKIP)
      IF (NVEC_SUB.GT.MAXVEC) THEN
        NVEC_SUB=MAXVEC
        NVEC_SKIP=NVEC-NVEC_SUB
      END IF
      IF (DAMP.NE.DAMP_OLD) NDOVEC=NVEC_SUB

      DAMP_OLD = DAMP

      IF (NTEST.GE.100) WRITE(6,*) 'NDOVEC, NVECSUB: ',NDOVEC, NVECSUB

*
* We will use vectors NVEC_SKIP+1 to NVEC to span the 
* subspace used for the projection to the subspace which is   
*
* P = sum_{ij} x_i s_{ij}-1 x^t_j
*
* The approximate Jacobian is thus (J0 is diagonal part of Jacobian)
*
* J_apr = J0 (1-P) + Jex P
*       = J0 + (Jex-J0) P
*       = J0 + sum_ij (s_i - s0_i) S_{ij}-1 x^t_j (s_i = Jex x_i,s0_i = J0 x_i)
*       = J0  + R S^{-1} T^t
*
* 
* where R is the collection of NVEC_SUB column vectors (s_i - s0_i)
* and T is the collection of NVEC_SUB column vectors x_i.
*
* the Sherman-Morrison formula gives the inverse of a rank-n updated matrix
* 
* A' = A + R X T^t 
*
* as
*
*   A'^{-1} = A^{-1} -  A^{-1}R U^{-1} T^t A^{-1}
*
* U = X^{-1} + T^t A^{-1} R
*
* 
* In our case the U matrix become
* U_ij =  S - T^t_i A^{-1} R_j = 
*          x_i^t J_0^{-1} J x_j 
*
*     First we update (full) U Matrix, where we can skip the contributions to
*     lines and column that have been expelled from the subspace
*     We will later extract the part of the U matrix within the non-singular
*     space.
C     a) update all NDOVEC columns j
      IF (NVEC_SUB.GT.0) THEN
        CALL SKPVCD(LUJVEC,NVEC-NDOVEC,VEC1,1,LBLK)
        DO J = NVEC-NDOVEC+1,NVEC
*. Diagonal times Jx_j
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
          CALL REWINO(LUJDIA)
          CALL REWINO(LUSCR)
          CALL DMTVCD(VEC1,VEC2,LUJDIA,LUJVEC,LUSCR,DAMP,0,1,LBLK)
          CALL REWINO(LUVEC)
          CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
c        DO I = 1, NVEC_SUB
          DO I = NVEC_SKIP+1, NVEC
*. x_i^t J0^{-1}Jx_j 
            CALL REWINO(LUSCR)
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
            UIJ = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
            IJDX = (J-1)*MAXVECMAX + I
            UMAT(IJDX) = UIJ
          END DO
        END DO
C     b) update all NDOVEC lines i
        CALL REWINO(LUVEC)
        CALL SKPVCD(LUVEC,NVEC-NDOVEC,VEC1,1,LBLK)
        DO I = NVEC-NDOVEC+1,NVEC
*. Diagonal times x_i
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
          CALL REWINO(LUJDIA)
          CALL REWINO(LUSCR)
          CALL DMTVCD(VEC1,VEC2,LUJDIA,LUVEC,LUSCR,DAMP,0,1,LBLK)
          CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
          DO J = NVEC_SKIP+1, NVEC-NDOVEC
*                                ^ do not doubly calc. lower right square
*. x_j^t J^T J0^{-1}x_i 
            CALL REWINO(LUSCR)
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
            UIJ = INPRDD(VEC1,VEC2,LUJVEC,LUSCR,0,LBLK)
            IJDX = (J-1)*MAXVECMAX + I
            UMAT(IJDX) = UIJ
          END DO
        END DO
      
      END IF ! NVEC_SUB.GT.0

      IF (NVEC_SUB.GT.0) THEN

      FACT_U=.FALSE.
c      IF (NVEC_SUB.EQ.0) FACT_U=.TRUE. ! skip section in this case
      DO WHILE(.NOT.FACT_U)
* pack U matrix in NVEC_SUB x NVEC_SUB array:
      DO J = 1, NVEC_SUB
        DO I = 1, NVEC_SUB
          IJDX1 = (J-1)*NVEC_SUB+I
          IJDX2 = (J+NVEC_SKIP-1)*MAXVECMAX+I+NVEC_SKIP
          SCR(KLSU-1+IJDX1) = UMAT(IJDX2)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Try to invert U matrix subspace ',
     &       NVEC_SKIP+1,' to ',NVEC
c        WRITE(6,*) ' Full U matrix '
c        CALL WRTMAT(UMAT,NVEC,NVEC,MAXVECMAX,NVEC_SUB)
        WRITE(6,*) ' U matrix '
        CALL WRTMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*. Inverse U matrix and keep in KLU
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)

C well: this routine takes only the lower diagonal part:      
c      CALL INVERT_BY_DIAG(SCR(KLSU),SCR(KLEIGVEC),SCR(KLSCR),
c     &                    SCR(KLVEC),NVEC_SUB)
C version with symmetrized U
c      CALL SYMMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB)
c      CALL INVERT_BY_DIAG(SCR(KLSU),SCR(KLEIGVEC),SCR(KLSCR),
c     &                    SCR(KLVEC),NVEC_SUB)

c VESION WITH LINPACK ROUTINES:
c factorize U matrix and get its condition on COND
      CALL DGECO(SCR(KLSU),NVEC_SUB,NVEC_SUB,IPIV,COND,SCR(KLVEC))
      IF(COND.LT.1D-10) THEN
        IF(NTEST.GE.10) THEN
          WRITE(6,'(x,a,e8.2,a)')
     &         'DGECO gave near singularity (COND=',COND,')'
          WRITE(6,'(x,a)') 'I increase NVEC_SKIP and retry!'
        END IF
        NVEC_SKIP = NVEC_SKIP+1
        NVEC_SUB  = NVEC_SUB-1
        IF (NVEC_SKIP.EQ.NVEC) STOP 'I''m fucked!'
      ELSE
        FACT_U = .TRUE.
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) 'Factorised matrix:'
          CALL WRTMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
          WRITE(6,*) 'Pivot list:'
          CALL WRTIMAT(IPIV,NVEC_SUB,1,NVEC_SUB,1)
        END IF
      END IF
      END DO ! WHILE(.NOT.FACT_U)
c and get its inverse
      CALL DGEDI(SCR(KLSU),NVEC_SUB,NVEC_SUB,
     &           IPIV,DUMMY,SCR(KLSCR),1)

      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverted U matrix '
        CALL WRTMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*
* And then J0{-1} R U{-1} T^t J0^{-1} VEC
*
*. 1 : J0^{-1} VEC
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,DAMP,1,1,LBLK)
*. 2 : T^t J0^{-1} VEC 
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO I = 1, NVEC_SUB
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
        CALL REWINO(LUSCR)
        SCR(KLVEC-1+I) = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*. 3 :   VEC' = U{-1} T^t J0^{-1} VEC (Vec' is a vector in the subspace)
      CALL MATVCB(SCR(KLSU),SCR(KLVEC),SCR(KLVEC2),NVEC_SUB,
     &            NVEC_SUB,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' U(-1) T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC2),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*  4 : J0^{-1} R U{-1} T^t J0^{-1} VEC
*      = J0^{-1} ( sum_j Jx_j vec_j)  - sum_j x_j vec_j
*
*. sum_j Jx_j vec_j on LUSCR
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  LUSCR before call to MVCSMD'
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
      CALL REWINO(LUSCR)
      CALL REWINO(LUSCR2)
      CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUJVEC,SCR(KLVEC2),LUSCR,LUSCR2,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*.  J0^{-1}  sum_j Jx_j vec_j on LUSCR2
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUSCR,LUSCR2,DAMP,1,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' J0(-1) (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR2,1,LBLK)
      END IF
*.  sum_j x_j vec_j on LUSCR (LUOUT is used as scratch)
      CALL REWINO(LUSCR)
      CALL REWINO(LUOUT)
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUVEC,SCR(KLVEC2),LUSCR,LUOUT,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' xU(-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*. and the synthesis in LUOUT
C  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK
      ONE = 1.0D0
      ONEM = -1.0D0
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR2,LUSCR,LUOUT,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Total low rank correction'
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
c trust factor
      TFAC=100d0
c
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' Current STEPNORM_LAST, TFAC: ',STEPNORM_LAST,TFAC
      END IF
* TEST NORM
      XNORM = SQRT(INPRDD(VEC1,VEC2,LUOUT,LUOUT,1,LBLK))
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' Norm of low rank correction: ',XNORM
      END IF
      ELSE  ! IF NVEC_SUB==0
        XNORM = 0d0
      END IF
*. And the remainders 
*. J0^{-1} Vecin on LUSCR
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,DAMP,1,1,LBLK)
      YNORM = SQRT(INPRDD(VEC1,VEC2,LUSCR,LUSCR,1,LBLK))
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' Norm of perturbation correction: ',YNORM
      END IF
      IF (NVEC_SUB.GT.0) THEN
        XYPRD = INPRDD(VEC1,VEC2,LUSCR,LUOUT,1,LBLK)
        IF (NTEST.GE.10) THEN
          WRITE(6,*) ' Overlap perturbation/LR correction: ',XYPRD
        END IF
      END IF
*. (J0^{-1} -  J0^{-1} R U{-1} T^t J0^{-1} ) VECIN
c
      IF (XNORM.LT.XMXSTP_APRJ.AND.
     &    XNORM.LE.STEPNORM_LAST.AND.
     &    XNORM.LT.YNORM*TFAC) THEN
        IF (NVEC_SUB.GT.0) THEN
          CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR,LUOUT,LUSCR2,1,LBLK)
*. Well the result ended up on LUSCR2, place it properly on LUOUT
C      COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
          XNORM2=SQRT(INPRDD(VEC1,VEC2,LUSCR2,LUSCR2,1,LBLK))
          IF (XNORM2.LE.XMXSTP) THEN
            CALL COPVCD(LUSCR2,LUOUT,VEC1,1,LBLK)
            STEPNORM_LAST=SQRT(INPRDD(VEC1,VEC2,LUOUT,LUOUT,1,LBLK))
          ELSE
            WRITE (6,*) 'Low rank correction failed(2)!' 
            WRITE (6,*) 'Reverting to perturbation step!'
            IF (YNORM.LT.XMXSTP_APRJ) THEN
              WRITE (6,*) 'Restarting subspace with current vector'
              NVEC_SKIP=NVEC
              STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
            ELSE
              WRITE (6,*) 'Deleting current subspace'
              NVEC_SKIP=NVEC+1
              STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
            END IF
            CALL COPVCD(LUSCR,LUOUT,VEC1,1,LBLK)
          END IF
        ELSE
          CALL COPVCD(LUSCR,LUOUT,VEC1,1,LBLK)
        END IF
      ELSE
        WRITE (6,*) 'Low rank correction failed!' 
        WRITE (6,*) 'Reverting to perturbation step!'
        IF (YNORM.LT.XMXSTP_APRJ) THEN
          WRITE (6,*) 'Restarting subspace with current vector'
          NVEC_SKIP=NVEC
          STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
        ELSE
          WRITE (6,*) 'Deleting current subspace'
          NVEC_SKIP=NVEC+1
          STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
        END IF
        CALL COPVCD(LUSCR,LUOUT,VEC1,1,LBLK)
      END IF
*
 1002 CONTINUE
      IF(NVEC_SUB.EQ.0) THEN
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUOUT,DAMP,1,1,LBLK)
      END IF
*
      IF (NTEST.GE.10) THEN
        XNORM=SQRT(INPRDD(VEC1,VEC2,LUOUT,LUOUT,1,LBLK))
        WRITE(6,*) ' Norm of new correction vector: ',XNORM
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector '
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUIN,1,LBLK)
        WRITE(6,*) ' Inverted approximate Jacobian times Vector '
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE APRJAC_TV3(NVEC,NNEWVEC,NDELVEC,LUIN,LUOUT,LUVEC,
     &           LUJVEC,LUJDIA,VEC1,VEC2,
     &           SMAT,UMAT,SCR,
     &           N_CC_AMP,LUSCR,LUSCR2,MAXVEC,MAXVECMAX,
     &           XMXSTP,XMXSTP_APRJ)
*
* An approximate Jacobian is given in the form of 
* a diagonal approximation, JACDIA, and  
* NVEC vectors(in LUVEC)  and Jacobian times these vectors(LUJVEC).
*
* Find Inverse approximate Jacobian times the vector in luin to
* obtain vector in luout.
* Largest number vectors in subspave of approximate in Jacobian is MAXVEC
*
* Jeppe Olsen, March 2000
*
* Modified version: reduced I/O to ord(3*NVEC)
*  Andreas, Jan 2004
*
      INCLUDE 'implicit.inc'
*. Local vectors
      DIMENSION VEC1(*),VEC2(*)
*. Scratch space : Should atleast be length :
      DIMENSION UMAT(*),SMAT(*),SCR(*)
*. Local mini-scratch-space:
      DIMENSION IPIV(NVEC)
*. A SAVE variable to remember skipped vectors:
      INTEGER, SAVE :: NVEC_SKIP=0
      REAL*8, SAVE :: DAMP_OLD=-999.999d0, STEPNORM_LAST=1000d0
*.
      LOGICAL FACT_U, TOO_LARGE_STEP
      REAL*8 INPRDD
*
      NTEST = 150
      LBLK = -1
*. 
*. Obtain overlap of input vectors
*
      KLFREE = 1
*
      KLSU = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLEIGVEC = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLSCR = KLFREE
      KLFREE = KLFREE + NVEC**2
*
      KLVEC = KLFREE
      KLFREE = KLFREE + NVEC 
*
      KLVEC2 = KLFREE
      KLFREE = KLFREE + NVEC 
*
      STOP 'under construction'

      IF(NNEWVEC.LE.0) THEN
        WRITE(6,'(/,x,a,i2,/,x,a,/)')
     &   'ARPRJA_TV: The input for NNEWVEC (',NNEWVEC,') is definitely',
     &   '           silly! It is actually hard to continue....'
        STOP 'NNEWVEC is silly'
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'APRJA_TV3 reporting to work '
        WRITE(6,*) '==========================='
        WRITE(6,*) ' LUIN,  LUOUT         ', LUIN, LUOUT
        WRITE(6,*) ' LUVEC, LUJVEC, LUJDIA ',LUVEC,LUJVEC,LUJDIA
        WRITE(6,*) ' LUSCR, LUSCR2         ',LUSCR,LUSCR2
        WRITE(6,*) ' NVEC, NNEWVEC, MAXVEC = ', NVEC,NNEWVEC,MAXVEC 
        WRITE(6,*) ' XMXSTP, XMXSTP_APRJ = ',XMXSTP,XMXSTP_APRJ
      END IF
      IF(NTEST.GE.150) THEN
        XNRM = SQRT(INPRDD(VEC1,VEC2,LUIN,LUIN,1,LBLK))
        WRITE(6,*) 'Contents of LUIN: norm = ',XNRM
        XNRM = SQRT(INPRDD(VEC1,VEC2,LUJDIA,LUJDIA,1,LBLK))
        WRITE(6,*) 'Contents of LUJDIA: norm = ',XNRM 
        WRITE(6,*) 'Contents of LUVEC:'
        CALL REWINO(LUVEC)
        DO II = 1, NVEC
          XNRM = SQRT(INPRDD(VEC1,VEC2,LUVEC,LUVEC,0,LBLK))
          WRITE(6,*) II,'. ',XNRM
        END DO  
        WRITE(6,*) 'Contents of LUJVEC:'
        CALL REWINO(LUJVEC)
        DO II = 1, NVEC
          XNRM = SQRT(INPRDD(VEC1,VEC2,LUJVEC,LUJVEC,0,LBLK))
          WRITE(6,*) II,'. ',XNRM
        END DO  
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vectors on file LUVEC '
        CALL REWINO(LUVEC)
        DO IVEC = 1, NVEC
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
          CALL WRTVCD(VEC1,LUVEC,0,LBLK)
        END DO
        WRITE(6,*) ' Vectors in file LUJVEC '
        CALL REWINO(LUJVEC)
        DO IVEC = 1, NVEC
          CALL WRTVCD(VEC1,LUJVEC,0,LBLK)
        END DO
      END IF
    
      ! we first probe the simple perturbation correction:
      DAMP = 0d0
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,DAMP,1,1,LBLK)
      XNORM = SQRT(INPRDD(VEC1,VEC2,LUSCR,LUSCR,1,LBLK))
      WRITE(6,*) 'NORM,DAMP : ',XNORM, DAMP

      ! normally only the new vectors' contribution to U is needed:
      NDOVEC = NNEWVEC 
      NVEC_SUB = NVEC
      IF (NVEC_SUB.GT.MAXVEC) THEN
        STOP 'OH NO!'
      END IF

      IF (NTEST.GE.100) WRITE(6,*) 'NDOVEC, NVECSUB: ',NDOVEC, NVECSUB

*
* We will use vectors NVEC_SKIP+1 to NVEC to span the 
* subspace used for the projection to the subspace which is   
*
* P = sum_{ij} x_i s_{ij}-1 x^t_j
*
* The approximate Jacobian is thus (J0 is diagonal part of Jacobian)
*
* J_apr = J0 (1-P) + Jex P
*       = J0 + (Jex-J0) P
*       = J0 + sum_ij (s_i - s0_i) S_{ij}-1 x^t_j (s_i = Jex x_i,s0_i = J0 x_i)
*       = J0  + R S^{-1} T^t
*
* 
* where R is the collection of NVEC_SUB column vectors (s_i - s0_i)
* and T is the collection of NVEC_SUB column vectors x_i.
*
* the Sherman-Morrison formula gives the inverse of a rank-n updated matrix
* 
* A' = A + R X T^t 
*
* as
*
*   A'^{-1} = A^{-1} -  A^{-1}R U^{-1} T^t A^{-1}
*
* U = X^{-1} + T^t A^{-1} R
*
* 
* In our case the U matrix become
* U_ij =  S - T^t_i A^{-1} R_j = 
*          x_i^t J_0^{-1} J x_j 
*
*     First we update (full) U Matrix, where we can skip the contributions to
*     lines and column that have been expelled from the subspace
*     We will later extract the part of the U matrix within the non-singular
*     space.
C     a) update all NDOVEC columns j
      IF (NVEC_SUB.GT.0) THEN
        CALL SKPVCD(LUJVEC,NVEC-NDOVEC,VEC1,1,LBLK)
        DO J = NVEC-NDOVEC+1,NVEC
*. Diagonal times Jx_j
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
          CALL REWINO(LUJDIA)
          CALL REWINO(LUSCR)
          CALL DMTVCD(VEC1,VEC2,LUJDIA,LUJVEC,LUSCR,DAMP,0,1,LBLK)
          CALL REWINO(LUVEC)
          CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
c        DO I = 1, NVEC_SUB
          DO I = NVEC_SKIP+1, NVEC
*. x_i^t J0^{-1}Jx_j 
            CALL REWINO(LUSCR)
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
            UIJ = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
            IJDX = (J-1)*MAXVECMAX + I
            UMAT(IJDX) = UIJ
          END DO
        END DO
C     b) update all NDOVEC lines i
        CALL REWINO(LUVEC)
        CALL SKPVCD(LUVEC,NVEC-NDOVEC,VEC1,1,LBLK)
        DO I = NVEC-NDOVEC+1,NVEC
*. Diagonal times x_i
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
          CALL REWINO(LUJDIA)
          CALL REWINO(LUSCR)
          CALL DMTVCD(VEC1,VEC2,LUJDIA,LUVEC,LUSCR,DAMP,0,1,LBLK)
          CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
          DO J = NVEC_SKIP+1, NVEC-NDOVEC
*                                ^ do not doubly calc. lower right square
*. x_j^t J^T J0^{-1}x_i 
            CALL REWINO(LUSCR)
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
            UIJ = INPRDD(VEC1,VEC2,LUJVEC,LUSCR,0,LBLK)
            IJDX = (J-1)*MAXVECMAX + I
            UMAT(IJDX) = UIJ
          END DO
        END DO
      
      END IF ! NVEC_SUB.GT.0

      IF (NVEC_SUB.GT.0) THEN

      FACT_U=.FALSE.
c      IF (NVEC_SUB.EQ.0) FACT_U=.TRUE. ! skip section in this case
      DO WHILE(.NOT.FACT_U)
* pack U matrix in NVEC_SUB x NVEC_SUB array:
      DO J = 1, NVEC_SUB
        DO I = 1, NVEC_SUB
          IJDX1 = (J-1)*NVEC_SUB+I
          IJDX2 = (J+NVEC_SKIP-1)*MAXVECMAX+I+NVEC_SKIP
          SCR(KLSU-1+IJDX1) = UMAT(IJDX2)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Try to invert U matrix subspace ',
     &       NVEC_SKIP+1,' to ',NVEC
c        WRITE(6,*) ' Full U matrix '
c        CALL WRTMAT(UMAT,NVEC,NVEC,MAXVECMAX,NVEC_SUB)
        WRITE(6,*) ' U matrix '
        CALL WRTMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*. Inverse U matrix and keep in KLU
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)

C well: this routine takes only the lower diagonal part:      
c      CALL INVERT_BY_DIAG(SCR(KLSU),SCR(KLEIGVEC),SCR(KLSCR),
c     &                    SCR(KLVEC),NVEC_SUB)
C version with symmetrized U
c      CALL SYMMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB)
c      CALL INVERT_BY_DIAG(SCR(KLSU),SCR(KLEIGVEC),SCR(KLSCR),
c     &                    SCR(KLVEC),NVEC_SUB)

c VESION WITH LINPACK ROUTINES:
c factorize U matrix and get its condition on COND
      CALL DGECO(SCR(KLSU),NVEC_SUB,NVEC_SUB,IPIV,COND,SCR(KLVEC))
      IF(COND.LT.1D-10) THEN
        IF(NTEST.GE.10) THEN
          WRITE(6,'(x,a,e8.2,a)')
     &         'DGECO gave near singularity (COND=',COND,')'
          WRITE(6,'(x,a)') 'I increase NVEC_SKIP and retry!'
        END IF
        NVEC_SKIP = NVEC_SKIP+1
        NVEC_SUB  = NVEC_SUB-1
        IF (NVEC_SKIP.EQ.NVEC) STOP 'I''m fucked!'
      ELSE
        FACT_U = .TRUE.
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) 'Factorised matrix:'
          CALL WRTMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
          WRITE(6,*) 'Pivot list:'
          CALL WRTIMAT(IPIV,NVEC_SUB,1,NVEC_SUB,1)
        END IF
      END IF
      END DO ! WHILE(.NOT.FACT_U)
c and get its inverse
      CALL DGEDI(SCR(KLSU),NVEC_SUB,NVEC_SUB,
     &           IPIV,DUMMY,SCR(KLSCR),1)

      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverted U matrix '
        CALL WRTMAT(SCR(KLSU),NVEC_SUB,NVEC_SUB,NVEC_SUB,NVEC_SUB)
      END IF
*
* And then J0{-1} R U{-1} T^t J0^{-1} VEC
*
*. 1 : J0^{-1} VEC
C            DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,DAMP,1,1,LBLK)
*. 2 : T^t J0^{-1} VEC 
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      DO I = 1, NVEC_SUB
C               INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
        CALL REWINO(LUSCR)
        SCR(KLVEC-1+I) = INPRDD(VEC1,VEC2,LUVEC,LUSCR,0,LBLK)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*. 3 :   VEC' = U{-1} T^t J0^{-1} VEC (Vec' is a vector in the subspace)
      CALL MATVCB(SCR(KLSU),SCR(KLVEC),SCR(KLVEC2),NVEC_SUB,
     &            NVEC_SUB,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' U(-1) T(t) J0(-1) Vecin '
        CALL WRTMAT(SCR(KLVEC2),1,NVEC_SUB,1,NVEC_SUB)
      END IF
*  4 : J0^{-1} R U{-1} T^t J0^{-1} VEC
*      = J0^{-1} ( sum_j Jx_j vec_j)  - sum_j x_j vec_j
*
*. sum_j Jx_j vec_j on LUSCR
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  LUSCR before call to MVCSMD'
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
      CALL REWINO(LUSCR)
      CALL REWINO(LUSCR2)
      CALL SKPVCD(LUJVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUJVEC,SCR(KLVEC2),LUSCR,LUSCR2,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*.  J0^{-1}  sum_j Jx_j vec_j on LUSCR2
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUSCR,LUSCR2,DAMP,1,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' J0(-1) (Jx)(U-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR2,1,LBLK)
      END IF
*.  sum_j x_j vec_j on LUSCR (LUOUT is used as scratch)
      CALL REWINO(LUSCR)
      CALL REWINO(LUOUT)
      CALL SKPVCD(LUVEC,NVEC_SKIP,VEC1,1,LBLK)
      CALL MVCSMD(LUVEC,SCR(KLVEC2),LUSCR,LUOUT,VEC1,VEC2,
     &             NVEC_SUB,0,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' xU(-1)x(T)J0(-1) Vecin '
        CALL WRTVCD(VEC1,LUSCR,1,LBLK)
      END IF
*. and the synthesis in LUOUT
C  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK
      ONE = 1.0D0
      ONEM = -1.0D0
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR2,LUSCR,LUOUT,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Total low rank correction'
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
c trust factor
      TFAC=100d0
c
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' Current STEPNORM_LAST, TFAC: ',STEPNORM_LAST,TFAC
      END IF
* TEST NORM
      XNORM = SQRT(INPRDD(VEC1,VEC2,LUOUT,LUOUT,1,LBLK))
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' Norm of low rank correction: ',XNORM
      END IF
      ELSE  ! IF NVEC_SUB==0
        XNORM = 0d0
      END IF
*. And the remainders 
*. J0^{-1} Vecin on LUSCR
      CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUSCR,DAMP,1,1,LBLK)
      YNORM = SQRT(INPRDD(VEC1,VEC2,LUSCR,LUSCR,1,LBLK))
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' Norm of perturbation correction: ',YNORM
      END IF
      IF (NVEC_SUB.GT.0) THEN
        XYPRD = INPRDD(VEC1,VEC2,LUSCR,LUOUT,1,LBLK)
        IF (NTEST.GE.10) THEN
          WRITE(6,*) ' Overlap perturbation/LR correction: ',XYPRD
        END IF
      END IF
*. (J0^{-1} -  J0^{-1} R U{-1} T^t J0^{-1} ) VECIN
c
      IF (XNORM.LT.XMXSTP_APRJ.AND.
     &    XNORM.LE.STEPNORM_LAST.AND.
     &    XNORM.LT.YNORM*TFAC) THEN
        IF (NVEC_SUB.GT.0) THEN
          CALL VECSMD(VEC1,VEC2,ONE,ONEM,LUSCR,LUOUT,LUSCR2,1,LBLK)
*. Well the result ended up on LUSCR2, place it properly on LUOUT
C      COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
          XNORM2=SQRT(INPRDD(VEC1,VEC2,LUSCR2,LUSCR2,1,LBLK))
          IF (XNORM2.LE.XMXSTP) THEN
            CALL COPVCD(LUSCR2,LUOUT,VEC1,1,LBLK)
            STEPNORM_LAST=SQRT(INPRDD(VEC1,VEC2,LUOUT,LUOUT,1,LBLK))
          ELSE
            WRITE (6,*) 'Low rank correction failed(2)!' 
            WRITE (6,*) 'Reverting to perturbation step!'
            IF (YNORM.LT.XMXSTP_APRJ) THEN
              WRITE (6,*) 'Restarting subspace with current vector'
              NVEC_SKIP=NVEC
              STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
            ELSE
              WRITE (6,*) 'Deleting current subspace'
              NVEC_SKIP=NVEC+1
              STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
            END IF
            CALL COPVCD(LUSCR,LUOUT,VEC1,1,LBLK)
          END IF
        ELSE
          CALL COPVCD(LUSCR,LUOUT,VEC1,1,LBLK)
        END IF
      ELSE
        WRITE (6,*) 'Low rank correction failed!' 
        WRITE (6,*) 'Reverting to perturbation step!'
        IF (YNORM.LT.XMXSTP_APRJ) THEN
          WRITE (6,*) 'Restarting subspace with current vector'
          NVEC_SKIP=NVEC
          STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
        ELSE
          WRITE (6,*) 'Deleting current subspace'
          NVEC_SKIP=NVEC+1
          STEPNORM_LAST=MIN(TFAC*YNORM,XMXSTP_APRJ)
        END IF
        CALL COPVCD(LUSCR,LUOUT,VEC1,1,LBLK)
      END IF
*
 1002 CONTINUE
      IF(NVEC_SUB.EQ.0) THEN
        CALL DMTVCD(VEC1,VEC2,LUJDIA,LUIN,LUOUT,DAMP,1,1,LBLK)
      END IF
*
      IF (NTEST.GE.10) THEN
        XNORM=SQRT(INPRDD(VEC1,VEC2,LUOUT,LUOUT,1,LBLK))
        WRITE(6,*) ' Norm of new correction vector: ',XNORM
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector '
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUIN,1,LBLK)
        WRITE(6,*) ' Inverted approximate Jacobian times Vector '
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE COP_SYMMAT(AIN,AOUT,I1OUT,NOUT)
*
* Extract symmetric packed submatrix AOUT from symmetric 
* packed submatrix AIN.
*
*. First row/column in output matrix is I1OUT, and
*. number of rows/columns in AOUT is NOUT
*
* symmetric matrix is packed in conventional form  with 
* row wise packing 
*
* Jeppe Olsen, March 4 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION AIN(*)
*. Output
      DIMENSION AOUT(*)
*
      IADD = I1OUT-1
      DO I = 1, NOUT 
        DO J = 1, I    
          IJ_OUT = I*(I-1)/2 + J
          IJ_IN  = (I+IADD)*(I+IADD-1)/2 + J + IADD
          AOUT(IJ_OUT) = AIN(IJ_IN)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output matrix from COP_SYMMAT'
        CALL PRSYM(AOUT,NOUT)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_JSELSUB
     &                    (NPOP,IPOP,LUDIA,LUJX,LUXJ,LUX,N_CC_AMP,TAMP,
     &                    VEC1,VEC2,VEC3,CCVEC1,CCVEC2,CCVEC3,UINV,SCR)
*
* Generate vectors and matrices for obtaining inverse of CC Jacobian    
* that is choosen to be exact in a subspace 
*
* The Approximate Jacobian is 
*                      ( 0    1       )
* J_apr = J_d + (J'X X)(              ) (J'^tX  X)^t
*                      ( 1  -J^{sub}  )
*
* where J' = J - J_d
* The Inverse of the approximate Jacobian is then  
*                           
* J_apr^{-1}  = J_d^{-1} - 
*                                         
*   (J_d^{-1}J'X  J_d^{-1}X ) U^{-|} (J_d^{-1}J'^tX J_d^{-1} X)^t
*
* Where        (J_{sub}   1 )
*       U =    (            ) + (J'^t X X)^t (J_d^{-1}J'X J_d^{-1} X)
*              (   1      0 )
*
* Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input
      INTEGER IPOP(NPOP)
      DIMENSION TAMP(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*),VEC3(*)
      DIMENSION CCVEC1(*), CCVEC2(*),CCVEC3(*) 
*. Output 
      DIMENSION UINV(2*NPOP,2*NPOP)
*. Local scratch
      DIMENSION SCR(2*NPOP,2*NPOP)
* 
      LBLK = -1
      NTEST = 00
*    
*
* Initialize U as 
*                   ( 0    1 )
*                   (        )
*                   ( 1    0 )
*
      ZERO = 0.0D0
      CALL SETVEC(UINV,ZERO,(2*NPOP) ** 2)
      DO IVEC = 1, NPOP
        UINV(IVEC,NPOP+IVEC) = 1.0D0
        UINV(NPOP+IVEC,IVEC) = 1.0D0
      END DO
*
* 1 : Create (J_d^{-1}J'X = (J_d^{-1}J -1 ) X and save on LUJX
* 2 : J_d^{-1}X on LUX
* 3 : J'_{sub}
*
      CALL REWINO(LUJX)
      CALL REWINO(LUX)
      DO IVEC = 1, NPOP
        ZERO = 0.0D0
        CALL SETVEC(CCVEC1,ZERO,N_CC_AMP)   
        CCVEC1(IPOP(IVEC)) = 1.0D0
* Jac x vector
C            JAC_T_VEC(L_OR_R,CC_AMP,JAC_VEC,TVEC,VEC1,VEC2,CCVEC)
        CALL JAC_T_VECF(2,CCVEC1,CCVEC2,TAMP,VEC1,VEC2,CCVEC3)
* J_d
C            VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
        LBLK = -1
        CALL VEC_FROM_DISC(CCVEC3,N_CC_AMP,1,LBLK,LUDIA)
*. J_{sub}_{ij} = x_i^t J' x_j in upper part of U
        DO JVEC = 1, NPOP
          UINV(JVEC,IVEC) = CCVEC2(IPOP(JVEC))
        END DO
        UINV(IVEC,IVEC) = UINV(IVEC,IVEC) - CCVEC3(IPOP(IVEC))

* J_d^{-1} Jac X
C            DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        CALL DIAVC2(CCVEC2,CCVEC2,CCVEC3,ZERO,N_CC_AMP)
        ONE = 1.0D0
        ONEM = -1.0D0
        CALL VECSUM(CCVEC2,CCVEC2,CCVEC1,ONE,ONEM,N_CC_AMP)
C            VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC1)
        CALL VEC_TO_DISC(CCVEC2,N_CC_AMP,0,LBLK,LUJX)
* J_d^{-1} X
        CALL DIAVC2(CCVEC1,CCVEC1,CCVEC3,ZERO,N_CC_AMP)
        CALL VEC_TO_DISC(CCVEC1,N_CC_AMP,0,LBLK,LUX)
      END DO
*
* 4 : J_d^{-1}J'^t X on LUXJ
* 5   (J'^tX  )^t (J_d^{-1} J'X J_d^{-1}X)
*
      CALL REWINO(LUXJ)
      DO IVEC = 1, NPOP
        ZERO = 0.0D0
        CALL SETVEC(CCVEC1,ZERO,N_CC_AMP)   
        CCVEC1(IPOP(IVEC)) = 1.0D0
* Jac^t x vector
C            JAC_T_VECF(L_OR_R,CC_AMP,JAC_VEC,TVEC,VEC1,VEC2,CCVEC)
        CALL JAC_T_VECF(1,CCVEC1,CCVEC2,TAMP,VEC1,VEC2,CCVEC3)
* Jac'^t x vector  in CCVEC2
        CALL VEC_FROM_DISC(CCVEC3,N_CC_AMP,1,LBLK,LUDIA)
        CCVEC2(IPOP(IVEC)) = CCVEC2(IPOP(IVEC)) - CCVEC3(IPOP(IVEC))
*. Update U matrix
        CALL REWINO(LUX)
        CALL REWINO(LUJX) 
        DO JVEC = 1, NPOP
          CALL VEC_FROM_DISC(CCVEC3,N_CC_AMP,0,LBLK,LUX)
          UINV(IVEC,JVEC+NPOP) =  UINV(IVEC,JVEC+NPOP) 
     &                      +  INPROD(CCVEC3,CCVEC2,N_CC_AMP)
          CALL VEC_FROM_DISC(CCVEC3,N_CC_AMP,0,LBLK,LUJX)
          UINV(IVEC,JVEC) =  UINV(IVEC,JVEC) 
     &                      +  INPROD(CCVEC3,CCVEC2,N_CC_AMP)
        END DO
          
* J_d^{-1}J'^t  X
C            VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
        CALL VEC_FROM_DISC(CCVEC3,N_CC_AMP,1,LBLK,LUDIA)
        CALL DIAVC2(CCVEC2,CCVEC2,CCVEC3,ZERO,N_CC_AMP)
C            VEC_TO_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUSC1)
        CALL VEC_TO_DISC(CCVEC2,N_CC_AMP,0,LBLK,LUXJ)
      END DO
*
* We are now just missing (    X )^t(J_d^{-1}J X J_d^{-1} X)
*
      DO IVEC = 1, NPOP
        CALL SETVEC(CCVEC1,ZERO,N_CC_AMP)
        CCVEC1(IPOP(IVEC)) = 1.0D0
        CALL REWINO(LUX)
        CALL REWINO(LUJX)
        DO JVEC = 1, NPOP
          CALL VEC_FROM_DISC(CCVEC2,N_CC_AMP,0,LBLK,LUX) 
          UINV(NPOP+IVEC,NPOP+JVEC) = UINV(NPOP+IVEC,NPOP+JVEC)
     &                           + INPROD(CCVEC1,CCVEC2,N_CC_AMP)
C            VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
          CALL VEC_FROM_DISC(CCVEC2,N_CC_AMP,0,LBLK,LUJX)
         
          UINV(NPOP+IVEC,JVEC) =  UINV(NPOP+IVEC,JVEC) 
     &                      +  INPROD(CCVEC1,CCVEC2,N_CC_AMP)
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The U matrix : '
        WRITE(6,*) ' ============== '
        CALL WRTMAT(UINV,2*NPOP,2*NPOP,2*NPOP,2*NPOP)
*
        WRITE(6,*) ' The LUX file containing J_d^{-1} X' 
        WRITE(6,*) ' =================================='
        CALL REWINO(LUX)
        DO IVEC = 1, NPOP
          CALL VEC_FROM_DISC(CCVEC1,N_CC_AMP,0,LBLK,LUX)
          CALL WRTMAT(CCVEC1,1,N_CC_AMP,1,N_CC_AMP)
        END DO
*
        WRITE(6,*) ' The LUX file containing J_d^{-1} JX' 
        WRITE(6,*) ' ====================================='
        CALL REWINO(LUJX)
        DO IVEC = 1, NPOP
          CALL VEC_FROM_DISC(CCVEC1,N_CC_AMP,0,LBLK,LUJX)
          CALL WRTMAT(CCVEC1,1,N_CC_AMP,1,N_CC_AMP)
        END DO
*
        WRITE(6,*) ' The LUX file containing J_d^{-1} J^tX' 
        WRITE(6,*) ' ======================================='
        CALL REWINO(LUXJ)
        DO IVEC = 1, NPOP
          CALL VEC_FROM_DISC(CCVEC1,N_CC_AMP,0,LBLK,LUXJ)
          CALL WRTMAT(CCVEC1,1,N_CC_AMP,1,N_CC_AMP)
        END DO
      END IF
*
* Invert U matrix 
      CALL INVMAT(UINV,SCR,2*NPOP,2*NPOP,ISING)
C          INVMAT(A,B,MATDIM,NDIM,ISING)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The inverted U matrix ' 
        CALL WRTMAT(UINV,2*NPOP,2*NPOP,2*NPOP,2*NPOP)
      END IF
*
      RETURN
      END 
      SUBROUTINE JSELSUB_TV(VECIN,VECOUT,LUX,LUJX,LUXJ,NVEC,
     &                      NVAR,LUDIA,UINV,VECSCR1,VECSCR2,
     &                      VECSCR3)
*
* Multiply a Vector in VECIN by Jacobian that is exact in a subspace 
* of NVEC vectors 
*
* The inverse Jacobian reads 
*                           
* J_apr^{-1}  = J_d^{-1} - 
*                                         
*   (J_d^{-1}J'X  J_d^{-1}X ) U^{-|} (J_d^{-1}J'^tX J_d^{-1} X)^t
*
*
* Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input 
      DIMENSION VECIN(*), UINV(2*NVEC,2*NVEC)
*. output
      DIMENSION VECOUT(*)
*. Scratch 
      DIMENSION VECSCR1(NVAR), VECSCR2(2*NVEC), VECSCR3(2*NVEC)
*
      ZERO = 0.0D0
      ONEM = -1.0D0
      LBLK = -1
*
*1 Diagonal^-1 * Vecin
*
      CALL VEC_FROM_DISC(VECSCR1,NVAR,1,LBLK,LUDIA) 
C          DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
      CALL DIAVC2(VECOUT,VECIN,VECSCR1,ZERO,NVAR)
*
* 2 : -(J_d^{-1}J'^tX J_d^{-1} X)^t VECIN
      CALL REWINO(LUXJ)
      CALL REWINO(LUX)
      DO IVEC = 1, NVEC
         CALL VEC_FROM_DISC(VECSCR1,NVAR,0,LBLK,LUXJ) 
         VECSCR2(IVEC) = INPROD(VECIN,VECSCR1,NVAR)
         CALL VEC_FROM_DISC(VECSCR1,NVAR,0,LBLK,LUX) 
         VECSCR2(IVEC+NVEC) = INPROD(VECIN,VECSCR1,NVAR)
      END DO
C?    WRITE(6,*) ' Overlap between LUXJ and LUX and input vector'
C?    CALL WRTMAT(VECSCR2,1,2*NVEC,1,2*NVEC)
*. U^{-1}  (J_d^{-1}J'^tX J_d^{-1} X)^t VECIN 
      
      CALL SCALVE(VECSCR2,ONEM,2*NVEC)
C MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(UINV,VECSCR2,VECSCR3,2*NVEC,2*NVEC,0)
C?    WRITE(6,*) ' Uinv times - overlap '
C?    CALL WRTMAT(VECSCR3,1,2*NVEC,1,2*NVEC)
*  (J_d^{-1}J'X  J_d^{-1}X ) VECSCR2
      CALL MVECSUM(VECSCR3,NVEC,NVAR,VECOUT,VECSCR1,LUJX,1,0)
      CALL MVECSUM(VECSCR3(1+NVEC),NVEC,NVAR,VECOUT,VECSCR1,LUX,1,0)
*
      RETURN
      END 
      SUBROUTINE MVECSUM(FAC,NVEC,NVAR,VECOUT,VECSCR,LU,IREW,INI_ZERO)
*
*
* A set of NVEC vectors reside on Disc. 
* Scale each vector with a factor and add
* 
* Version where three vectors may be in core
*
* IF INI_ZERO .ne. 0. VECOUT is initialized as the zero vector
*
* Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION FAC(NVEC)
*. Output 
      DIMENSION VECOUT(NVAR)
*. Scratch 
      DIMENSION VECSCR(NVAR)
*
      ZERO = 0.0D0
      ONE = 1.0D0
      LBLK = -1
*
      IF(INI_ZERO.NE.0) THEN
        ZERO = 0.0D0
        CALL SETVEC(VECOUT,ZERO,NVAR)
      END IF
*
      IF(IREW.NE.0) THEN
        CALL REWINO(LU)
      END IF
*
      DO IVEC = 1, NVEC
C        CALL VEC_FROM_DISC(VECSCR1,NVAR,0,LBLK,LUXJ) 
        CALL VEC_FROM_DISC(VECSCR,NVAR,0,LBLK,LU)
        CALL VECSUM(VECOUT,VECOUT,VECSCR,ONE,FAC(IVEC),NVAR)
      END DO
*
      RETURN
      END 
      SUBROUTINE TEST_DEC15
*
* Test new routines for exact Jacobian in subspace 
*
*.        
      INCLUDE 'implicit.inc'
      PARAMETER(NDIM = 100)
      DIMENSION UINV(2*NDIM,2*NDIM)
      DIMENSION SCR(2*NDIM,2*NDIM)
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
      DIMENSION VEC4(2*NDIM),VEC5(2*NDIM)
      DIMENSION XJACIN1(NDIM,NDIM),XJACIN2(NDIM,NDIM)
      INTEGER IPOP(NDIM)
*
      COMMON/JACCOM/XJAC(NDIM,NDIM)
*
      LBLK = -1
*. Initialize Jacobian 
      ZERO = 0.0D0
      CALL SETVEC(XJAC,ZERO,NDIM**2) 
*. Diagonal
      DO I = 1, NDIM
        XJAC(I,I) = I
      END DO
*. A bit of off-diagonal 
      DO I = 2, NDIM
        XJAC(I,1) = 0.1
        XJAC(1,I) = 0.2
      END DO
      DO I = 3, NDIM
        XJAC(I,2) = 0.2
        XJAC(2,I) = 0.4
      END DO
      DO I = 4, NDIM
        XJAC(I,3) = 0.3
        XJAC(3,I) = 0.6
      END DO
*. Diagonal on say LU45
      DO I = 1, NDIM
        VEC1(I) = XJAC(I,I) 
      END DO
      LUDIA = 45
      CALL VEC_TO_DISC(VEC1,NDIM,1,LBLK,LUDIA)
C     CALL VEC_TO_DISC(CCVEC2,N_CC_AMP,0,LBLK,LUJX)
*
      NPOP = 3
      IPOP(1) = 1
      IPOP(2) = 2
      IPOP(3) = 3
*
      LUJX = 46
      LUX  = 47
      LUXJ = 48
      XDUM = 0.0D0
      CALL GEN_JSELSUB(NPOP,IPOP,LUDIA,LUJX,LUXJ,LUX,NDIM,XDUM,
     &     XDUM,XDUM,XDUM,VEC1,VEC2,VEC3,UINV,SCR)
C     GEN_JSELSUB
C    &                    (NPOP,IPOP,LUDIA,LUJX,LUXJ,LUX,N_CC_AMP,TAMP,
C    &                    VEC1,VEC2,VEC3,CCVEC1,CCVEC2,CCVEC3,UINV,SCR)
*
* Direct inversion of Jacobian
C          INVMAT(A,B,MATDIM,NDIM,ISING)
      CALL COPVEC(XJAC,XJACIN1,NDIM**2)
      CALL INVMAT(XJACIN1,SCR,NDIM,NDIM,ISING)
      WRITE(6,*) ' Inverted Jacobian '
      CALL WRTMAT(XJACIN1,NDIM,NDIM,NDIM,NDIM)
*
* Invert by using constructed subspace approw
      DO IVEC = 1, NDIM
        CALL SETVEC(VEC1,0.0D0,NDIM)
        VEC1(IVEC) = 1.0D0
        CALL JSELSUB_TV(VEC1,VEC2,LUX,LUJX,LUXJ,NPOP,NDIM,LUDIA,
     &       UINV,VEC3,VEC4,VEC5)
        CALL COPVEC(VEC2,XJACIN2(1,IVEC),NDIM)
C     JSELSUB_TV(VECIN,VECOUT,LUX,LUJX,LUXJ,NVEC,
C    &                      NVAR,LUDIA,UINV,VECSCR1,VECSCR2,
C    &                      VECSCR3)
      END DO
*
      WRITE(6,*) ' Inverted JAcobian obtained from JSELSUB : '
      CALL WRTMAT(XJACIN2,NDIM,NDIM,NDIM,NDIM)
C  CMP2VC(VEC1,VEC2,NDIM,THRES)
      WRITE(6,*) ' Compare the two inverted Jacs '
      CALL CMP2VC(XJACIN1,XJACIN2,NDIM,1.0D-10)
*
      STOP
      END 
C            JAC_T_VECF(L_OR_R,CC_AMP,JAC_VEC,TVEC,VEC1,VEC2,CCVEC)
      SUBROUTINE JAC_T_VECF(L_OR_R,CC_AMP,XJAC_VEC,TVEC,VEC1,VEC2,CCVEC)
*
* Fusk routine for Jac times vector
*
      INCLUDE 'implicit.inc'
*
      PARAMETER(NDIM = 100)
      COMMON/JACCOM/XJAC(NDIM,NDIM)
*. Input vector 
      DIMENSION CC_AMP(NDIM)
*. Output vector 
      DIMENSION XJAC_VEC(NDIM)
*
      IF(L_OR_R.EQ.1) THEN
C       CALL MATVCB(UINV,VECSCR2,VECSCR3,2*NVEC,2*NVEC,0)
        CALL MATVCB(XJAC,CC_AMP,XJAC_VEC,NDIM,NDIM,1)
      ELSE 
        CALL MATVCB(XJAC,CC_AMP,XJAC_VEC,NDIM,NDIM,0)
      END IF
*
      NTEST = 00 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output from JAC_T_VECF fusk '
        CALL WRTMAT(CC_AMP,1,NDIM,1,NDIM)
        CALL WRTMAT(XJAC_VEC,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END 
*----------------------------------------------------------------------*
      subroutine update_dia(nfpos,luvec,luomg,ludia,luscr,luscr2,
     &                      vec1,vec2,n_cc_amp,imode)
*----------------------------------------------------------------------*
      implicit none
*----------------------------------------------------------------------*
* purpose: update the diagonal of the Jacobian via
*          imode = 1    BFGS
*          imode = 2    PSB
*          imode = 3    simple (Powell ?) update
*----------------------------------------------------------------------*
      integer, intent(in) ::
     &     imode
* input:   nfpos : position of latest delta T and delta Omega in files
*          luvec : unit with latest delta T
*          luomg : unit with latest delta Omega
*          ludia : latest diagonal approximation
      integer, intent(in) ::
     &     luvec, luomg, ludia, nfpos
* output:  ludia : updated diagonal approximation
*
* memory and file scratch:
      integer, intent(in)   ::
     &     luscr, luscr2, n_cc_amp
      real*8, intent(inout) ::
     &     vec1(*), vec2(*)

* local parameters
      integer, parameter ::
     &     ntest = 1000, lblk = -1
* local variables:
      integer ::
     &     icont, ncont
      real*8 ::
     &     denom(3), xnrm
* external function
      real*8 ::
     &     inprdd, inprod


      if (ntest.ge.100) then
        write(6,'(/x,a,/,x,a,/,x,a,3i4,/,x,a,i4,/)')
     &       'UPDATE_DIA messing around',
     &       '=========================',
     &       'LUVEC, LUOMG, LUDIA : ', luvec, luomg, ludia,
     &       'NFPOS               : ', nfpos
        if (imode.eq.1) write(6,*) 'mode: BFGS'
        if (imode.eq.2) write(6,*) 'mode: PSB'
        if (imode.eq.3) write(6,*) 'mode: simple'
      end if

      if (imode.lt.1.or.imode.gt.3) then
        write (6,'(/,x,a,i3,/,x,a,/)')
     &       'imode has non-implemented value: ',imode,
     &       'That''s too much, I quit!'
        stop 'update_dia'
      end if

      if (ntest.ge.1000) then
        write(6,*) 'A^0 on entry:', n_cc_amp
        xnrm = inprdd(vec1,vec2,ludia,ludia,1,lblk)
        write(6,*) 'Norm : ', xnrm
        call vec_from_disc(vec1,n_cc_amp,1,lblk,ludia)
        xnrm = inprod(vec1,vec1,n_cc_amp)
        write(6,*) 'Norm(2) : ', xnrm
c        call list_small(vec1,n_cc_amp,vec2,10)
        call wrtmat(vec2,1,10,1,10)
c        call wrtvcd(vec1,ludia,1,lblk)
        
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        call skpvcd(luomg,nfpos-1,vec1,1,lblk)
        xnrm = inprdd(vec1,vec2,luvec,luvec,0,lblk)
        write(6,*) 'Norm(delta T)   : ', xnrm
        xnrm = inprdd(vec1,vec2,luomg,luomg,0,lblk)
        write(6,*) 'Norm(delta Omg) : ', xnrm
      end if

* store old A^0 as first vector on luscr
      call copvcd(ludia,luscr,vec1,1,lblk)
      denom(1) = 1d0
      
      if (imode.eq.1) then
* BFGS: 2 contributions
        ncont = 2
*----------------------------------------------------------------------*
* Number one:
*
*         (A^0_ii delta T_i)^2 
*    - -------------------------
*        delta T * A^0 delta T
*
*----------------------------------------------------------------------*
* calculate A^0,old delta t -> luscr      
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        call rewino(ludia)
        call rewino(luscr2)
c DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
        call dmtvcd(vec1,vec2,ludia,luvec,luscr2,0d0,0,0,lblk)

* get delta T * A^0,old delta t -> denom(2)
        call rewino(luscr2)
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        denom(2) = -1d0/inprdd(vec1,vec2,luvec,luscr2,0,lblk)

* square each element of A^0,old delta t on luscr
        call skpvcd(luscr,1,vec1,1,lblk)
        call rewino(luscr2)
        call dmtvcd(vec1,vec2,luscr2,luscr2,luscr,0d0,0,0,lblk)

*----------------------------------------------------------------------*
* Number two:
*
*           (delta Omega_i)^2
*      + -----------------------
*         delta Omega * delta T
*
*----------------------------------------------------------------------*
* square each element of delta Omega -> luscr
        call skpvcd(luomg,nfpos-1,vec1,1,lblk)
        call skpvcd(luscr,2,vec1,1,lblk)
        call dmtvcd(vec1,vec2,luomg,luomg,luscr,0d0,0,0,lblk)

* get delta T * delta Omega -> denom2
        call skpvcd(luomg,nfpos-1,vec1,1,lblk)
        denom(3) = 1d0/inprdd(vec1,vec2,luomg,luomg,0,lblk)
      else if (imode.eq.2) then
        stop 'PSB not here'
      else if (imode.eq.3) then
        ncont = 2
* delta T * delta T -> -denom(2) = denom(3)
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        denom(2) = -1d0/inprdd(vec1,vec2,luvec,luvec,0,lblk)
        denom(3) = -1d0*denom(2)
* A^0 delta T -> luscr2
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        call rewino(ludia)
        call rewino(luscr2)
        call dmtvcd(vec1,vec2,luvec,ludia,luscr2,0d0,0,0,lblk)

* delta T A^0 delta T ->  luscr
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        call rewino(luscr2)
        call skpvcd(luscr,1,vec1,1,lblk)
        call dmtvcd(vec1,vec2,luvec,luscr2,luscr,0d0,0,0,lblk)

* delta Omg delta T -> luscr
        call skpvcd(luvec,nfpos-1,vec1,1,lblk)
        call skpvcd(luomg,nfpos-1,vec1,1,lblk)
        call skpvcd(luscr,2,vec1,1,lblk)
        call dmtvcd(vec1,vec2,luvec,luomg,luscr,0d0,0,0,lblk)

      end if
* add contents on luscr times denom array to ludia
c MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)      
      if (ntest.ge.100) then
        write(*,*) 'ncont = ', ncont
        write(*,*) 'denom = ', denom(1:ncont+1)
        call rewino(luscr)
        do icont = 1, ncont+1
          xnrm = inprdd(vec1,vec2,luscr,luscr,0,lblk)
          write(6,*) 'Norm : ', icont, xnrm
        end do
      end if
      call rewino(ludia)
      call rewino(luscr)
      call mvcsmd(luscr,denom,ludia,luscr2,vec1,vec2,1+ncont,1,lblk)

      if (ntest.ge.1000) then
        write(6,*) 'A^0 on exit:'
        xnrm = inprdd(vec1,vec2,ludia,ludia,1,lblk)
        write(6,*) 'Norm : ', xnrm
        call vec_from_disc(vec1,n_cc_amp,1,lblk,ludia)
c        call list_small(vec1,n_cc_amp,vec2,10)
        call wrtmat(vec2,1,10,1,10)
c        call wrtvcd(vec1,ludia,1,lblk)
      end if

      return
      end
c $Id$
