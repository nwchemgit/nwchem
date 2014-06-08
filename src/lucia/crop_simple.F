      SUBROUTINE CROP(NVEC,NVECP,MXNVEC,NDIM,LUE,LUP,A,
     &                EOUT,POUT,SCR,LUSCR,IDIRDEL)
*
*. CROP version of DIIS
*
* LUEVEC contains a list of errorvectors, LUPVEC contains
* the corresponding parameter lists.
* 
* Perform a CROP norm minimization to obtained improved
* error vector(EOOUT) and corresponding parametervector(POUT)
*
*
* Notice: the routine does modify the files LUEVEC, LUPVEC
*
*. Jeppe Olsen, Aug 2009
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION A(MXNVEC*(MXNVEC+1)/2)
*. Output
      DIMENSION EOUT(NDIM),POUT(NDIM)
*. Length of SCR: 3*MXNVEC**2 + 3*MAX(NDIM,MXNVEC)
      DIMENSION SCR(*)
*
C?    WRITE(6,*) ' NVEC at start of CROP', NVEC
      KLAUP = 1
      KLFREE = KLAUP + MXNVEC*MXNVEC
*
      NDIM_E = MAX(NDIM,MXNVEC)
*
      KLVEC1 = KLFREE
      KLFREE = KLVEC1 + NDIM_E
*
      KLVEC2 = KLFREE
      KLFREE = KLFREE + NDIM_E
*
      KLVEC3 = KLFREE
      KLFREE = KLFREE + NDIM_E
*
      KLS1 = KLFREE
      KLFREE = KLFREE + MXNVEC**2
*
      KLS2 = KLFREE
      KLFREE = KLFREE + MXNVEC**2
      CALL CROP_INNER(NVEC,NVECP,MXNVEC,NDIM,LUE,LUP,IDIRDEL,
     &     POUT,EOUT,A,SCR(KLAUP),
     &     SCR(KLVEC1),SCR(KLVEC2),SCR(KLVEC3),SCR(KLS1),SCR(KLS2),
     &     LUSCR)
*
      RETURN
      END
      SUBROUTINE CROP_INNER
     &           (NVEC,NVECP,MXNVEC,NDIM,LUE,LUP,IDIRDEL,
     &            POUT,EOUT,
     &            A, AUP, VEC1,VEC2,VEC3,SMAT1,SMAT2,LUSCR)
*
*. CROP version of DIIS, inner routine
*
* LUEVEC contains a list of errorvectors, LUPVEC contains
* the corresponding parameter lists.
* 
* Perform a CROP norm minimization to obtained improved
* error vector(EOOUT) and corresponding parametervector(POUT)
*
*
* Notice: the routine does modify the files LUEVEC, LUPVEC
*         as well as NVEC, NVECP
*
*. Jeppe Olsen, Aug 2009
*
*. Parameters
* NVEC: Number of vectors currenty residing on LUEVEC,LUPVEC
* NVECP: Number of vectors for which the overlap matrix of 
*        errors A(IJ) = <E(I)!E(J)> already has been constructed
* MXNVEC: Allowed number of vectors in subspace: If  NVEC > MXNVEC
* the first NVEC-MXNVEC on LUEVEC,LUPVEC are eliminated
* NDIM: Dimension of parameterspace
* LUEVEC: File containing errorvectors, each as single record
* LUPVEC: File containing parametervectors, each as single record
* IDIRDEL: The matrix <E(I)-E(NVEC)!E(J)-E(NVEC)> is set up directly
* POUT   : Output parameter set
* EOUT   : Output error vector
* A: A(IJ) = <E(I)!E(J)> has been set up for the first NVECP 
* AUP: Space for unpacked matrix over subspace vectors
* VEC1,VEC2,VEC3: vectors must hold either subspace vectors
*                 or parametervectors
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
      DIMENSION A(MXNVEC*(MXNVEC+1)/2),AUP(MXNVEC**2)
*. Output
      DIMENSION EOUT(NDIM),POUT(NDIM)
*. Scratch
      DIMENSION VEC1(*),VEC2(*),VEC3(*)
      DIMENSION  SMAT1(MXNVEC*MXNVEC),SMAT2(MXNVEC*MXNVEC)
*
      NTEST = 5
      IF(NTEST.GE.10) THEN
        WRITE(6,*)' Entering CROP_SIMPLE_INNER'
        WRITE(6,*) ' LUE, LUP = ', LUE, LUP
        WRITE(6,*) ' NVEC, NVECP, MXNVEC = ', NVEC,NVECP,MXNVEC
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The error vectors entering CROP'
        REWIND(LUE)
        DO I = 1, NVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUE)
          WRITE(6,*) ' Vector:', I
          CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END DO
      END IF
*
* ----------------------------------------
* 1: Reduce number of vectors if required
* ----------------------------------------
*
      IF(NVEC.GT.MXNVEC) THEN
*. Remove the first NVEC-MXNVEC  vectors from files
*. First parametervectors
        REWIND(LUP)
*. a: Skip NVEC-MXNVEC vectors
        DO I = 1,  NVEC-MXNVEC
C VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUP)
        END DO
*. b: copy MXNVEC vectors to LUSCR
        REWIND(LUSCR)
        DO I = 1, MXNVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUP)
          CALL VEC_TO_DISC(VEC1,NDIM,0,-1,LUSCR)
        END DO
*. c: Copy the vectors back  
        REWIND(LUSCR)
        REWIND(LUP)
        DO I = 1, MXNVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUSCR)
          CALL VEC_TO_DISC(VEC1,NDIM,0,-1,LUP)
        END DO
*. and error vectors
        REWIND(LUE)
*. a: Skip NVEC-MXNVEC vectors
        DO I = 1,  NVEC-MXNVEC
C VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUE)
        END DO
*. b: copy MXNVEC vectors to LUSCR
        REWIND(LUSCR)
        DO I = 1, MXNVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUE)
          CALL VEC_TO_DISC(VEC1,NDIM,0,-1,LUSCR)
        END DO
*. c: Copy the vectors back  
        REWIND(LUSCR)
        REWIND(LUE)
        DO I = 1, MXNVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUSCR)
          CALL VEC_TO_DISC(VEC1,NDIM,0,-1,LUE)
        END DO
*. And remove the eliminated vectors from A
        IJ = 0
        DO I = 1, MXNVEC
          DO J = 1, I
            IJ = IJ + 1
            I_ORIG = I + (NVEC-MXNVEC)
            J_ORIG = J + (NVEC-MXNVEC)
            AUP(IJ) = A(I_ORIG*(I_ORIG-1)/2+J_ORIG)
          END DO
        END DO
        CALL COPVEC(AUP,A,MXNVEC*(MXNVEC+1)/2)
*. And update the dimensions
        NVEC = MXNVEC
        NVECP = NVECP-(NVEC-MXNVEC)
        IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Updated A-matrix'
         CALL PRSYM(A,NVEC)
        END IF
      END IF
*     ^ End if vectors should be reset
*
* -----------------------------------
*. 2: Construct the error matrices
* -----------------------------------
*
*. A(IJ) = <E(I)-E(N)!E(J)-E(N)
*. RHS(I) = <E(I)-E(N)!E(N)>
*
      IF(IDIRDEL.EQ.1) THEN
        IF(NTEST.GE.1000) WRITE(6,*) ' Start construction of A'
*. Direct construction of error matrix and RHS without any reuse 
* 1: Obtain E(N) = E(NVEC) in VEC1
        REWIND(LUE)
        DO I = 1, NVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUE)
        END DO
        IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Error vector NVEC '
           CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END IF
*
        DO J = 1, NVEC-1
          REWIND(LUE)
          DO K = 1, J
            CALL VEC_FROM_DISC(VEC2,NDIM,0,-1,LUE)
          END DO
          IF(NTEST.GE.100) WRITE(6,*) ' J-VECTOR read in for J= ', J
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Errorvector J for J = ', J
            CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
          END IF
 
          ONE = 1.0D0
          ONEM = -1.0D0
* E(J) - E(N) in Vec2
          CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
          IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' Error vector J - NVEC for J =  ', J
             CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
          END IF
*. Element JJ
          JJ = J*(J+1)/2
          A(JJ) = INPROD(VEC2,VEC2,NDIM)
*. RHS(J) in AUP(J)
          AUP(J) = INPROD(VEC2,VEC1,NDIM)
          DO I = J+1,NVEC-1
*. E(I) - E(N) in VEC3A
            IJ = I*(I-1)/2 + J
            CALL VEC_FROM_DISC(VEC3,NDIM,0,-1,LUE)
            IF(NTEST.GE.100) WRITE(6,*) ' I-vector read in for I= ', I
            CALL VECSUM(VEC3,VEC3,VEC1,ONE,ONEM,NDIM)
            IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' Error vector I - NVEC for I =  ', I
               CALL WRTMAT(VEC3,1,NDIM,1,NDIM)
            END IF
* A(IJ) =<E(I) - E(N)!E(J)-E(N)
            A(IJ) = INPROD(VEC2,VEC3,NDIM)
          END DO
        END DO
        CALL COPVEC(AUP,VEC1,NVEC-1)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The A-matrix: '
          CALL PRSYM(A,NVEC-1)
        END IF
*. Outpacked matrix in AUP
        CALL TRIPAK(AUP,A,2,NVEC-1,NVEC-1)
C  TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      ELSE
*. <E(I)!E(J)> is in A for I = I,J leq NVECP
*. Obtain <E(I)!E(J) for I,J > NVECP
        REWIND(LUE)
        DO I = NVECP+1, NVEC
*. Read EVEC(I) in VEC3
          REWIND(LUE)
          DO II = 1, I
           CALL VEC_FROM_DISC(VEC3,NDIM,0,-1,LUE)
          END DO
          REWIND(LUE)
          DO J = 1, I-1
            IJ = I*(I-1)/2 + J
            CALL VEC_FROM_DISC(VEC2,NDIM,0,-1,LUE)
            A(IJ) = INPROD(VEC2,VEC3,NDIM)
          END DO
          II = I*(I-1)/2 + I
          A(II) = INPROD(VEC3,VEC3,NDIM)
        END DO
*. Construct <E(I)-E(N)!E(J)-E(N)> in AUP in unpacked form,
        DO I = 1, NVEC-1
          DO J = 1,  I
            IJ=(J-1)*(NVEC-1) + I
            JI=(I-1)*(NVEC-1) + J
*
            IJ_SYM = I*(I-1)/2+J
            NI_SYM = NVEC*(NVEC-1)/2+I
            NJ_SYM = NVEC*(NVEC-1)/2+J
            NN_SYM = NVEC*(NVEC-1)/2+NVEC
            AUP(IJ) = A(IJ_SYM)+A(NN_SYM)-A(NI_SYM)-A(NJ_SYM)
            AUP(JI) = AUP(IJ)
          END DO
          NI_SYM = NVEC*(NVEC-1)/2+I
          NN_SYM = NVEC*(NVEC-1)/2+NVEC
          VEC1(I) = A(NI_SYM)-A(NJ_SYM)
        END DO
      END IF
*     ^ End if <E(I)-E(N)!E(J)-E(N)> should be constructed direct
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' The matrix <E(I)-E(N)!E(I)-E(N)>'
        CALL WRTMAT(AUP,NVEC-1,NVEC-1,NVEC-1,NVEC-1)
        WRITE(6,*) ' The vector <E(I)-E(N)!E(N)>'
        CALL WRTMAT(VEC1,1,NVEC-1,1,NVEC-1)
      END IF
*
* ------------------------------------------------------
* 3: Solve the equations identifying mixing coefficients
*    AUP C = -RHS
* ------------------------------------------------------
*
      ONEM = -1.0D0
      CALL SCALVE(VEC1,ONEM,NVEC)
*
      IF(NVEC.GT.1) THEN
        CALL INVERT_BY_DIAG2(AUP,SMAT1,SMAT2,VEC2,NVEC-1)
        CALL MATVCB(AUP,VEC1,VEC2,NVEC-1,NVEC-1,0)
        VEC2(NVEC) = 1.0D0 - ELSUM(VEC2,NVEC-1)
      ELSE
        VEC2(1) = 1.0D0
      END IF
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' The mixing coefficients from CROP'
        CALL WRTMAT(VEC2,1,NVEC,1,NVEC)
      END IF
*. Assemble in EOUT the modified error and overwrite vector NVEC
C (FAC,NVEC,NVAR,VECOUT,VECSCR,LU,IREW,INI_ZERO)
      CALL MVECSUM(VEC2,NVEC,NDIM,EOUT,VEC3,LUE,1,1)
*. and overwrite
      REWIND LUE
      DO I = 1, NVEC-1
        CALL VEC_FROM_DISC(VEC3,NDIM,0,-1,LUE)
      END DO
      CALL VEC_TO_DISC(EOUT,NDIM,0,-1,LUE)
*. and the parameter vector
      CALL MVECSUM(VEC2,NVEC,NDIM,POUT,VEC3,LUP,1,1)
*. and overwrite
      REWIND LUP
      DO I = 1, NVEC-1
        CALL VEC_FROM_DISC(VEC3,NDIM,0,-1,LUP)
      END DO
      CALL VEC_TO_DISC(POUT,NDIM,0,-1,LUP)
*. The  files LUP and LUE are left at the end of files
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The error vectors leaving CROP'
        REWIND(LUE)
        DO I = 1, NVEC
          CALL VEC_FROM_DISC(VEC1,NDIM,0,-1,LUE)
          WRITE(6,*) ' Vector:', I
          CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END DO
      END IF
      RETURN
      END
*
c $Id$
