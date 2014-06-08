      SUBROUTINE DIIS_SIMPLE(LUEVEC,NVEC,NDIM,C)
*
* NVEC error-vectors are given on LUEVEC
* Use DIIS to find optimal combination of these vectors
*
*. Simple and generic DIIS code - more simple than efficient 
*. All subspace matrices constructed without reuse
*
*. Jeppe Olsen, August 2005
*
      INCLUDE 'wrkspc.inc'
*. Output : Optimial combination of directions 
       REAL*8 INPROD
       DIMENSION C(NVEC)
*. Non-redundant formalism : X = X(NVEC) + sum(i=1,NVEC-1)C(I)(X(I)-X(N))
*  corresponding error       E = E(NVEC) + sum(i=1,NVEC-1)C(I)(E(I)-E(N))
*
* Minimizing norm of gradient corresponds to solving linear set of equations
*
* A C = RHS
* A_ij = <E(I)-E(N)!E(J)-E(N)>, RHS_i = -<E(I)-E(N)!E(N)>
*
*
      NTEST = 100
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Entering DIIS '
        WRITE(6,*) ' Number of vectors ', NVEC
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DIISS')
*. A few scratch vectors
      LEN = MAX(NDIM,NVEC-1)
      CALL MEMMAN(KLVEC1,LEN,'ADDL ',2,'VC1_DI')
      CALL MEMMAN(KLVEC2,NDIM,'ADDL ',2,'VC2_DI')
      CALL MEMMAN(KLVEC3,NDIM,'ADDL ',2,'VC3_DI')
*. Space for A and RHS and for inverting 
      CALL MEMMAN(KLA,(NVEC-1)**2,'ADDL ',2,'AMAT  ')
      CALL MEMMAN(KLRHS,(NVEC-1), 'ADDL ',2,'RHS   ')
      CALL MEMMAN(KLSCRM1,(NVEC-1)**2,'ADDL ',2,'SCRM1 ')
      CALL MEMMAN(KLSCRM2,(NVEC-1)**2,'ADDL ',2,'SCRM2 ')
  
*. Read E(NVEC) in VEC1
      CALL REWINO(LUEVEC)
      DO IVEC = 1, NVEC
C?      WRITE(6,*) ' Initial loop, IVEC = ', IVEC
C  VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
        CALL VEC_FROM_DISC(WORK(KLVEC1),NDIM,0,-1,LUEVEC)
      END DO
*
      ONE = 1.0D0
      ONEM = -1.0D0
      DO IVEC = 1, NVEC-1
C?     WRITE(6,*) ' IVEC = ', IVEC
*. Read vector E(IVEC) in VEC2
       CALL REWINO(LUEVEC)
       DO KVEC = 1, IVEC
         CALL VEC_FROM_DISC(WORK(KLVEC2),NDIM,0,-1,LUEVEC)
       END DO
*. E(IVEC)-E(NVEC) in VEC2
       CALL VECSUM(WORK(KLVEC2),WORK(KLVEC2),WORK(KLVEC1),
     &             ONE,ONEM,NDIM)
       WORK(KLRHS-1+IVEC) = INPROD(WORK(KLVEC2),WORK(KLVEC1),NDIM)
       EEII = INPROD(WORK(KLVEC2),WORK(KLVEC2),NDIM)
       WORK(KLA-1+(IVEC-1)*(NVEC-1)+IVEC) = EEII
       DO JVEC = IVEC+1, NVEC-1
C?     WRITE(6,*) ' JVEC = ', JVEC
*. Read vector E(JVEC) on VEC3
         CALL VEC_FROM_DISC(WORK(KLVEC3),NDIM,0,-1,LUEVEC)
*. E(JVEC) - E(NVEC) in VEC3
         CALL VECSUM(WORK(KLVEC3),WORK(KLVEC3),WORK(KLVEC1),
     &               ONE,ONEM,NDIM)
         EEIJ = INPROD(WORK(KLVEC2),WORK(KLVEC3),NDIM)
         WORK(KLA-1+(IVEC-1)*(NVEC-1)+JVEC) = EEIJ
         WORK(KLA-1+(JVEC-1)*(NVEC-1)+IVEC) = EEIJ
       END DO
      END DO
*. Well, the RHS should contain a - so
      CALL SCALVE(WORK(KLRHS),ONEM,NVEC-1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'The <E(I)-E(N)!E(J)-E(N)> matrix in DIIS '
        CALL WRTMAT(WORK(KLA),NVEC-1,NVEC-1,NVEC-1,NVEC-1)
        WRITE(6,*) 'The <E(I)-E(N)!E(N)> vector in DIIS '
        CALL WRTMAT(WORK(KLRHS),1,NVEC-1,1,NVEC-1)
      END IF
*. Solve linear equations 
*. Invert EE
      IF(NVEC.GE.1) THEN
      CALL INVERT_BY_DIAG2(WORK(KLA),WORK(KLSCRM1),WORK(KLSCRM2),
     &                    WORK(KLVEC1),NVEC-1)
      END IF
      CALL MATVCB(WORK(KLA),WORK(KLRHS),C,NVEC-1,NVEC-1,0)
*. IN C, we have coefficients for C(I),I.ne.NVEC, obtain C(NVEC)
      C(NVEC) = 1.0D0-ELSUM(C,NVEC-1)
*. Clean up : leave the file at end of vector NVEC 
      IF(NVEC.GE.2) THEN
         CALL VEC_FROM_DISC(WORK(KLVEC2),NDIM,0,-1,LUEVEC)
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' The C-coefficients obtained from DIIS '
        CALL WRTMAT(C,1,NVEC,1,NVEC)
      END IF
*
C     WRITE(6,*) ' Leaving DIIS '
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'DIISS')
      RETURN
      END


        
c $Id$
