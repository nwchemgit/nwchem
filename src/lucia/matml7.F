      SUBROUTINE MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(T) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(T) FOR ITRNSP = 2
C
C     C =  FACTORC*C + FACTORAB* A(T) * B(T) FOR ITRNSP = 3
*
* Warning ITRNSP = 3 should only be used for small matrices, 
* as this path involves notunit strides in inner loops. 
* As Lasse points out, it is better to calculate C(T) = BA
* and then transpose C
*
C... JEPPE OLSEN,
*
* ITRNSP = 3 added, march 2003
C
*. Notice : If the summation index has dimension zero nothing
*           is performed
      IMPLICIT REAL*8           (A-H,O-Z)
      include 'akbk.inc'
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
*
      COMMON/XXTEST/ISETVECOPS(10)
*
c     CALL LUCIAQENTER('MATML')
*
!      cpu_MATML7  = util_cpusec()
!      wall_MATML7 = util_wallsec()
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL ' 
      WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      WRITE(6,*) ' FACTORC, FACTORAB, ITRNSP : ', 
     &             FACTORC, FACTORAB, ITRNSP
      WRITE(6,*)
      WRITE(6,*) ' A, B and initial C MATRIX FROM MATML7 ' 
      WRITE(6,*)
      CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
      WRITE(6,*)
      CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
*. Statistics
cVOG      XNCALL = XNCALL + 1
cVOG      XLCROW = XLCROW + NCROW
cVOG      XLCCOL = XLCCOL + NCCOL
cVOG      XLCROWCOL = XLCROWCOL + DBLE(NCROW)*DBLE(NCCOL)
cVOG      T_INI = 0.01*LTIMEX()
*
cVOG      IF(ITRNSP.EQ.1) THEN
cVOG        XNFLOP = XNFLOP + 2*DBLE(NCROW)*DBLE(NCCOL)*DBLE(NAROW)
cVOG      ELSE
cVOG        XNFLOP = XNFLOP + 2*DBLE(NCROW)*DBLE(NCCOL)*DBLE(NACOL)
cVOG      END IF
*. Hvilken vej skal jeg vælge ??
cVOG      IESSL = 0
cVOG      ICONVEX = 0
!*
!      IF(NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL .EQ. 0 .OR.
!     &   FACTORC.EQ.0.0D0) THEN
!        IZERO = 1
!      ELSE
!        IZERO = 0
!      END IF
!*
!      IF(IZERO.EQ.1.AND.NCROW*NCCOL.NE.0) THEN 
!        IF(FACTORC.NE.0.0D0) THEN
!         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
!        ELSE IF(.NOT.(IESSL.EQ.0.AND.ITRNSP.EQ.2) ) THEN
!         IF(ITRNSP.NE.1) THEN
!           ISETVECOPS(1) = ISETVECOPS(1) + NCROW*NCCOL
!           ZERO = 0.0D0
!           CALL SETVEC(C,ZERO,NCROW*NCCOL)
!         END IF
!        END IF
!      END IF
!*
!      IF ((ICONVEX.EQ.1 .OR.IESSL.EQ.1).AND. IZERO.EQ.0 ) THEN
*. DGEMM from CONVEX/ESSL  lib
!        IF (ga_nodeid().eq.0)
!     &   write(*,'(a,4i4)') 'DGEMM,ICONVEX,IESSL,IZERO,ITRNSP',
!     &                             ICONVEX,IESSL,IZERO,ITRNSP
        LDA = MAX(1,NAROW)
        LDB = MAX(1,NBROW)
* 
        LDC = MAX(1,NCROW)
        IF(ITRNSP.EQ.0) THEN
        CALL DGEMM('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        END IF
!      ELSE
!        IF (ga_nodeid().eq.0)
!     &   write(*,'(a,4i4)') 'Jeppe,ICONVEX,IESSL,IZERO,ITRNSP',
!     &                             ICONVEX,IESSL,IZERO,ITRNSP
!* Use Jeppes version ( it should be working )
!        IF( ITRNSP .EQ. 0 ) THEN     
!* ======
!* C=A*B
!* ======
!*
!          IROLL = 6
!          NBREST = MOD(NBROW,IROLL)
!C?        WRITE(6,*) ' IROLL, NBROW, NBREST ', IROLL, NBROW,NBREST
!*. The 'remainder'
!          IF(NBREST.EQ.0) THEN
!            DO J = 1, NCCOL
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) 
!              END DO
!            END DO
!*. Just scale 
!          ELSE IF(NBREST.EQ.1) THEN
!* C(I,J) = FACTOR*C(I,J) + FACTORAB*A(I,1)*B(1,J)
!            DO J = 1, NCCOL
!              B1J = FACTORAB*B(1,J)
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1)
!              END DO
!            END DO
!          ELSE IF (NBREST.EQ.2) THEN
!            DO J = 1, NCCOL
!             B1J =  FACTORAB*B(1,J)
!             B2J =  FACTORAB*B(2,J)
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1) + B2J*A(I,2)
!              END DO
!            END DO
!          ELSE IF (NBREST .EQ.3) THEN
!            DO J = 1, NCCOL
!             B1J =  FACTORAB*B(1,J)
!             B2J =  FACTORAB*B(2,J)
!             B3J =  FACTORAB*B(3,J)
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1) + B2J*A(I,2)
!     &                 + B3J*A(I,3)
!              END DO
!            END DO
!          ELSE IF (NBREST.EQ.4) THEN
!            DO J = 1, NCCOL
!             B1J =  FACTORAB*B(1,J)
!             B2J =  FACTORAB*B(2,J)
!             B3J =  FACTORAB*B(3,J)
!             B4J =  FACTORAB*B(4,J)
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1) + B2J*A(I,2)
!     &                 + B3J*A(I,3) + B4J*A(I,4)
!              END DO
!            END DO
!          ELSE IF (NBREST.EQ.5) THEN
!            DO J = 1, NCCOL
!             B1J =  FACTORAB*B(1,J)
!             B2J =  FACTORAB*B(2,J)
!             B3J =  FACTORAB*B(3,J)
!             B4J =  FACTORAB*B(4,J)
!             B5J =  FACTORAB*B(5,J)
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1) + B2J*A(I,2)
!     &                 + B3J*A(I,3) + B4J*A(I,4) + B5J*A(I,5)
!              END DO
!            END DO
!          ELSE IF (NBREST.EQ.6) THEN
!            DO J = 1, NCCOL
!             B1J =  FACTORAB*B(1,J)
!             B2J =  FACTORAB*B(2,J)
!             B3J =  FACTORAB*B(3,J)
!             B4J =  FACTORAB*B(4,J)
!             B5J =  FACTORAB*B(5,J)
!             B6J =  FACTORAB*B(6,J)
!              DO I = 1, NCROW
!                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1) + B2J*A(I,2)
!     &                 + B3J*A(I,3) + B4J*A(I,4) + B5J*A(I,5)
!     &                 + B6J*A(I,6) 
!              END DO
!            END DO
!          END IF
!*. And then the remaining part of the loop
!          IF(NBREST.NE.NBROW) THEN
!            DO J = 1, NCCOL
!C?           write(6,*) ' nbrest, nbrow, iroll ',
!C?   &                    nbrest, nbrow, iroll
!             DO KOFF = NBREST+1,NBROW-IROLL+1,IROLL
!C?            WRITE(6,*) ' KOFF = ', KOFF
!              B1J =  FACTORAB*B(KOFF,J)
!              B2J =  FACTORAB*B(KOFF+1,J)
!              B3J =  FACTORAB*B(KOFF+2,J)
!              B4J =  FACTORAB*B(KOFF+3,J)
!              B5J =  FACTORAB*B(KOFF+4,J)
!              B6J =  FACTORAB*B(KOFF+5,J)
!              DO I = 1, NCROW
!                C(I,J) = C(I,J) 
!     &                 + B1J*A(I,KOFF) 
!     &                 + B2J*A(I,KOFF+1)
!     &                 + B3J*A(I,KOFF+2) 
!     &                 + B4J*A(I,KOFF+3) 
!     &                 + B5J*A(I,KOFF+4)
!     &                 + B6J*A(I,KOFF+5) 
!              END DO
!            END DO
!           END DO
!          END IF
!*         ^ End of NBREST.NE.NBROW
!        END IF
!*
!* =========
!* C=A(T)*B
!* =========
!*
!        IF ( ITRNSP .EQ. 1 ) THEN      
!          IF(FACTORC.NE.0) THEN
!            DO J = 1, NCCOL
!              DO I = 1, NCROW
!                T = 0.0D0         
!                DO K = 1, NBROW
!                  T = T  + A(K,I)*B(K,J)
!                END DO   
!                C(I,J) = FACTORC*C(I,J) + FACTORAB*T
!              END DO   
!            END DO   
!          ELSE
!            DO J = 1, NCCOL
!              DO I = 1, NCROW
!                T = 0.0D0         
!                DO K = 1, NBROW
!                  T = T  + A(K,I)*B(K,J)
!                END DO   
!                C(I,J) = FACTORAB*T
!              END DO   
!            END DO   
!          END IF
!*         ^ End of FACTORC = 0
!        END IF   
!C
!        IF ( ITRNSP .EQ. 2 ) THEN     
!C ===========
!C. C = A*B(T)
!C ===========
!* Playing around with batching over row index I to 
!* reduce cache overflow
!         LENIBATCH = 200
!         NIBATCH = NCROW/LENIBATCH
!         IF(NIBATCH*LENIBATCH.LT.NCROW) NIBATCH = NIBATCH + 1
!         I_START = 0
!         DO IBATCH = 1, NIBATCH
!*
!          IF(IBATCH.EQ.1) THEN
!            I_START = 1
!            I_STOP = MIN(LENIBATCH,NCROW)
!          ELSE 
!            I_START = I_START + LENIBATCH
!            I_STOP  = MIN(I_START-1+LENIBATCH,NCROW)
!          END IF
!*
!          DO J = 1,NCCOL
!*. Initialization
!            IF(FACTORC.NE.0.0D0) THEN
!              IF(NBCOL.GE.1) THEN
!                BJ1 = FACTORAB*B(J,1)
!                DO I = I_START, I_STOP
!                  C(I,J) = FACTORC*C(I,J) + BJ1*A(I,1)
!                END DO   
!              END IF
!            ELSE IF (FACTORC.EQ.0.0D0) THEN
!              IF(NBCOL.GE.1) THEN
!                BJ1 = FACTORAB*B(J,1)
!                DO I = I_START, I_STOP
!                  C(I,J) =  BJ1*A(I,1)
!                END DO   
!              END IF
!            END IF
!*. And the rest
!            DO K = 2,NBCOL
!              BJK = FACTORAB*B(J,K)
!              DO I = I_START, I_STOP
!                C(I,J) = C(I,J) + BJK*A(I,K)
!              END DO   
!            END DO   
!          END DO   
!          END DO
!*         ^ End of loop over I-batches
!        END IF   
!      END IF
!*     ^ End of switch to home-made
!*
!      IF(ITRNSP.EQ.3) THEN
!C ================
!C. C = A(T)*B(T)
!C ================
!* C(I,J) = FACTORC*C(I,J) + FACTORAB*sum(K) A(K,I)*B(J,K)
!        CALL SCALVE(C,FACTORC,NCROW*NCCOL)
!        DO I = 1, NCROW
!          DO K = 1, NAROW
!            AKI = FACTORAB*A(K,I)
!            DO J = 1,NBROW
!              C(I,J) = C(I,J) + AKI*B(J,K)
!            END DO
!          END DO
!        END DO
!      END IF
C
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' C MATRIX FROM MATML7 ' 
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
!      T_END = 0.01*LTIMEX()
!      TMULT = TMULT + T_END - T_INI 
C
C     CALL LUCIAQEXIT('MATML')
!      cpu_MATML7=cpu_MATML7+util_cpusec()-cpu_MATML7
!      wall_MATML7=wall_MATML7+util_wallsec()-wall_MATML7
      RETURN
      END
c $Id$
