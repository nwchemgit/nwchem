      SUBROUTINE DUCC(RTDB,D_T1,K_T1_OFFSET,D_T2,K_T2_OFFSET,
     &        D_V2,K_V2_OFFSET,ENE_ORB,
     &        NOS,NVS,NOAS,NOBS,NVAS,NVBS,NACTV,
     &        TRANSFORM,
     &        ACTINDEX1,PRINTINDEX1,ACTINDEX2,PRINTINDEX2,
     &        H,HORB,V,
     &        T1,T2,
     &        HT,VT,DUCC_MODEL,NONHF)
#ifndef OPTIMIZE_DUCC
CDIR$ OPTIMIZE: 2
#endif
C
C  ONLY FOR RHF
C
C  ALL VARIABLES BELOW HAVE TO BE DEFINED IN THE TCE_ENERGY.F
C  NACT SHOULD BE DEFINED IN THE INPUT SET TCE:NACTV
C
C  PLEASE DO NOT FREEZE THE OCCUPIED ORBITALS FOR NOW
C  ____________________________________________________________________
C
C  NACTV - NUMBER OF ACTIVE VIRTUAL ORBITALS
C  NOS   - NUMBER OF OCCUPIED SPIN ORBITALS
C  NVS   - NUMBER OF VIRTUAL SPIN ORBITALS
C  NOAS  - NUMBER OF OCCUPIED ALPHA SPIN ORBITALS
C  NOBS  - NUMBER OF OCCUPIED BETA SPIN ORBITALS
C  NVAS  - NUMBER OF VIRTUAL ALPHA SPIN ORBITALS
C  NVBS  - NUMBER OF VIRTUAL BETA SPIN ORBITALS
C
C  SPIN ORBITAL CONVENTION
C
C  | NOAS | NOBS | NVAS | NVBS |
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C  CHANGE ACTINDEX and PRINT_INDEX to logic in tce_energy.f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
C
      INTEGER          :: D_T1, K_T1_OFFSET
      INTEGER          :: D_T2, K_T2_OFFSET
      INTEGER          :: D_V2, K_V2_OFFSET
C
      INTEGER          :: NOS, NVS     ! # OF OCCUPIED/VIRTUAL SPINORBITALS
      INTEGER          :: NOAS, NVAS   ! # OF OCCUPIED/VIRTUAL ALPHA SPINORBITALS
      INTEGER          :: NOBS, NVBS   ! # OF OCCUPIED/VIRTUAL BETA SPINORBITALS
      INTEGER          :: NACTV        ! # OF ACTIVE OCCUPIED SPINORBITALS
C
      INTEGER          :: I, J, K, L   ! FIXED OCCUPIED SPIN ORBITAL INDICIES
      INTEGER          :: A, B, C, D   ! FIXED VIRTUAL SPIN ORBITAL INDICIES
      INTEGER          :: E, F, G, H   ! FREE SPIN ORBITAL INDICIES
      INTEGER          :: M, N, O, P   ! FREE SPIN ORBITAL INDICIES
      INTEGER          :: Q, R, S, T
C
      DOUBLE PRECISION :: PTHRESH  !PRINTING THRESHOLD
      INTEGER          :: TRANSFORM(NOS+NVS)  !RETURNS THE CORRESPONDING ORBITAL LABEL FOR A GIVEN SPIN ORBITAL
C
C     MATRICES WITH ORIGINAL INTEGRALS
      DOUBLE PRECISION :: Fock(NOS+NVS,NOS+NVS)
      DOUBLE PRECISION :: HORB((NOS+NVS)/2,(NOS+NVS)/2) ! VALID ONLY FOR RHF CASE
      DOUBLE PRECISION :: V(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)
C
C     CLUSTER OPERATORS: OCCUPIED IS NUMBERED 1--NOS AND VIRTUAL IS NUMBERED (NOS+1)--(NOS+NVS)
      DOUBLE PRECISION :: T1(NOS,NOS+1:NOS+NVS)                   ! (OCCUPIED,VIRTUAL)
      DOUBLE PRECISION :: T2(NOS,NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS) ! (OCCUPIED,OCCUPIED,VIRTUAL,VIRTUAL)
C
C     TRANSFORMED MATRICES
      DOUBLE PRECISION :: HT(NOS+NVS,NOS+NVS)                 ! TRANSFORMED 1E- INTEGRALS
      DOUBLE PRECISION :: VT(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS) ! TRANSFORMED 2E- INTEGRALS
      DOUBLE PRECISION :: FOCKT(NOS+NVS,NOS+NVS)              ! TRANSFORMED FOCK INTEGRALS
C
      INTEGER          :: ACTINDEX1(NOS+NVS,NOS+NVS)
      INTEGER          :: PRINTINDEX1(NOS+NVS,NOS+NVS)
      LOGICAL*1        :: ACTINDEX2(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)
      LOGICAL*1        :: PRINTINDEX2(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)
      INTEGER          :: DUCC_MODEL
      LOGICAL          :: NONHF
C
      DOUBLE PRECISION :: ENE_ORB(NOS+NVS)  ! A COPY OF DBL_MB(K_SORTED)
C
      DOUBLE PRECISION :: ZZZ
      INTEGER          :: VALUE1
      INTEGER          :: VALUE2
C
      REAL :: TIME_START, TIME_FINISH
C
      INTEGER RTDB
      LOGICAL NODEZERO
C
      DOUBLE PRECISION :: sum, sum_total
C ___________________________________________________________________________________________________________
C     INITIALIZATION
C
      NODEZERO = (ga_nodeid().EQ.0)
C
      IF(NODEZERO) THEN
         WRITE(6,*)'DUCC MODEL = ',DUCC_MODEL
         WRITE(6,*)'NONHF = ',NONHF
         CALL UTIL_FLUSH(6)
      ENDIF
C
C     INITIALIZE/ZERO OUT HORB, H, AND HT
      DO Q = 1, (NOS+NVS)/2
         DO R = 1, (NOS+NVS)/2
            HORB(Q,R) = 0.0D0
      ENDDO; ENDDO
C
      DO Q = 1, NOS+NVS
         DO R = 1, NOS+NVS
            Fock(Q,R)  = 0.0D0
            HT(Q,R) = 0.0D0
      ENDDO; ENDDO
C
C     CONSTRUCT H FROM HORB
      CALL KINETIC_HCORE_1(RTDB,HORB,FOCK,NOS,NVS,NOAS,NOBS,NVAS,NVBS)
C
C     INITIALIZE/ZERO OUT V AND VT
      DO Q = 1, NOS+NVS
         DO R = 1, NOS+NVS
            DO S = 1, NOS+NVS
               DO T = 1, NOS+NVS
                  V(Q,R,S,T)  = 0.0D0
                  VT(Q,R,S,T) = 0.0D0
      ENDDO; ENDDO; ENDDO; ENDDO
C
C     CONSTRUCT V
      CALL MAPPING_V2_M(RTDB,D_V2,K_V2_OFFSET,V,NOS,NVS)
C
C     INITIALIZE/ZERO T1 AND T2
      DO I = 1, NOS
         DO A = NOS+1, NOS+NVS
            T1(I,A) = 0.0D0
      ENDDO; ENDDO
C
      DO I = 1, NOS
         DO J = 1, NOS
            DO A = NOS+1, NOS+NVS
               DO B = NOS+1, NOS+NVS
                  T2(I,J,A,B) = 0.0D0
      ENDDO; ENDDO; ENDDO; ENDDO
C
C     CONSTRUCT/FILL T1 AND T2
      CALL MAPPING_T1(RTDB,D_T1,K_T1_OFFSET,T1,NOS,NVS)
      CALL MAPPING_T2(RTDB,D_T2,K_T2_OFFSET,T2,NOS,NVS)
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHANGE 1e INTEGRALS TO FOCK MATRIX

      DO A = NOS+1, NOS+NVS
        DO B = NOS+1, NOS+NVS
          DO M = 1, NOS
            Fock(A,B) = Fock(A,B) + V(M,A,M,B)
      ENDDO; ENDDO; ENDDO

      DO I = 1, NOS
        DO A = NOS+1, NOS+NVS
          DO M = 1, NOS
            Fock(I,A) = Fock(I,A) + V(M,I,M,A)
      ENDDO; ENDDO; ENDDO

      DO A = NOS+1, NOS+NVS
        DO I = 1, NOS
          DO M = 1, NOS
            Fock(A,I) = Fock(A,I) + V(M,A,M,I)
      ENDDO; ENDDO; ENDDO

      DO I = 1, NOS
        DO J = 1, NOS
          DO M = 1, NOS
            Fock(I,J) = Fock(I,J) + V(M,I,M,J)
      ENDDO; ENDDO; ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sum_total = 0.0
      sum = 0.0
      DO I = 1, NOS
         DO A = NOS+1, NOS+NVS
            sum = sum + FOCK(A,I)*T1(I,A)
      ENDDO; ENDDO
      write(*,*)"F-T1 correlation    = ",sum
      sum_total = sum_total + sum
      sum = 0.0
      DO I = 1, NOS
         DO J = 1, NOS
            DO A = NOS+1, NOS+NVS
               DO B = NOS+1, NOS+NVS
                  sum = sum + 0.5*V(I,J,A,B)*T1(I,A)*T1(J,B)
      ENDDO; ENDDO; ENDDO; ENDDO
      write(*,*)"V-T1-T1 correlation = ",sum
      sum_total = sum_total + sum
      sum = 0.0
      DO I = 1, NOS
         DO J = 1, NOS
            DO A = NOS+1, NOS+NVS
               DO B = NOS+1, NOS+NVS
                  sum = sum + 0.25*V(I,J,A,B)*T2(I,J,A,B)
      ENDDO; ENDDO; ENDDO; ENDDO
      write(*,*)"V-T2 correlation    = ",sum
      sum_total = sum_total + sum
      write(*,*)"Total correlation   = ", sum_total
C ___________________________________________________________________________________________________________
C     FORM EXTERNAL CLUSTER AMPLITUDES (T1E AND T2E)
C     WE ARE USING THE SAME T1 AND T2 TENSORS TO STORE T1E AND T2E
C
C     ATTENTION: T1E AND T2E CALCULATED FOR SPECIFIC SES G(NO,NACTV) !!!!!!!!

C     TRANSFORMATION T1-->T1E
      DO I = 1, NOS
C      ALPHA
         DO A=NOS+1,NOS+NACTV
            T1(I,A)=0.0D0
         ENDDO
C      BETA
         DO A=NOS+NVAS+1,NOS+NVAS+NACTV
            T1(I,A)=0.0D0
         ENDDO
      ENDDO
C
C TRANSFORMATION T2-->T2E  !
      DO I=1,NOS
         DO J=1,NOS
C           VIRT ALPHA ALPHA
            DO A=NOS+1,NOS+NACTV
               DO B=NOS+1,NOS+NACTV
                  T2(I,J,A,B)=0.0D0
            ENDDO; ENDDO
C
C           VIRT ALPHA BETA  & BETA ALPHA
            DO A=NOS+1,NOS+NACTV
               DO B=NOS+NVAS+1,NOS+NVAS+NACTV
                  T2(I,J,A,B)=0.0D0
                  T2(I,J,B,A)=0.0D0
            ENDDO; ENDDO
C
C           VIRT BETA BETA
            DO A=NOS+NVAS+1,NOS+NVAS+NACTV
               DO B=NOS+NVAS+1,NOS+NVAS+NACTV
                  T2(I,J,A,B)=0.0D0
            ENDDO; ENDDO
C
          ENDDO
       ENDDO
C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sum_total = 0.0
      sum = 0.0
      DO I = 1, NOS
         DO A = NOS+1, NOS+NVS
            sum = sum + FOCK(A,I)*T1(I,A)
      ENDDO; ENDDO
      write(*,*)"F-T1(Ext) correlation         = ",sum
      sum_total = sum_total + sum
      sum = 0.0
      DO I = 1, NOS
         DO J = 1, NOS
            DO A = NOS+1, NOS+NVS
               DO B = NOS+1, NOS+NVS
                  sum = sum + 0.5*V(I,J,A,B)*T1(I,A)*T1(J,B)
      ENDDO; ENDDO; ENDDO; ENDDO
      write(*,*)"V-T1(Ext)-T1(Ext) correlation = ",sum
      sum_total = sum_total + sum
      sum = 0.0
      DO I = 1, NOS
         DO J = 1, NOS
            DO A = NOS+1, NOS+NVS
               DO B = NOS+1, NOS+NVS
                  sum = sum + 0.25*V(I,J,A,B)*T2(I,J,A,B)
      ENDDO; ENDDO; ENDDO; ENDDO
      write(*,*)"V-T2(Ext) correlation         = ",sum
      sum_total = sum_total + sum
      write(*,*)"Total Ext correlation         = ", sum_total


      CALL CPU_TIME(TIME_START)
      CALL HH(NOS,NVS,T1,T2,FOCK,V,DUCC_MODEL,NONHF,HT)
      CALL HP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,HT)
      CALL PP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,HT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("1-BODY TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START
      CALL CPU_TIME(TIME_START)
      CALL HHHH(NOS,NVS,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("  HHHH TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START
      CALL CPU_TIME(TIME_START)
      CALL PPPP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("  PPPP TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START
      CALL CPU_TIME(TIME_START)
      CALL HHHP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("  HHHP TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START
      CALL CPU_TIME(TIME_START)
      CALL PPPH(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("  PPPH TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START
      CALL CPU_TIME(TIME_START)
      CALL HPHP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("  HPHP TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START
      CALL CPU_TIME(TIME_START)
      CALL HHPP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
      CALL CPU_TIME(TIME_FINISH)
      PRINT '("  HHPP TIME = ",f8.3," SECONDS.")',TIME_FINISH-TIME_START

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     FINAL CHECK FOR UNITARY and ANTISYMMETRY ****
C     UNCOMMENTING THE FOLLOWING LOOPS WILL TEST THE UNITARY AND ANTISYMETRY
C     CHARACTER OF THE TRANSFORMED MATRIX ELEMENTS. IDEALLY, THE FOLLOWING OUGHT
C     TO BE ZERO OR NEAR ZERO. IT WILL PRINT OUT OTHERWISE. THE THRESHOLD IS SO
C     THAT IT PRINT OUT 'SOMETHING', BUT THIS WILL BE SMALL IF THERE IS NO PROBLEM.
C     DEBUG
C
      IF(NODEZERO) THEN
C
         WRITE(6,*) 'ONE-BODY UNITARY CHECK'
         DO P = 1, NOS+NVS
         IF(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS))CYCLE
         IF(P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS))CYCLE
            DO Q = 1, NOS+NVS
         IF(Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS))CYCLE
         IF(Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS))CYCLE
               ZZZ=HT(P,Q)-HT(Q,P)
               IF(abs(ZZZ).GT.(1.0d-12)) THEN
                  WRITE(6,*)P,Q,ZZZ
               ENDIF
         ENDDO; ENDDO
C
         WRITE(6,*) 'TWO-BODY UNITARY CHECK'
C
         DO P = 1, NOS+NVS
         IF(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS))CYCLE
         IF(P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS))CYCLE
            DO Q = 1, NOS+NVS
         IF(Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS))CYCLE
         IF(Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS))CYCLE
               DO R = 1, NOS+NVS
         IF(R.GT.(NOS+NACTV).AND.R.LE.(NOS+NVAS))CYCLE
         IF(R.GT.(NOS+NVAS+NACTV).AND.R.LE.(NOS+NVS))CYCLE
                  DO S = 1, NOS+NVS
         IF(S.GT.(NOS+NACTV).AND.S.LE.(NOS+NVAS))CYCLE
         IF(S.GT.(NOS+NVAS+NACTV).AND.S.LE.(NOS+NVS))CYCLE
                     ZZZ=VT(P,Q,R,S)-VT(R,S,P,Q)
                     IF(abs(ZZZ).GT.(1.0d-12)) THEN
                        WRITE(6,*)P,Q,R,S,ZZZ
                     ENDIF
         ENDDO; ENDDO; ENDDO; ENDDO
C
         WRITE(6,*) 'TWO-BODY ANTISYMMETRY CHECK'
         DO P = 1, NOS+NVS
         IF(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS))CYCLE
         IF(P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS))CYCLE
            DO Q = 1, NOS+NVS
         IF(Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS))CYCLE
         IF(Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS))CYCLE
               DO R = 1, NOS+NVS
         IF(R.GT.(NOS+NACTV).AND.R.LE.(NOS+NVAS))CYCLE
         IF(R.GT.(NOS+NVAS+NACTV).AND.R.LE.(NOS+NVS))CYCLE
                  DO S = 1, NOS+NVS
         IF(S.GT.(NOS+NACTV).AND.S.LE.(NOS+NVAS))CYCLE
         IF(S.GT.(NOS+NVAS+NACTV).AND.S.LE.(NOS+NVS))CYCLE
                     ZZZ=VT(P,Q,R,S)+VT(Q,P,R,S)
                     IF(abs(ZZZ).GT.(1.0d-12)) THEN
                        WRITE(6,*)P,Q,R,S,ZZZ
                     ENDIF
         ENDDO; ENDDO; ENDDO; ENDDO
C
      ENDIF
C     END DEBUG
C
C ___________________________________________________________________________________________________________
C
      FOCKT=HT

      DO A = NOS+1, NOS+NVS
        DO B = NOS+1, NOS+NVS
          DO M = 1, NOS
            HT(A,B) = HT(A,B) - VT(M,A,M,B)
      ENDDO; ENDDO; ENDDO

      DO I = 1, NOS
        DO A = NOS+1, NOS+NVS
          DO M = 1, NOS
            HT(I,A) = HT(I,A) - VT(M,I,M,A)
      ENDDO; ENDDO; ENDDO

      DO A = NOS+1, NOS+NVS
        DO I = 1, NOS
          DO M = 1, NOS
            HT(A,I) = HT(A,I) - VT(M,A,M,I)
      ENDDO; ENDDO; ENDDO

      DO I = 1, NOS
        DO J = 1, NOS
          DO M = 1, NOS
            HT(I,J) = HT(I,J) - VT(M,I,M,J)
      ENDDO; ENDDO; ENDDO
C
C ___________________________________________________________________________________________________________
C  _______   _______   ______  __    __  ________
C |       \ |       \ |      \|  \  |  \|        \
C | $$$$$$$\| $$$$$$$\ \$$$$$$| $$\ | $$ \$$$$$$$$
C | $$__/ $$| $$__| $$  | $$  | $$$\| $$   | $$
C | $$    $$| $$    $$  | $$  | $$$$\ $$   | $$
C | $$$$$$$ | $$$$$$$\  | $$  | $$\$$ $$   | $$
C | $$      | $$  | $$ _| $$_ | $$ \$$$$   | $$
C | $$      | $$  | $$|   $$ \| $$  \$$$   | $$
C  \$$       \$$   \$$ \$$$$$$ \$$   \$$    \$$

C  THIS IS SET UP WITH ALL OCCUPIED ORBITALS BEING ACTIVE ALONG WITH
C  A SUBSET OF VIRTUALS, WHICH IS INDICATED BY THE INPUT VARIABLE NACTV
C  AND COUNTED STARTING AT THE LOWEST ENERGY UNOCCUPIED ORBITAL.

! 900   CONTINUE
C PRINTING THRESHOLD
      PTHRESH = 1.0d-10

C SPIN ORBITAL TO ORBITAL TRANSFORM. USED WHEN PRINTING.

      Q = 0

      DO P = 1, NOAS
        Q=Q+1
        TRANSFORM(Q) = P
      ENDDO

      DO P = NOAS+1, NOS
        Q=Q+1
        TRANSFORM(Q) = P-NOAS
      ENDDO

      DO P = NOS+1, NOS+NVAS
        Q=Q+1
        TRANSFORM(Q) = P-NOS+NOAS
      ENDDO

      DO P = NOS+NVAS+1, NOS+NVS
        Q=Q+1
        TRANSFORM(Q) = P-NOS-NVAS+NOAS
      ENDDO

C ACTINDEX1 IS AN ARRAY TO TEST IF A SET OF INDICES {P,Q} CONTAIN
C AT LEAST ONE EXTERNAL(INACTIVE) INDEX. IF SO, THEN THAT SET OF
C INDICES IS 'EXTERNAL' AND ACTINDEX1=1, ELSE ACTINDEX1=0.

      ACTINDEX1 = 0

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          IF(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS))ACTINDEX1(P,Q)=1
          IF(Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS))ACTINDEX1(P,Q)=1
          IF(P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS))ACTINDEX1(P,Q)=1
          IF(Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS))ACTINDEX1(P,Q)=1
      ENDDO; ENDDO

      PRINTINDEX1 = 0

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
C SKIP PRINTING IF P OR Q IS EXTERNAL
          IF(ACTINDEX1(P,Q).EQ.1)PRINTINDEX1(P,Q)=1
C SKIP PRINTING IF P OR Q IS BETA
          IF(P.GT.NOAS.AND.P.LE.NOS)PRINTINDEX1(P,Q)=1
          IF(Q.GT.NOAS.AND.Q.LE.NOS)PRINTINDEX1(P,Q)=1
          IF(P.GT.NOS+NVAS)PRINTINDEX1(P,Q)=1
          IF(Q.GT.NOS+NVAS)PRINTINDEX1(P,Q)=1
      ENDDO; ENDDO

      IF(NODEZERO) THEN
        WRITE(6,*)'TRANSFORMED ~FOCK~ MATRIX'
        CALL UTIL_FLUSH(6)
      ENDIF

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
          IF(abs(FOCKT(P,Q)).GT.PTHRESH)THEN
             IF(NODEZERO) THEN
               WRITE(6,'(2I5,3X,F18.10)')
     &         TRANSFORM(P),TRANSFORM(Q),FOCKT(P,Q)
               CALL UTIL_FLUSH(6)
             ENDIF
          ENDIF
      ENDDO; ENDDO

      IF(NODEZERO) THEN
        WRITE(6,*)'begin_one_electron_integrals'
        CALL UTIL_FLUSH(6)
      ENDIF

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
          IF(abs(HT(P,Q)).GT.PTHRESH)THEN
            IF(NODEZERO) THEN
              WRITE(6,'(2I5,3X,F18.10)')
     &        TRANSFORM(P),TRANSFORM(Q),HT(P,Q)
              CALL UTIL_FLUSH(6)
            ENDIF
          ENDIF
      ENDDO; ENDDO

      IF(NODEZERO) THEN
        WRITE(6,*)'end_one_electron_integrals'
        CALL UTIL_FLUSH(6)
      ENDIF

C ACTINDEX2 IS AN ARRAY TO TEST IF A SET OF INDICES {P,Q,R,S} CONTAIN
C AT LEAST ONE EXTERNAL(INACTIVE) INDEX. IF SO, THEN THAT SET OF
C INDICES IS 'EXTERNAL' AND ACTINDEX2=1, ELSE ACTINDEX2=0.

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          DO R = 1, NOS+NVS
            DO S = 1, NOS+NVS
      ACTINDEX2(P,Q,R,S)=(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS)).OR.
     O              (Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS)).OR.
     O              (R.GT.(NOS+NACTV).AND.R.LE.(NOS+NVAS)).OR.
     O              (S.GT.(NOS+NACTV).AND.S.LE.(NOS+NVAS)).OR.
     O              (P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS)).OR.
     O              (Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS)).OR.
     O              (R.GT.(NOS+NVAS+NACTV).AND.R.LE.(NOS+NVS)).OR.
     O              (S.GT.(NOS+NVAS+NACTV).AND.S.LE.(NOS+NVS))
      ENDDO; ENDDO; ENDDO; ENDDO


      DO P = 1, NOS+NVS
        DO R = 1, NOS+NVS
          DO Q = 1, NOS+NVS
            DO S = 1, NOS+NVS
C SKIP PRINTING IF P, Q, R, OR S IS EXTERNAL
               PRINTINDEX2(P,Q,R,S)=ACTINDEX2(P,Q,R,S).OR.
C SKIP PRINTING IF P OR R IS BETA
     O              (P.GT.NOAS.AND.P.LE.NOS).OR.
     O              (R.GT.NOAS.AND.R.LE.NOS).OR.
     O              (P.GT.NOS+NVAS).OR.
     O              (R.GT.NOS+NVAS).OR.
C SKIP PRINTING IF Q OR S IS ALPHA
     O              (Q.LE.NOAS).OR.
     O              (S.LE.NOAS).OR.
     O              (Q.GT.NOS.AND.Q.LE.NOS+NVAS).OR.
     O              (S.GT.NOS.AND.S.LE.NOS+NVAS)
      ENDDO; ENDDO; ENDDO; ENDDO

      IF(NODEZERO) THEN
        WRITE(6,*)'begin_two_electron_integrals'
        CALL UTIL_FLUSH(6)
      ENDIF

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          DO R = 1, NOS+NVS
            DO S = 1, NOS+NVS
              IF(PRINTINDEX2(P,Q,R,S))CYCLE
              IF(abs(VT(P,Q,R,S)).GT.PTHRESH)THEN
                 IF(NODEZERO) THEN
                   WRITE(6,'(4I5,2X,F18.10)')TRANSFORM(P),
     &             TRANSFORM(R),TRANSFORM(Q),TRANSFORM(S),
     &             VT(P,Q,R,S)
                   CALL UTIL_FLUSH(6)
                 ENDIF
              ENDIF
      ENDDO; ENDDO; ENDDO; ENDDO

      IF(NODEZERO) THEN
        WRITE(6,*)'end_two_electron_integrals'
        CALL UTIL_FLUSH(6)
      ENDIF


C ADDITIONAL PRINTING
C ___________________________________________________________________________________________________________
C PRINTING FOR EXTERNAL CCSDTQ
C   !    IF(NODEZERO) THEN
C   !      WRITE(6,*)'onebody.inp'
C   !      CALL UTIL_FLUSH(6)
C   !    ENDIF
C
C      I=0
C      DO P = 1, NOS+NVS
C        DO Q = 1, P
C          IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
C          I=I+1
C          IF(NODEZERO) THEN
C            WRITE(990,'(F18.10,3X,I3)')
C     &      HT(P,Q),I
C            CALL UTIL_FLUSH(6)
C          ENDIF
C      ENDDO; ENDDO
C
C   !    IF(NODEZERO) THEN
C   !      WRITE(6,*)'twobody.inp'
C   !      CALL UTIL_FLUSH(6)
C   !    ENDIF
C      DO P = 1, NOS+NVS
C        DO R = 1, NOS+NVS
C          DO Q = 1, NOS+NVS
C            DO S = 1, NOS+NVS
C              IF(PRINTINDEX2(P,Q,R,S))CYCLE
C              IF(abs(VT(P,Q,R,S)).GT.PTHRESH)THEN
C                IF(NODEZERO) THEN
C                  WRITE(991,'(4I5,2X,F18.10)')TRANSFORM(P),
C     &            TRANSFORM(R),TRANSFORM(Q),TRANSFORM(S),VT(P,Q,R,S)
C                  CALL UTIL_FLUSH(6)
C                ENDIF
C              ENDIF
C      ENDDO; ENDDO; ENDDO; ENDDO
c
      return
      end
c
