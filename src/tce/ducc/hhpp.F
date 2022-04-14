      SUBROUTINE HHPP(NOS,NVS,NVAS,NACTV,T1,T2,
     &                FOCK,V,DUCC_MODEL,NONHF,VT)
#ifndef OPTIMIZE_DUCC
CDIR$ OPTIMIZE: 2
#endif
      IMPLICIT NONE
      INTEGER          :: A, B, C, D   ! FIXED VIRTUAL SPIN ORBITAL INDICIES
      INTEGER          :: E, F, G, H   ! FREE VIRTUAL SPIN ORBITAL INDICIES
      INTEGER          :: I, J, K, L   ! FIXED OCCUPIED SPIN ORBITAL INDICIES
      INTEGER          :: M, N, O, P   ! FREE OCCUPIED SPIN ORBITAL INDICIES
      INTEGER          :: NOS, NVS     ! # OF OCCUPIED/VIRTUAL SPINORBITALS
      INTEGER          :: NVAS         ! # OF ALPHA VIRTUAL SPINORBITALS
      INTEGER          :: NACTV        ! # OF ACTIVE VIRTUAL ORBITALS
      INTEGER          :: DUCC_MODEL
      LOGICAL          :: NONHF
      DOUBLE PRECISION :: FOCK(NOS+NVS,NOS+NVS)
      DOUBLE PRECISION :: V(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)
      DOUBLE PRECISION :: VT(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)     ! TRANSFORMED 2E- INTEGRALS
C
C     CLUSTER OPERATORS: OCCUPIED IS NUMBERED 1--NOS AND VIRTUAL IS NUMBERED (NOS+1)--(NOS+NVS)
      DOUBLE PRECISION :: T1(NOS,NOS+1:NOS+NVS)                   ! (OCCUPIED,VIRTUAL)
      DOUBLE PRECISION :: T2(NOS,NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS) ! (OCCUPIED,OCCUPIED,VIRTUAL,VIRTUAL)
C
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVV(:,:,:,:)

!
c     1*V{aibj}
c      Commutator = 0, PT Order = 1
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) + V(a,b,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*F{mi}*T2{ajbm}
! c      Commutator = 1, PT Order = 1
! !
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      & FOCK(m,i)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     1*F{be}*T2{aiej}
c      Commutator = 1, PT Order = 1
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     & FOCK(b,e)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*F{ae}*T2{biej}
c      Commutator = 1, PT Order = 1
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) -
     &          (1.0)*FOCK(a,e)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*F{mj}*T2{aibm}
! c      Commutator = 1, PT Order = 1
! !
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) -
!     & (1.0)*FOCK(m,j)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     -1*V{ajme}*T2{biem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) -
     &          (1.0)*V(a,m,j,e)*T2(i,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{aime}*T2{bjem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     & V(a,m,i,e)*T2(j,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*V{bime}*T2{ajem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) -
     &          (1.0)*V(b,m,i,e)*T2(j,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*V{aebf}*T2{eifj}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     & (1.0/2.0)*V(a,b,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*V{minj}*T2{ambn}
! c      Commutator = 1, PT Order = 2
! !
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      & (1.0/2.0)*V(m,n,i,j)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     1*V{bjme}*T2{aiem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     & V(b,m,j,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*V{bimj}*T1{am}
! c      Commutator = 1, PT Order = 3
! !
!         IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      & V(b,m,i,j)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     -1*V{ajbe}*T1{ei}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) -
     &          (1.0)*V(a,b,j,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{aibe}*T1{ej}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     & V(a,b,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*V{aimj}*T1{bm}
! c      Commutator = 1, PT Order = 3
! !
!         IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) -
!     & (1.0)*V(a,m,i,j)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     -1/2*( (V{menf}*T2{bmfn}) * T2{aiej})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,b) = INTERMVV(e,b) -
     &          (1.0/2.0)*V(m,n,e,f)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(e,b)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1*( (V{menf}*T2{aiem}) * T2{bjfn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,i) = INTERMOVVO(n,f,a,i) +
     & V(m,n,e,f)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,a,i)*T2(j,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (V{menf}*T2{eifm}) * T2{ajbn})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) -
!     & (1.0/2.0)*V(m,n,e,f)*T2(i,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,i)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( (T2+{menf}*V{eifm}) * T2{ajbn})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) -
!     & (1.0/4.0)*T2(m,n,e,f)*V(e,f,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,i)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -1*( (V{menf}*T2{biem}) * T2{ajfn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,i) = INTERMOVVO(n,f,b,i) -
     &          (1.0)*V(m,n,e,f)*T2(i,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,b,i)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1/2*( (T2+{menf}*V{ajem}) * T2{bifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,j) = INTERMOVVO(n,f,a,j) -
     &          (1.0/2.0)*T2(m,n,e,f)*V(a,e,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,a,j)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{bmfn}) * V{aiej})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,b) = INTERMVV(e,b) -
     &          (1.0/4.0)*T2(m,n,e,f)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(e,b)*V(a,e,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( (T2+{menf}*V{ejfm}) * T2{aibn})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      & (1.0/4.0)*T2(m,n,e,f)*V(e,f,j,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,j)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     1/2*( (T2+{menf}*V{aiem}) * T2{bjfn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,i) = INTERMOVVO(n,f,a,i) +
     & (1.0/2.0)*T2(m,n,e,f)*V(a,e,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,a,i)*T2(j,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1/2*( (T2+{menf}*V{biem}) * T2{ajfn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,i) = INTERMOVVO(n,f,b,i) -
     &          (1.0/2.0)*T2(m,n,e,f)*V(b,e,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,b,i)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1/2*( (T2+{menf}*V{bjem}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,j) = INTERMOVVO(n,f,b,j) +
     & (1.0/2.0)*T2(m,n,e,f)*V(b,e,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,b,j)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (V{menf}*T2{ejfm}) * T2{aibn})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      & (1.0/2.0)*V(m,n,e,f)*T2(j,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,j)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     1/4*( (T2+{menf}*V{bmen}) * T2{aifj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,b) = INTERMVV(f,b) +
     & (1.0/4.0)*T2(m,n,e,f)*V(b,e,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(f,b)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/2*( (V{menf}*T2{amfn}) * T2{biej})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) +
     & (1.0/2.0)*V(m,n,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(e,a)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/4*( (T2+{menf}*T2{eifn}) * V{ajbm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) +
     & (1.0/4.0)*T2(m,n,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOO(m,i)*V(a,b,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{menf}*V{amen}) * T2{bifj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,a) = INTERMVV(f,a) -
     &          (1.0/4.0)*T2(m,n,e,f)*V(a,e,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(f,a)*T2(i,j,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{ejfn}) * V{aibm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,j) = INTERMOO(m,j) -
     &          (1.0/4.0)*T2(m,n,e,f)*T2(j,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOO(m,j)*V(a,b,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/8*( (T2+{menf}*T2{ambn}) * V{eifj})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) +
!      & (1.0/8.0)*T2(m,n,e,f)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(e,f,a,b)*V(e,f,i,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( (V{menf}*T2{ambn}) * T2{eifj})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) +
!      & (1.0/4.0)*V(m,n,e,f)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(e,f,a,b)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!
c     1/4*( (T2+{menf}*T2{amfn}) * V{biej})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) +
     & (1.0/4.0)*T2(m,n,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(e,a)*V(b,e,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/8*( (T2+{menf}*V{ambn}) * T2{eifj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) +
     & (1.0/8.0)*T2(m,n,e,f)*V(a,b,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVVV(e,f,a,b)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/2*( (T1+{me}*F{bm}) * T2{aiej})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(e,b) = INTERMVV(e,b) -
     &          (1.0/2.0)*T1(m,e)*FOCK(b,m)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(e,b)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (F{me}*T1{bm}) * T2{aiej})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) -
!     & (1.0)*FOCK(m,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(e,b)*T2(i,j,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (F{me}*T1{ei}) * T2{ajbm})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) +
!      & FOCK(m,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(m,i)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     1/2*( (T1+{me}*F{am}) * T2{biej})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) +
     & (1.0/2.0)*T1(m,e)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(e,a)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (F{me}*T1{am}) * T2{biej})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      & FOCK(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE

!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(e,a)*T2(i,j,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T1+{me}*F{ei}) * T2{ajbm})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) +
!      & (1.0/2.0)*T1(m,e)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(m,i)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{me}*F{ej}) * T2{aibm})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,j) = INTERMOO(m,j) -
!     & (1.0/2.0)*T1(m,e)*FOCK(e,j)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(m,j)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (F{me}*T1{ej}) * T2{aibm})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,j) = INTERMOO(m,j) -
!     & (1.0)*FOCK(m,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(m,j)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     -1*( (T1+{me}*V{bmef}) * T2{aifj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,f) = INTERMVV(b,f) -
     &          (1.0)*T1(m,e)*V(b,e,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,f)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1*( (T1+{me}*V{bief}) * T2{ajfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,i,f) = INTERMOVOV(m,b,i,f) -
     &          (1.0)*T1(m,e)*V(b,e,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,b,i,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{me}*V{einj}) * T2{ambn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
!         INTERMOOOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOO(m,i,n,j) = INTERMOOOO(m,i,n,j) +
!      & T1(m,e)*V(e,n,i,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOO(m,i,n,j)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOO)
!         ENDIF
!
c     1/2*( (T1+{me}*T2{eifm}) * V{ajbf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) +
     & (1.0/2.0)*T1(m,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOV(i,f)*V(a,b,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1*( (T1+{me}*T2{eifj}) * V{ambf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,i,f,j) = INTERMOOVO(m,i,f,j) -
     &          (1.0)*T1(m,e)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOOVO(m,i,f,j)*V(a,b,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1*( (V{bemf}*T1{ei}) * T2{ajfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(b,m,f,i) = INTERMVOVO(b,m,f,i) -
     &          (1.0)*V(b,m,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(b,m,f,i)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1*( (T1+{me}*V{aief}) * T2{bjfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,i,f) = INTERMOVOV(m,a,i,f) +
     & T1(m,e)*V(a,e,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,a,i,f)*T2(j,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/2*( (T1+{me}*T2{amen}) * V{binj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(a,n) = INTERMVO(a,n) -
     &          (1.0/2.0)*T1(m,e)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVO(a,n)*V(b,n,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1*( (T1+{me}*V{amef}) * T2{bifj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(a,f) = INTERMVV(a,f) +
     & T1(m,e)*V(a,e,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,f)*T2(i,j,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1*( (T1+{me}*T2{bjen}) * V{ainm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,b,j,n) = INTERMOVOO(m,b,j,n) -
     &          (1.0)*T1(m,e)*T2(j,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOO(m,b,j,n)*V(a,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{mjne}*T1{am}) * T2{bien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(j,n,e,a) = INTERMOOVV(j,n,e,a) +
!      & V(m,n,j,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOVV(j,n,e,a)*T2(i,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
!
c     1*( (V{aemf}*T1{ei}) * T2{bjfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,m,f,i) = INTERMVOVO(a,m,f,i) +
     & V(a,m,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(a,m,f,i)*T2(j,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (T1+{me}*T2{bmen}) * V{ainj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(b,n) = INTERMVO(b,n) +
     & (1.0/2.0)*T1(m,e)*T2(m,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVO(b,n)*V(a,n,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{mine}*T1{em}) * T2{ajbn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(i,n) = INTERMOO(i,n) -
!     & (1.0)*V(m,n,i,e)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(i,n)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -1/2*( (T1+{me}*T2{ejfm}) * V{aibf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) -
     &          (1.0/2.0)*T1(m,e)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOV(j,f)*V(a,b,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (V{mine}*T1{ej}) * T2{ambn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
!         INTERMOOOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOO(m,i,n,j) = INTERMOOOO(m,i,n,j) +
!      & (1.0/2.0)*V(m,n,i,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOO(m,i,n,j)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOO)
!         ENDIF
!
c     -1*( (T1+{me}*T2{aien}) * V{bjnm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,a,i,n) = INTERMOVOO(m,a,i,n) -
     &          (1.0)*T1(m,e)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOO(m,a,i,n)*V(b,n,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1*( (V{bemf}*T1{ej}) * T2{aifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(b,m,f,j) = INTERMVOVO(b,m,f,j) +
     & V(b,m,e,f)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(b,m,f,j)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1*( (V{aemf}*T1{ej}) * T2{bifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,m,f,j) = INTERMVOVO(a,m,f,j) -
     &          (1.0)*V(a,m,e,f)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(a,m,f,j)*T2(i,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1*( (T1+{me}*T2{bien}) * V{ajnm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,b,i,n) = INTERMOVOO(m,b,i,n) +
     & T1(m,e)*T2(i,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOO(m,b,i,n)*V(a,n,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{mine}*T1{am}) * T2{bjen})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(i,n,e,a) = INTERMOOVV(i,n,e,a) -
!     & (1.0)*V(m,n,i,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOVV(i,n,e,a)*T2(j,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
!
c     1*( (T1+{me}*T2{ajen}) * V{binm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,a,j,n) = INTERMOVOO(m,a,j,n) +
     & T1(m,e)*T2(j,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOO(m,a,j,n)*V(b,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (T1+{me}*V{einm}) * T2{ajbn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(i,n) = INTERMOO(i,n) -
!     & (1.0)*T1(m,e)*V(e,n,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(i,n)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{mjne}*T1{em}) * T2{aibn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(j,n) = INTERMOO(j,n) +
!      & V(m,n,j,e)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(j,n)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -1*( (V{bemf}*T1{em}) * T2{aifj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,f) = INTERMVV(b,f) -
     &          (1.0)*V(b,m,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,f)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (V{mjne}*T1{ei}) * T2{ambn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
!         INTERMOOOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOO(m,j,n,i) = INTERMOOOO(m,j,n,i) -
!     & (1.0/2.0)*V(m,n,j,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOO(m,j,n,i)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{me}*V{ejnm}) * T2{aibn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(j,n) = INTERMOO(j,n) +
!      & T1(m,e)*V(e,n,j,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(j,n)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     1*( (V{aemf}*T1{em}) * T2{bifj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(a,f) = INTERMVV(a,f) +
     & V(a,m,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,f)*T2(i,j,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1*( (T1+{me}*V{ajef}) * T2{bifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,j,f) = INTERMOVOV(m,a,j,f) -
     &          (1.0)*T1(m,e)*V(a,e,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,a,j,f)*T2(i,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{mine}*T1{bm}) * T2{ajen})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(i,n,e,b) = INTERMOOVV(i,n,e,b) +
!      & V(m,n,i,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOVV(i,n,e,b)*T2(j,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{mjne}*T1{bm}) * T2{aien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(j,n,e,b) = INTERMOOVV(j,n,e,b) -
!     & (1.0)*V(m,n,j,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOVV(j,n,e,b)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
!
c     1*( (T1+{me}*V{bjef}) * T2{aifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,j,f) = INTERMOVOV(m,b,j,f) +
     & T1(m,e)*V(b,e,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,b,j,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (V{bemf}*T1{am}) * T2{eifj})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(b,e,f,a) = INTERMVVVV(b,e,f,a) +
!      & (1.0/2.0)*V(b,m,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(b,e,f,a)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (V{aemf}*T1{bm}) * T2{eifj})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(a,e,f,b) = INTERMVVVV(a,e,f,b) -
!     & (1.0/2.0)*V(a,m,e,f)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(a,e,f,b)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{minj}*T1{am}) * T1{bn})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOOOV(i,n,j,a) = INTERMOOOV(i,n,j,a) +
!      & V(m,n,i,j)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOV(i,n,j,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
!
c     1*( (V{aebf}*T1{ei}) * T1{fj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  1:NOS))
        INTERMVVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVO(a,b,f,i) = INTERMVVVO(a,b,f,i) +
     & V(a,b,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVVO(a,b,f,i)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T1+{me}*T1{am}) * V{biej})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      & (1.0/2.0)*T1(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(e,a)*V(b,e,i,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
c     -1/2*( (T1+{me}*T1{ej}) * V{aibm})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,j) = INTERMOO(m,j) -
     &          (1.0/2.0)*T1(m,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOO(m,j)*V(a,b,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/2*( (T1+{me}*T1{ei}) * V{ajbm})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) +
     & (1.0/2.0)*T1(m,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOO(m,i)*V(a,b,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{me}*T1{bm}) * V{aiej})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) -
!     & (1.0/2.0)*T1(m,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(e,b)*V(a,e,i,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{bime}*T1{am}) * T1{ej})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(b,i,e,a) = INTERMVOVV(b,i,e,a) +
!      & V(b,m,i,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVOVV(b,i,e,a)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{aime}*T1{bm}) * T1{ej})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(a,i,e,b) = INTERMVOVV(a,i,e,b) -
!     & (1.0)*V(a,m,i,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVOVV(a,i,e,b)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{bjme}*T1{am}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(b,j,e,a) = INTERMVOVV(b,j,e,a) -
!     & (1.0)*V(b,m,j,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVOVV(b,j,e,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{ajme}*T1{bm}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(a,j,e,b) = INTERMVOVV(a,j,e,b) +
!      & V(a,m,j,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVOVV(a,j,e,b)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/24*( ( (T2+{menf}*T2{ambn}) * T2{eifo}) * F{oj})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) -
!     & (1.0/24.0)*T2(m,n,e,f)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,i,o) = INTERMVVOO(a,b,i,o) +
!      &  INTERMVVVV(e,f,a,b)*T2(i,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVOO(a,b,i,o)*FOCK(o,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!
c     -1/6*( ( (F{bg}*T2{fjgn}) * T2+{menf}) * T2{aiem})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(b,f,j,n) = INTERMVVOO(b,f,j,n) -
     &          (1.0/6.0)*FOCK(b,g)*T2(j,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(b,j,m,e) = INTERMVOOV(b,j,m,e) +
     &  INTERMVVOO(b,f,j,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOOV(b,j,m,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/24*( ( (T2+{menf}*T2{ambn}) * T2{ejfo}) * F{oi})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) +
!      & (1.0/24.0)*T2(m,n,e,f)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,j,o) = INTERMVVOO(a,b,j,o) +
!      &  INTERMVVVV(e,f,a,b)*T2(j,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVOO(a,b,j,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!
c     1/6*( ( (T2+{menf}*T2{ameo}) * T2{bjfn}) * F{oi})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,o) = INTERMOVVO(n,f,a,o) +
     & (1.0/6.0)*T2(m,n,e,f)*T2(m,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,o,b,j) = INTERMVOVO(a,o,b,j) +
     &  INTERMOVVO(n,f,a,o)*T2(j,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(a,o,b,j)*FOCK(o,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( ( (F{om}*T2{anbo}) * T2+{menf}) * T2{eifj})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVOV(m,a,n,b) = INTERMOVOV(m,a,n,b) +
!      & (1.0/4.0)*FOCK(o,m)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(a,b,e,f) = INTERMVVVV(a,b,e,f) +
!      &  INTERMOVOV(m,a,n,b)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(a,b,e,f)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!
c     1/2*( ( (T2+{menf}*T2{aien}) * F{om}) * T2{bjfo})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,f,a,i) = INTERMOVVO(m,f,a,i) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(f,a,i,o) = INTERMVVOO(f,a,i,o) +
     &  INTERMOVVO(m,f,a,i)*FOCK(o,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(f,a,i,o)*T2(j,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (T2+{menf}*T2{enfo}) * F{oi}) * T2{ajbm})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,o) = INTERMOO(m,o) +
!      & (1.0/12.0)*T2(m,n,e,f)*T2(n,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOO2(m,i) = INTERMOO2(m,i) +
!      &  INTERMOO(m,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(m,i)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
!
c     1/2*( ( (T2+{menf}*T2{bjfn}) * F{eg}) * T2{aigm})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,e,b,j) = INTERMOVVO(m,e,b,j) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(j,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,j,g) = INTERMOVOV(m,b,j,g) +
     &  INTERMOVVO(m,e,b,j)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,b,j,g)*T2(i,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( ( (T2+{menf}*F{om}) * T2{ejfo}) * T2{aibn})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(e,n,f,o) = INTERMVOVO(e,n,f,o) -
!     & (1.0/4.0)*T2(m,n,e,f)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      &  INTERMVOVO(e,n,f,o)*T2(j,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,j)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -1/2*( ( (F{om}*T2{bnfo}) * T2+{menf}) * T2{aiej})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,n,f) = INTERMOVOV(m,b,n,f) -
     &          (1.0/2.0)*FOCK(o,m)*T2(n,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,e) = INTERMVV(b,e) +
     &  INTERMOVOV(m,b,n,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +

     &  INTERMVV(b,e)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/6*( ( (F{ag}*T2{fjgn}) * T2+{menf}) * T2{biem})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,f,j,n) = INTERMVVOO(a,f,j,n) +
     & (1.0/6.0)*FOCK(a,g)*T2(j,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(a,j,m,e) = INTERMVOOV(a,j,m,e) +
     &  INTERMVVOO(a,f,j,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOOV(a,j,m,e)*T2(i,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     -1/6*( ( (T2+{menf}*T2{bmeo}) * T2{ajfn}) * F{oi})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,o) = INTERMOVVO(n,f,b,o) -
     &          (1.0/6.0)*T2(m,n,e,f)*T2(m,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(b,o,a,j) = INTERMVOVO(b,o,a,j) +
     &  INTERMOVVO(n,f,b,o)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(b,o,a,j)*FOCK(o,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1/2*( ( (T2+{menf}*T2{bien}) * F{om}) * T2{ajfo})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,f,b,i) = INTERMOVVO(m,f,b,i) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(i,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(f,b,i,o) = INTERMVVOO(f,b,i,o) +
     &  INTERMOVVO(m,f,b,i)*FOCK(o,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(f,b,i,o)*T2(j,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{amen}) * T2{bjfo}) * F{oi})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,a) = INTERMVV(f,a) -
     &          (1.0/12.0)*T2(m,n,e,f)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,j,o) = INTERMVVOO(a,b,j,o) +
     &  INTERMVV(f,a)*T2(j,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(a,b,j,o)*FOCK(o,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( ( (T2+{menf}*T2{ambn}) * F{eg}) * T2{figj})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) -
!     & (1.0/4.0)*T2(m,n,e,f)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV2 = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMVVVV2(f,a,b,g) = INTERMVVVV2(f,a,b,g) +
!      &  INTERMVVVV(e,f,a,b)*FOCK(e,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV2(f,a,b,g)*T2(i,j,f,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV2)
!         ENDIF
!
c     -1/6*( ( (F{ag}*T2{eigm}) * T2+{menf}) * T2{bjfn})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,e,i,m) = INTERMVVOO(a,e,i,m) -
     &          (1.0/6.0)*FOCK(a,g)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(a,i,n,f) = INTERMVOOV(a,i,n,f) +
     &  INTERMVVOO(a,e,i,m)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOOV(a,i,n,f)*T2(j,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     1/12*( ( (T2+{menf}*T2{bmfn}) * T2{eigj}) * F{ag})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,b) = INTERMVV(e,b) +
     & (1.0/12.0)*T2(m,n,e,f)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(b,i,g,j) = INTERMVOVO(b,i,g,j) +
     &  INTERMVV(e,b)*T2(i,j,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(b,i,g,j)*FOCK(a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( ( (F{om}*T2{anfo}) * T2+{menf}) * T2{biej})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,n,f) = INTERMOVOV(m,a,n,f) +
     & (1.0/2.0)*FOCK(o,m)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(a,e) = INTERMVV(a,e) +
     &  INTERMOVOV(m,a,n,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,e)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( ( (T2+{menf}*F{om}) * T2{eifo}) * T2{ajbn})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(e,n,f,o) = INTERMVOVO(e,n,f,o) +
!      & (1.0/4.0)*T2(m,n,e,f)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) +
!      &  INTERMVOVO(e,n,f,o)*T2(i,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,i)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( ( (T2+{menf}*T2{eifn}) * F{om}) * T2{ajbo})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) -
!     & (1.0/4.0)*T2(m,n,e,f)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOO2(i,o) = INTERMOO2(i,o) +
!      &  INTERMOO(m,i)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(i,o)*T2(j,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
!
c     1/12*( ( (T2+{menf}*T2{bmen}) * T2{ajfo}) * F{oi})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,b) = INTERMVV(f,b) +
     & (1.0/12.0)*T2(m,n,e,f)*T2(m,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,j,o) = INTERMVVOO(b,a,j,o) +
     &  INTERMVV(f,b)*T2(j,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(b,a,j,o)*FOCK(o,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/6*( ( (T2+{menf}*T2{ameo}) * T2{bifn}) * F{oj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,o) = INTERMOVVO(n,f,a,o) -
     &          (1.0/6.0)*T2(m,n,e,f)*T2(m,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,o,b,i) = INTERMVOVO(a,o,b,i) +
     &  INTERMOVVO(n,f,a,o)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(a,o,b,i)*FOCK(o,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/6*( ( (T2+{menf}*T2{bmeo}) * T2{aifn}) * F{oj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,o) = INTERMOVVO(n,f,b,o) +
     & (1.0/6.0)*T2(m,n,e,f)*T2(m,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(b,o,a,i) = INTERMVOVO(b,o,a,i) +
     &  INTERMOVVO(n,f,b,o)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(b,o,a,i)*FOCK(o,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/6*( ( (F{bg}*T2{fign}) * T2+{menf}) * T2{ajem})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(b,f,i,n) = INTERMVVOO(b,f,i,n) +
     & (1.0/6.0)*FOCK(b,g)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(b,i,m,e) = INTERMVOOV(b,i,m,e) +
     &  INTERMVVOO(b,f,i,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOOV(b,i,m,e)*T2(j,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/12*( ( (T2+{menf}*T2{eifm}) * T2{anbo}) * F{oj})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) -
!     & (1.0/12.0)*T2(m,n,e,f)*T2(i,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVVO(i,a,b,o) = INTERMOVVO(i,a,b,o) +
!      &  INTERMOO(n,i)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(i,a,b,o)*FOCK(o,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!
c     1/12*( ( (F{bg}*T2{fmgn}) * T2+{menf}) * T2{aiej})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(b,f,m,n) = INTERMVVOO(b,f,m,n) +
     & (1.0/12.0)*FOCK(b,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,e) = INTERMVV(b,e) +
     &  INTERMVVOO(b,f,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,e)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/2*( ( (T2+{menf}*T2{ajen}) * F{om}) * T2{bifo})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,f,a,j) = INTERMOVVO(m,f,a,j) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(j,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(f,a,j,o) = INTERMVVOO(f,a,j,o) +
     &  INTERMOVVO(m,f,a,j)*FOCK(o,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(f,a,j,o)*T2(i,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{eifm}) * T2{ajgn}) * F{bg})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,i) = INTERMOO(n,i) -
     &          (1.0/12.0)*T2(m,n,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(i,a,j,g) = INTERMOVOV(i,a,j,g) +
     &  INTERMOO(n,i)*T2(j,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(i,a,j,g)*FOCK(b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (T2+{menf}*T2{ejfm}) * T2{anbo}) * F{oi})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      & (1.0/12.0)*T2(m,n,e,f)*T2(j,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVVO(j,a,b,o) = INTERMOVVO(j,a,b,o) +
!      &  INTERMOO(n,j)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(j,a,b,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T2+{menf}*F{eg}) * T2{figm}) * T2{ajbn})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,f,g) = INTERMOOVV(m,n,f,g) +
!      & (1.0/2.0)*T2(m,n,e,f)*FOCK(e,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) +
!      &  INTERMOOVV(m,n,f,g)*T2(i,m,f,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,i)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     1/2*( ( (T2+{menf}*T2{bjen}) * F{om}) * T2{aifo})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,f,b,j) = INTERMOVVO(m,f,b,j) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(j,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(f,b,j,o) = INTERMVVOO(f,b,j,o) +
     &  INTERMOVVO(m,f,b,j)*FOCK(o,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(f,b,j,o)*T2(i,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{amfn}) * T2{eigj}) * F{bg})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) -
     &          (1.0/12.0)*T2(m,n,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(a,i,g,j) = INTERMVOVO(a,i,g,j) +
     &  INTERMVV(e,a)*T2(i,j,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVOVO(a,i,g,j)*FOCK(b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/12*( ( (T2+{menf}*T2{eifm}) * T2{bjgn}) * F{ag})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,i) = INTERMOO(n,i) +
     & (1.0/12.0)*T2(m,n,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(i,b,j,g) = INTERMOVOV(i,b,j,g) +
     &  INTERMOO(n,i)*T2(j,n,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS

        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(i,b,j,g)*FOCK(a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T2+{menf}*F{eg}) * T2{fjgm}) * T2{aibn})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,f,g) = INTERMOOVV(m,n,f,g) -
!     & (1.0/2.0)*T2(m,n,e,f)*FOCK(e,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      &  INTERMOOVV(m,n,f,g)*T2(j,m,f,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,j)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/12*( ( (T2+{menf}*T2{enfo}) * F{oj}) * T2{aibm})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,o) = INTERMOO(m,o) -
!     & (1.0/12.0)*T2(m,n,e,f)*T2(n,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOO2(m,j) = INTERMOO2(m,j) +
!      &  INTERMOO(m,o)*FOCK(o,j)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(m,j)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
!
c     -1/2*( ( (T2+{menf}*T2{bifn}) * F{eg}) * T2{ajgm})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,e,b,i) = INTERMOVVO(m,e,b,i) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,i,g) = INTERMOVOV(m,b,i,g) +
     &  INTERMOVVO(m,e,b,i)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,b,i,g)*T2(j,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1/12*( ( (T2+{menf}*T2{amen}) * T2{bifo}) * F{oj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,a) = INTERMVV(f,a) +
     & (1.0/12.0)*T2(m,n,e,f)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,i,o) = INTERMVVOO(a,b,i,o) +
     &  INTERMVV(f,a)*T2(i,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(a,b,i,o)*FOCK(o,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{bmen}) * T2{aifo}) * F{oj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,b) = INTERMVV(f,b) -
     &          (1.0/12.0)*T2(m,n,e,f)*T2(m,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,i,o) = INTERMVVOO(b,a,i,o) +
     &  INTERMVV(f,b)*T2(i,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOO(b,a,i,o)*FOCK(o,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/12*( ( (F{ag}*T2{fmgn}) * T2+{menf}) * T2{biej})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,f,m,n) = INTERMVVOO(a,f,m,n) -
     &          (1.0/12.0)*FOCK(a,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(a,e) = INTERMVV(a,e) +
     &  INTERMVVOO(a,f,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,e)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{ejfm}) * T2{bign}) * F{ag})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,j) = INTERMOO(n,j) -
     &          (1.0/12.0)*T2(m,n,e,f)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(j,b,i,g) = INTERMOVOV(j,b,i,g) +
     &  INTERMOO(n,j)*T2(i,n,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(j,b,i,g)*FOCK(a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/24*( ( (F{ag}*T2{bmgn}) * T2+{menf}) * T2{eifj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,m,n) = INTERMVVOO(a,b,m,n) -
     &          (1.0/24.0)*FOCK(a,g)*T2(m,n,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(a,b,e,f) = INTERMVVVV(a,b,e,f) +
     &  INTERMVVOO(a,b,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVVV(a,b,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/4*( ( (T2+{menf}*T2{bmfn}) * F{eg}) * T2{aigj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,b) = INTERMVV(e,b) -
     &          (1.0/4.0)*T2(m,n,e,f)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV2 = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVV2(b,g) = INTERMVV2(b,g) +
     &  INTERMVV(e,b)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV2(b,g)*T2(i,j,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV2)
        ENDIF
!
c     -1/2*( ( (T2+{menf}*T2{ajfn}) * F{eg}) * T2{bigm})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,e,a,j) = INTERMOVVO(m,e,a,j) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,j,g) = INTERMOVOV(m,a,j,g) +
     &  INTERMOVVO(m,e,a,j)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,a,j,g)*T2(i,m,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/4*( ( (F{eg}*T2{amgn}) * T2+{menf}) * T2{bifj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(e,a,m,n) = INTERMVVOO(e,a,m,n) -
     &          (1.0/4.0)*FOCK(e,g)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(a,f) = INTERMVV(a,f) +
     &  INTERMVVOO(e,a,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,f)*T2(i,j,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/12*( ( (T2+{menf}*T2{ejfm}) * T2{aign}) * F{bg})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,j) = INTERMOO(n,j) +
     & (1.0/12.0)*T2(m,n,e,f)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(j,a,i,g) = INTERMOVOV(j,a,i,g) +
     &  INTERMOO(n,j)*T2(i,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(j,a,i,g)*FOCK(b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( ( (T2+{menf}*T2{ejfn}) * F{om}) * T2{aibo})
! c      Commutator = 3, PT Order = 3
! !
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,j) = INTERMOO(m,j) +
!      & (1.0/4.0)*T2(m,n,e,f)*T2(j,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOO2(j,o) = INTERMOO2(j,o) +
!      &  INTERMOO(m,j)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(j,o)*T2(i,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
!
c     1/2*( ( (T2+{menf}*T2{aifn}) * F{eg}) * T2{bjgm})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,e,a,i) = INTERMOVVO(m,e,a,i) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,i,g) = INTERMOVOV(m,a,i,g) +
     &  INTERMOVVO(m,e,a,i)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,a,i,g)*T2(j,m,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1/4*( ( (F{eg}*T2{bmgn}) * T2+{menf}) * T2{aifj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(e,b,m,n) = INTERMVVOO(e,b,m,n) +
     & (1.0/4.0)*FOCK(e,g)*T2(m,n,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,f) = INTERMVV(b,f) +
     &  INTERMVVOO(e,b,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,f)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/24*( ( (F{bg}*T2{amgn}) * T2+{menf}) * T2{eifj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,m,n) = INTERMVVOO(b,a,m,n) +
     & (1.0/24.0)*FOCK(b,g)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(b,a,e,f) = INTERMVVVV(b,a,e,f) +
     &  INTERMVVOO(b,a,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVVV(b,a,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     1/4*( ( (T2+{menf}*T2{amfn}) * F{eg}) * T2{bigj})
c      Commutator = 3, PT Order = 3
!
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) +
     & (1.0/4.0)*T2(m,n,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVV2(a,g) = INTERMVV2(a,g) +
     &  INTERMVV(e,a)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV2(a,g)*T2(i,j,b,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV2)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menf}*T1{fm}) * F{ei}) * T2{ajbn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,n) = INTERMVO(e,n) -
!     & (1.0/3.0)*T2(m,n,e,f)*T1(m,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) +
!      &  INTERMVO(e,n)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,i)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     -1/6*( ( (F{nf}*T2{bmfn}) * T1+{me}) * T2{aiej})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(b,m) = INTERMVO(b,m) -
     &          (1.0/6.0)*FOCK(n,f)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(b,e) = INTERMVV(b,e) +
     &  INTERMVO(b,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,e)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (F{nf}*T2{ejfn}) * T1+{me}) * T2{aibm})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,j) = INTERMVO(e,j) -
!     & (1.0/6.0)*FOCK(n,f)*T2(j,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(j,m) = INTERMOO(j,m) +
!      &  INTERMVO(e,j)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(j,m)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     -2/3*( ( (T1+{me}*T2{bjen}) * F{nf}) * T2{aifm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,b,j,n) = INTERMOVOO(m,b,j,n) -
     &          (2.0/3.0)*T1(m,e)*T2(j,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))

        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,j,f) = INTERMOVOV(m,b,j,f) +
     &  INTERMOVOO(m,b,j,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,b,j,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*F{em}) * T1{fi}) * T2{ajbn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,f) = INTERMOV(n,f) +
!      & (1.0/3.0)*T2(m,n,e,f)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) +
!      &  INTERMOV(n,f)*T1(i,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,i)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (F{nf}*T2{eifn}) * T1+{me}) * T2{ajbm})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,i) = INTERMVO(e,i) +
!      & (1.0/6.0)*FOCK(n,f)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(i,m) = INTERMOO(i,m) +
!      &  INTERMVO(e,i)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(i,m)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     2/3*( ( (T1+{me}*T2{ajen}) * F{nf}) * T2{bifm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,a,j,n) = INTERMOVOO(m,a,j,n) +
     & (2.0/3.0)*T1(m,e)*T2(j,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,j,f) = INTERMOVOV(m,a,j,f) +
     &  INTERMOVOO(m,a,j,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,a,j,f)*T2(i,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T1+{me}*T2{aien}) * F{nf}) * T2{bjfm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,a,i,n) = INTERMOVOO(m,a,i,n) -
     &          (2.0/3.0)*T1(m,e)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,i,f) = INTERMOVOV(m,a,i,f) +
     &  INTERMOVOO(m,a,i,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,a,i,f)*T2(j,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{me}*T2{ejfm}) * F{nf}) * T2{aibn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(j,f) = INTERMOV(j,f) +
!      & (1.0/2.0)*T1(m,e)*T2(j,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(j,n) = INTERMOO(j,n) +
!      &  INTERMOV(j,f)*FOCK(n,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(j,n)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     2/3*( ( (T1+{me}*T2{eifj}) * F{nf}) * T2{ambn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO(m,i,f,j) = INTERMOOVO(m,i,f,j) +
!      & (2.0/3.0)*T1(m,e)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
!         INTERMOOOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOOO(m,i,j,n) = INTERMOOOO(m,i,j,n) +
!      &  INTERMOOVO(m,i,f,j)*FOCK(n,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOO(m,i,j,n)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*T1{fm}) * F{ej}) * T2{aibn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,n) = INTERMVO(e,n) +
!      & (1.0/3.0)*T2(m,n,e,f)*T1(m,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      &  INTERMVO(e,n)*FOCK(e,j)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,j)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     1/2*( ( (T1+{me}*T2{bmen}) * F{nf}) * T2{aifj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(b,n) = INTERMVO(b,n) +
     & (1.0/2.0)*T1(m,e)*T2(m,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,f) = INTERMVV(b,f) +
     &  INTERMVO(b,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,f)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{nf}*T2{amfn}) * T1+{me}) * T2{biej})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(a,m) = INTERMVO(a,m) +
     & (1.0/6.0)*FOCK(n,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(a,e) = INTERMVV(a,e) +
     &  INTERMVO(a,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,e)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T1+{me}*T2{eifm}) * F{nf}) * T2{ajbn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(i,f) = INTERMOV(i,f) -
!     & (1.0/2.0)*T1(m,e)*T2(i,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(i,n) = INTERMOO(i,n) +
!      &  INTERMOV(i,f)*FOCK(n,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(i,n)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     -1/2*( ( (T1+{me}*T2{amen}) * F{nf}) * T2{bifj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(a,n) = INTERMVO(a,n) -
     &          (1.0/2.0)*T1(m,e)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(a,f) = INTERMVV(a,f) +
     &  INTERMVO(a,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,f)*T2(i,j,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
c     2/3*( ( (T1+{me}*T2{bien}) * F{nf}) * T2{ajfm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,b,i,n) = INTERMOVOO(m,b,i,n) +
     & (2.0/3.0)*T1(m,e)*T2(i,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,i,f) = INTERMOVOV(m,b,i,f) +
     &  INTERMOVOO(m,b,i,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVOV(m,b,i,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (T2+{menf}*T1{en}) * F{am}) * T2{bifj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,f) = INTERMOV(m,f) -
     &          (1.0/3.0)*T2(m,n,e,f)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,a) = INTERMVV(f,a) +
     &  INTERMOV(m,f)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(f,a)*T2(i,j,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{menf}*F{ei}) * T1{fj}) * T2{ambn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO(m,n,f,i) = INTERMOOVO(m,n,f,i) +
!      & (1.0/6.0)*T2(m,n,e,f)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
!         INTERMOOOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOOO(m,n,i,j) = INTERMOOOO(m,n,i,j) +
!      &  INTERMOOVO(m,n,f,i)*T1(j,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOO(m,n,i,j)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menf}*F{em}) * T1{fj}) * T2{aibn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,f) = INTERMOV(n,f) -
!     & (1.0/3.0)*T2(m,n,e,f)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(n,j) = INTERMOO(n,j) +
!      &  INTERMOV(n,f)*T1(j,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(n,j)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     1/3*( ( (T2+{menf}*T1{en}) * F{bm}) * T2{aifj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,f) = INTERMOV(m,f) +
     & (1.0/3.0)*T2(m,n,e,f)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(f,b) = INTERMVV(f,b) +
     &  INTERMOV(m,f)*FOCK(b,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(f,b)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T2+{menf}*F{ej}) * T1{fi}) * T2{ambn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO(m,n,f,j) = INTERMOOVO(m,n,f,j) -
!     & (1.0/6.0)*T2(m,n,e,f)*FOCK(e,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
!         INTERMOOOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOOO(m,n,j,i) = INTERMOOOO(m,n,j,i) +
!      &  INTERMOOVO(m,n,f,j)*T1(i,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOOOO(m,n,j,i)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{menf}*F{am}) * T1{bn}) * T2{eifj})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) +
!      & (1.0/6.0)*T2(m,n,e,f)*FOCK(a,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,a,b) = INTERMVVVV(e,f,a,b) +
!      &  INTERMVOVV(e,n,f,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(e,f,a,b)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T2+{menf}*F{bm}) * T1{an}) * T2{eifj})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,b) = INTERMVOVV(e,n,f,b) -
!     & (1.0/6.0)*T2(m,n,e,f)*FOCK(b,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVV(e,f,b,a) = INTERMVVVV(e,f,b,a) +
!      &  INTERMVOVV(e,n,f,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVVV(e,f,b,a)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menf}*F{em}) * T1{bn}) * T2{aifj})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,f) = INTERMOV(n,f) -
!     & (1.0/3.0)*T2(m,n,e,f)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(f,b) = INTERMVV(f,b) +
!      &  INTERMOV(n,f)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(f,b)*T2(i,j,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*F{em}) * T1{an}) * T2{bifj})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,f) = INTERMOV(n,f) +
!      & (1.0/3.0)*T2(m,n,e,f)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(f,a) = INTERMVV(f,a) +
!      &  INTERMOV(n,f)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(f,a)*T2(i,j,b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
!
c     -1/3*( ( (T2+{menf}*F{bm}) * T1{ej}) * T2{aifn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,f,b) = INTERMVOVV(e,n,f,b) -
     &          (1.0/3.0)*T2(m,n,e,f)*FOCK(b,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,j) = INTERMOVVO(n,f,b,j) +
     &  INTERMVOVV(e,n,f,b)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,b,j)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
        ENDIF
!
c     1/3*( ( (T2+{menf}*F{am}) * T1{ej}) * T2{bifn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) +
     & (1.0/3.0)*T2(m,n,e,f)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,j) = INTERMOVVO(n,f,a,j) +
     &  INTERMVOVV(e,n,f,a)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,a,j)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*T1{am}) * F{ej}) * T2{bifn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) +
!      & (1.0/3.0)*T2(m,n,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVO(n,f,a,j) = INTERMOVVO(n,f,a,j) +
!      &  INTERMVOVV(e,n,f,a)*FOCK(e,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE

!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(n,f,a,j)*T2(i,n,b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!         ENDIF
!
c     1/3*( ( (T2+{menf}*F{bm}) * T1{ei}) * T2{ajfn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,f,b) = INTERMVOVV(e,n,f,b) +
     & (1.0/3.0)*T2(m,n,e,f)*FOCK(b,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,b,i) = INTERMOVVO(n,f,b,i) +
     &  INTERMVOVV(e,n,f,b)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,b,i)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menf}*T1{am}) * F{ei}) * T2{bjfn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) -
!     & (1.0/3.0)*T2(m,n,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVO(n,f,a,i) = INTERMOVVO(n,f,a,i) +
!      &  INTERMVOVV(e,n,f,a)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(n,f,a,i)*T2(j,n,b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!         ENDIF
!
c     -1/3*( ( (T2+{menf}*F{am}) * T1{ei}) * T2{bjfn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) -
     &          (1.0/3.0)*T2(m,n,e,f)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(n,f,a,i) = INTERMOVVO(n,f,a,i) +
     &  INTERMVOVV(e,n,f,a)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMOVVO(n,f,a,i)*T2(j,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*T1{bm}) * F{ei}) * T2{ajfn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,b) = INTERMVOVV(e,n,f,b) +
!      & (1.0/3.0)*T2(m,n,e,f)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVO(n,f,b,i) = INTERMOVVO(n,f,b,i) +
!      &  INTERMVOVV(e,n,f,b)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(n,f,b,i)*T2(j,n,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menf}*T1{bm}) * F{ej}) * T2{aifn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,b) = INTERMVOVV(e,n,f,b) -
!     & (1.0/3.0)*T2(m,n,e,f)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVO(n,f,b,j) = INTERMOVVO(n,f,b,j) +
!      &  INTERMVOVV(e,n,f,b)*FOCK(e,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(n,f,b,j)*T2(i,n,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{nm}*T1{bn}) * T1+{me}) * T2{aiej})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOV(m,b) = INTERMOV(m,b) +
!      & (1.0/2.0)*FOCK(n,m)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(b,e) = INTERMVV(b,e) +
!      &  INTERMOV(m,b)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(b,e)*T2(i,j,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*T1{am}) * T2{bien}) * F{nj})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/6.0)*T1(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,i,n) = INTERMVVOO(a,b,i,n) +
!      &  INTERMVV(e,a)*T2(i,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVOO(a,b,i,n)*FOCK(n,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{bm}) * T2{eifj}) * F{af})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/6.0)*T1(m,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(b,i,f,j) = INTERMVOVO(b,i,f,j) +
!      &  INTERMVV(e,b)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVOVO(b,i,f,j)*FOCK(a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*T1{am}) * T2{eifj}) * F{bf})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/6.0)*T1(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(a,i,f,j) = INTERMVOVO(a,i,f,j) +
!      &  INTERMVV(e,a)*T2(i,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVOVO(a,i,f,j)*FOCK(b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{ei}) * T2{ambn}) * F{nj})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) +
!      & (1.0/6.0)*T1(m,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVVO(i,a,b,n) = INTERMOVVO(i,a,b,n) +
!      &  INTERMOO(m,i)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(i,a,b,n)*FOCK(n,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!
c     -1/6*( ( (F{bf}*T1{fm}) * T1+{me}) * T2{aiej})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(b,m) = INTERMVO(b,m) -
     &          (1.0/6.0)*FOCK(b,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(b,e) = INTERMVV(b,e) +
     &  INTERMVO(b,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(b,e)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T1+{me}*T1{bm}) * F{ef}) * T2{aifj})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) -
!     & (1.0/2.0)*T1(m,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV2 = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV2(b,f) = INTERMVV2(b,f) +
!      &  INTERMVV(e,b)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV2(b,f)*T2(i,j,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV2)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (F{nm}*T1{an}) * T1+{me}) * T2{biej})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOV(m,a) = INTERMOV(m,a) -
!     & (1.0/2.0)*FOCK(n,m)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(a,e) = INTERMVV(a,e) +
!      &  INTERMOV(m,a)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV(a,e)*T2(i,j,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{am}) * T2{bjen}) * F{ni})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      & (1.0/6.0)*T1(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,j,n) = INTERMVVOO(a,b,j,n) +
!      &  INTERMVV(e,a)*T2(j,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVOO(a,b,j,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*T1{en}) * F{ni}) * T2{ajbm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,n) = INTERMOO(m,n) -
!     & (1.0/6.0)*T1(m,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOO2(m,i) = INTERMOO2(m,i) +
!      &  INTERMOO(m,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(m,i)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*T1{ej}) * T2{ambn}) * F{ni})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,j) = INTERMOO(m,j) -
!     & (1.0/6.0)*T1(m,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVVO(j,a,b,n) = INTERMOVVO(j,a,b,n) +
!      &  INTERMOO(m,j)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOVVO(j,a,b,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{me}*T1{am}) * F{ef}) * T2{bifj})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      & (1.0/2.0)*T1(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV2 = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV2(a,f) = INTERMVV2(a,f) +
!      &  INTERMVV(e,a)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVV2(a,f)*T2(i,j,b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV2)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T1+{me}*T1{ei}) * F{nm}) * T2{ajbn})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) -
!     & (1.0/2.0)*T1(m,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOO2(i,n) = INTERMOO2(i,n) +
!      &  INTERMOO(m,i)*FOCK(n,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(i,n)*T2(j,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{bm}) * T2{aien}) * F{nj})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/6.0)*T1(m,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,i,n) = INTERMVVOO(b,a,i,n) +
!      &  INTERMVV(e,b)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVOO(b,a,i,n)*FOCK(n,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*T1{bm}) * T2{ajen}) * F{ni})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) -
!     & (1.0/6.0)*T1(m,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,j,n) = INTERMVVOO(b,a,j,n) +
!      &  INTERMVV(e,b)*T2(j,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMVVOO(b,a,j,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{me}*T1{ej}) * F{nm}) * T2{aibn})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,j) = INTERMOO(m,j) +
!      & (1.0/2.0)*T1(m,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOO2(j,n) = INTERMOO2(j,n) +
!      &  INTERMOO(m,j)*FOCK(n,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(j,n)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
!
c     1/6*( ( (F{af}*T1{fm}) * T1+{me}) * T2{biej})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(a,m) = INTERMVO(a,m) +
     & (1.0/6.0)*FOCK(a,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(a,e) = INTERMVV(a,e) +
     &  INTERMVO(a,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVV(a,e)*T2(i,j,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{me}*F{ef}) * T1{fi}) * T2{ajbm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(m,f) = INTERMOV(m,f) +
!      & (1.0/2.0)*T1(m,e)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) +
!      &  INTERMOV(m,f)*T1(i,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(m,i)*T2(j,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{en}) * F{nj}) * T2{aibm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,n) = INTERMOO(m,n) +
!      & (1.0/6.0)*T1(m,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO2(1:NOS,1:NOS))
!         INTERMOO2 = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOO2(m,j) = INTERMOO2(m,j) +
!      &  INTERMOO(m,n)*FOCK(n,j)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO2(m,j)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO2)
!         ENDIF
!
c     1/6*( ( (F{af}*T2{bifm}) * T1+{me}) * T1{ej})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN

        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,i,m) = INTERMVVOO(a,b,i,m) +
     & (1.0/6.0)*FOCK(a,f)*T2(i,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(a,b,i,e) = INTERMVVOV(a,b,i,e) +
     &  INTERMVVOO(a,b,i,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOV(a,b,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
!
c     -1/6*( ( (F{bf}*T2{aifm}) * T1+{me}) * T1{ej})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,i,m) = INTERMVVOO(b,a,i,m) -
     &          (1.0/6.0)*FOCK(b,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) +
     &  INTERMVVOO(b,a,i,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOV(b,a,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
!
c     1/6*( ( (F{bf}*T2{ajfm}) * T1+{me}) * T1{ei})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,j,m) = INTERMVVOO(b,a,j,m) +
     & (1.0/6.0)*FOCK(b,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(b,a,j,e) = INTERMVVOV(b,a,j,e) +
     &  INTERMVVOO(b,a,j,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOV(b,a,j,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T1+{me}*F{ef}) * T1{fj}) * T2{aibm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(m,f) = INTERMOV(m,f) -
!     & (1.0/2.0)*T1(m,e)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,j) = INTERMOO(m,j) +
!      &  INTERMOV(m,f)*T1(j,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,j) = VT(a,b,i,j) +
!      &  INTERMOO(m,j)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -1/6*( ( (F{af}*T2{bjfm}) * T1+{me}) * T1{ei})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,j,m) = INTERMVVOO(a,b,j,m) -
     &          (1.0/6.0)*FOCK(a,f)*T2(j,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(a,b,j,e) = INTERMVVOV(a,b,j,e) +
     &  INTERMVVOO(a,b,j,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,j) = VT(a,b,i,j) +
     &  INTERMVVOV(a,b,j,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF

        DO  a  = NOS+1, NOS+NVS
        DO  b  = NOS+1, NOS+NVS
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        VT(i,j,a,b) = VT(a,b,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO

      END
