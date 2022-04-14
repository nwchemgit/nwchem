      SUBROUTINE PPPH(NOS,NVS,NVAS,NACTV,T1,T2,
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
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV2(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVV(:,:,:,:)

!
c     1*V{aibc}
c      Commutator = 0, PT Order = 1
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) + V(a,b,i,c)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*F{mc}*T2{bmai}
! c      Commutator = 1, PT Order = 1
! !
!         IF(NONHF)THEN
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) -
!     & (1.0)*FOCK(m,c)*T2(m,i,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*V{minc}*T2{bnam}
! c      Commutator = 1, PT Order = 2
! !
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      & (1.0/2.0)*V(m,n,i,c)*T2(n,m,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     1*V{bcme}*T2{aiem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     & V(b,m,c,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*V{acme}*T2{biem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) -
     &          (1.0)*V(a,m,c,e)*T2(i,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*T1+{mc}*V{aibm}
! c      Commutator = 1, PT Order = 3
! !
!         IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) -
!     & (1.0)*T1(m,c)*V(a,b,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     1*V{aebc}*T1{ei}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     & V(a,b,e,c)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*V{aimc}*T1{bm}
! c      Commutator = 1, PT Order = 3
! !
!         IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) -
!     & (1.0)*V(a,m,i,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*V{mibc}*T1{am}
! c      Commutator = 1, PT Order = 3
! !
!         IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) -
!     & (1.0)*V(m,b,i,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T2+{menc}*F{em}) * T2{bnai})
! c      Commutator = 2, PT Order = 2
! !
!         IF(NONHF)THEN
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) -
!     & (1.0/2.0)*T2(m,n,e,c)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(n,c)*T2(n,i,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( (T2+{menc}*F{ei}) * T2{bnam})
! c      Commutator = 2, PT Order = 2
! !
!         IF(NONHF)THEN
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVO(m,n,c,i) = INTERMOOVO(m,n,c,i) +
!      & (1.0/4.0)*T2(m,n,e,c)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOVO(m,n,c,i)*T2(n,m,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         ENDIF
!         ENDIF
!
c     1/2*( (T2+{menc}*F{bm}) * T2{aien})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,c,b) = INTERMVOVV(e,n,c,b) +
     & (1.0/2.0)*T2(m,n,e,c)*FOCK(b,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVOVV(e,n,c,b)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     -1/2*( (T2+{menc}*F{am}) * T2{bien})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,c,a) = INTERMVOVV(e,n,c,a) -
     &          (1.0/2.0)*T2(m,n,e,c)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVOVV(e,n,c,a)*T2(i,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T1+{mc}*T2{bman}) * F{ni})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVVO(c,b,a,n) = INTERMVVVO(c,b,a,n) +
!      & (1.0/2.0)*T1(m,c)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(c,b,a,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{me}*F{ec}) * T2{bmai})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,c) = INTERMOV(m,c) -
!     & (1.0/2.0)*T1(m,e)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(m,c)*T2(m,i,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{mc}*F{nm}) * T2{bnai})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVO(c,n) = INTERMVO(c,n) +
!      & T1(m,c)*FOCK(n,m)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,n)*T2(n,i,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (F{be}*T2{aiem}) * T1+{mc})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,i,m) = INTERMVVOO(b,a,i,m) -
!     & (1.0/2.0)*FOCK(b,e)*T2(i,m,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(b,a,i,m)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (F{ae}*T2{biem}) * T1+{mc})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,i,m) = INTERMVVOO(a,b,i,m) +
!      & (1.0/2.0)*FOCK(a,e)*T2(i,m,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(a,b,i,m)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T2+{menc}*V{oiem}) * T2{bnao})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVOO(n,c,o,i) = INTERMOVOO(n,c,o,i) +
!      & T2(m,n,e,c)*V(o,e,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVOO(n,c,o,i)*T2(n,o,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
!
c     1*( (T2+{menc}*T2{fien}) * V{bmaf})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,c,f,i) = INTERMOVVO(m,c,f,i) +
     & T2(m,n,e,c)*T2(i,n,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVO(m,c,f,i)*V(b,a,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1*( (T2+{menc}*T2{bneo}) * V{aiom})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,c,b,o) = INTERMOVVO(m,c,b,o) +
     & T2(m,n,e,c)*T2(n,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVO(m,c,b,o)*V(a,o,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1*( (T2+{menc}*T2{aneo}) * V{biom})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,c,a,o) = INTERMOVVO(m,c,a,o) -
     &          (1.0)*T2(m,n,e,c)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVO(m,c,a,o)*V(b,o,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1*( (T2+{menc}*V{bmef}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(n,c,b,f) = INTERMOVVV(n,c,b,f) +
     & T2(m,n,e,c)*V(b,e,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(n,c,b,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     -1*( (T2+{menc}*V{amef}) * T2{bifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(n,c,a,f) = INTERMOVVV(n,c,a,f) -
     &          (1.0)*T2(m,n,e,c)*V(a,e,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(n,c,a,f)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T2+{menc}*V{emon}) * T2{aibo})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(c,o) = INTERMVO(c,o) +
!      & (1.0/2.0)*T2(m,n,e,c)*V(e,o,m,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,o)*T2(i,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
c     1/2*( (T2+{menf}*V{bcem}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(n,f,b,c) = INTERMOVVV(n,f,b,c) +
     & (1.0/2.0)*T2(m,n,e,f)*V(b,e,c,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(n,f,b,c)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     -1/2*( (T2+{menf}*V{acem}) * T2{bifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(n,f,a,c) = INTERMOVVV(n,f,a,c) -
     &          (1.0/2.0)*T2(m,n,e,f)*V(a,e,c,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(n,f,a,c)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     -1/2*( (T2+{menc}*V{bmon}) * T2{aieo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  1:NOS))
        INTERMVVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVVO(e,c,b,o) = INTERMVVVO(e,c,b,o) -
     &          (1.0/2.0)*T2(m,n,e,c)*V(b,o,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVVO(e,c,b,o)*T2(i,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
!
c     1/2*( (T2+{menc}*V{amon}) * T2{bieo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  1:NOS))
        INTERMVVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVVO(e,c,a,o) = INTERMVVVO(e,c,a,o) +
     & (1.0/2.0)*T2(m,n,e,c)*V(a,o,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVVO(e,c,a,o)*T2(i,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
!
c     1/2*( (T2+{menc}*T2{bmfn}) * V{aief})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(e,c,b,f) = INTERMVVVV(e,c,b,f) +
     & (1.0/2.0)*T2(m,n,e,c)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVVV(e,c,b,f)*V(a,e,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/2*( (T2+{menc}*T2{amfn}) * V{bief})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE

        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(e,c,a,f) = INTERMVVVV(e,c,a,f) -
     &          (1.0/2.0)*T2(m,n,e,c)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVVV(e,c,a,f)*V(b,e,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/4*( (T2+{menc}*T2{emfn}) * V{aibf})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(c,f) = INTERMVV(c,f) -
     &          (1.0/4.0)*T2(m,n,e,c)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVV(c,f)*V(a,b,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{eifn}) * V{ambc})
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
        INTERMOO(m,i) = INTERMOO(m,i) -
     &          (1.0/4.0)*T2(m,n,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOO(m,i)*V(a,b,m,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( (T2+{menf}*V{emfc}) * T2{aibn})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) -
!     & (1.0/4.0)*T2(m,n,e,f)*V(e,f,m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(n,c)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
c     -1/4*( (T2+{menf}*T2{bmfn}) * V{aiec})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVV(e,b)*V(a,e,i,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/4*( (T2+{menf}*T2{amfn}) * V{biec})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVV(e,a)*V(b,e,i,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/8*( (T2+{menf}*T2{bman}) * V{fiec})
! c      Commutator = 2, PT Order = 3
! !
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
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
!         INTERMVVVV(e,f,b,a) = INTERMVVVV(e,f,b,a) +
!      & (1.0/8.0)*T2(m,n,e,f)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(e,f,b,a)*V(f,e,i,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{menc}*T1{em}) * T2{aibn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) -
!     & (1.0)*V(m,n,e,c)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(n,c)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
c     1/4*( (T2+{menc}*T1{ei}) * V{ambn})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,c,i) = INTERMOOVO(m,n,c,i) +
     & (1.0/4.0)*T2(m,n,e,c)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOOVO(m,n,c,i)*V(a,b,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( (T1+{mc}*T2{eifm}) * V{aebf})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOV(c,e,i,f) = INTERMVVOV(c,e,i,f) -
!     & (1.0/4.0)*T1(m,c)*T2(i,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOV(c,e,i,f)*V(a,b,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{me}*V{einc}) * T2{ambn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOV(m,i,n,c) = INTERMOOOV(m,i,n,c) +
!      & T1(m,e)*V(e,n,i,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOOV(m,i,n,c)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (T1+{me}*V{emnc}) * T2{aibn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) -
!     & (1.0)*T1(m,e)*V(e,n,m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(n,c)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{menc}*T1{bm}) * T2{aien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,b) = INTERMVOVV(e,n,c,b) +
!      & V(m,n,e,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVOVV(e,n,c,b)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{menc}*T1{am}) * T2{bien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,a) = INTERMVOVV(e,n,c,a) -
!     & (1.0)*V(m,n,e,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVOVV(e,n,c,a)*T2(i,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!
c     1*( (T1+{me}*V{bcef}) * T2{aifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(m,b,c,f) = INTERMOVVV(m,b,c,f) +
     & T1(m,e)*V(b,e,c,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(m,b,c,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     -1*( (T1+{me}*V{acef}) * T2{bifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(m,a,c,f) = INTERMOVVV(m,a,c,f) -
     &          (1.0)*T1(m,e)*V(a,e,c,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(m,a,c,f)*T2(i,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     -1*( (T1+{me}*T2{aien}) * V{bcnm})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVOO(m,a,i,n)*V(b,n,c,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1*( (T1+{me}*T2{bien}) * V{acnm})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVOO(m,b,i,n)*V(a,n,c,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (T1+{mc}*V{bmne}) * T2{aien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(c,b,n,e) = INTERMVVOV(c,b,n,e) -
!     & (1.0)*T1(m,c)*V(b,n,m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOV(c,b,n,e)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{mc}*V{amne}) * T2{bien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(c,a,n,e) = INTERMVVOV(c,a,n,e) +
!      & T1(m,c)*V(a,n,m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOV(c,a,n,e)*T2(i,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (V{menc}*T1{ei}) * T2{ambn})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVO(m,n,c,i) = INTERMOOVO(m,n,c,i) +
!      & (1.0/2.0)*V(m,n,e,c)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOVO(m,n,c,i)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         ENDIF
!
c     -1/2*( (T1+{me}*T2{eifm}) * V{afbc})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) -
     &          (1.0/2.0)*T1(m,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOV(i,f)*V(a,b,f,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T2+{menc}*T1{en}) * V{aibm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(m,c) = INTERMOV(m,c) +
     & (1.0/2.0)*T2(m,n,e,c)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOV(m,c)*V(a,b,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{mc}*V{niom}) * T2{anbo})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOO(NOS+1:NOS+NVS,1:NOS,1:NOS,1:NOS))
!         INTERMVOOO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOOO(c,n,i,o) = INTERMVOOO(c,n,i,o) -
!     & (1.0/2.0)*T1(m,c)*V(n,o,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVOOO(c,n,i,o)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T2+{menc}*T1{bn}) * V{aiem})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,c,b) = INTERMOVVV(m,e,c,b) -
!     & (1.0/2.0)*T2(m,n,e,c)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVVV(m,e,c,b)*V(a,e,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T2+{menc}*T1{an}) * V{biem})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,c,a) = INTERMOVVV(m,e,c,a) +
!      & (1.0/2.0)*T2(m,n,e,c)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVVV(m,e,c,a)*V(b,e,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{mc}*T2{bmen}) * V{aine})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(c,b,e,n) = INTERMVVVO(c,b,e,n) -
!     & (1.0/2.0)*T1(m,c)*T2(m,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(c,b,e,n)*V(a,n,i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T1+{mc}*T2{amen}) * V{bine})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS

!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(c,a,e,n) = INTERMVVVO(c,a,e,n) +
!      & (1.0/2.0)*T1(m,c)*T2(m,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(c,a,e,n)*V(b,n,i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
c     -1/2*( (T1+{me}*T2{bnem}) * V{ainc})
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
        INTERMVO(b,n) = INTERMVO(b,n) -
     &          (1.0/2.0)*T1(m,e)*T2(n,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVO(b,n)*V(a,n,i,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1/2*( (T1+{me}*T2{anem}) * V{binc})
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
        INTERMVO(a,n) = INTERMVO(a,n) +
     & (1.0/2.0)*T1(m,e)*T2(n,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVO(a,n)*V(b,n,i,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{minc}*T1{am}) * T1{bn})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOOVV(i,n,c,a) = INTERMOOVV(i,n,c,a) +
!      & V(m,n,i,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOVV(i,n,c,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{mc}*T1{em}) * V{aibe})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(c,e) = INTERMVV(c,e) -
!     & (1.0/2.0)*T1(m,c)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVV(c,e)*V(a,b,i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{mc}*V{ambe}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVV(c,a,b,e) = INTERMVVVV(c,a,b,e) +
!      & T1(m,c)*V(a,b,m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(c,a,b,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!
c     1/2*( (T1+{me}*T1{ei}) * V{acbm})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOO(m,i)*V(a,b,c,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{mc}*V{ainm}) * T1{bn})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOO(c,a,i,n) = INTERMVVOO(c,a,i,n) +
!      & T1(m,c)*V(a,n,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(c,a,i,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (T1+{mc}*V{binm}) * T1{an})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOO(c,b,i,n) = INTERMVVOO(c,b,i,n) -
!     & (1.0)*T1(m,c)*V(b,n,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(c,b,i,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{bcme}*T1{am}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVV(b,c,e,a) = INTERMVVVV(b,c,e,a) -
!     & (1.0)*V(b,m,c,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(b,c,e,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{acme}*T1{bm}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVV(a,c,e,b) = INTERMVVVV(a,c,e,b) +
!      & V(a,m,c,e)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(a,c,e,b)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{me}*T1{bm}) * V{aiec})
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVV(e,b)*V(a,e,i,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T1+{me}*T1{am}) * V{biec})
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVV(e,a)*V(b,e,i,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -2/3*( ( (F{of}*T2{fien}) * T2+{menc}) * T2{bmao})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO(o,i,e,n) = INTERMOOVO(o,i,e,n) -
!     & (2.0/3.0)*FOCK(o,f)*T2(i,n,f,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOV(o,i,m,c) = INTERMOOOV(o,i,m,c) +
!      &  INTERMOOVO(o,i,e,n)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOOV(o,i,m,c)*T2(m,o,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( ( (T2+{menc}*T2{emfn}) * F{of}) * T2{aibo})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(c,f) = INTERMVV(c,f) +
!      & (1.0/4.0)*T2(m,n,e,c)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(c,o) = INTERMVO(c,o) +
!      &  INTERMVV(c,f)*FOCK(o,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,o)*T2(i,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/24*( ( (T2+{menf}*T2{bman}) * T2{eifo}) * F{oc})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
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
!         INTERMVVVV(e,f,b,a) = INTERMVVVV(e,f,b,a) +
!      & (1.0/24.0)*T2(m,n,e,f)*T2(m,n,b,a)
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
!         INTERMVVOO(b,a,i,o) = INTERMVVOO(b,a,i,o) +
!      &  INTERMVVVV(e,f,b,a)*T2(i,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(b,a,i,o)*FOCK(o,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (F{of}*T2{enfo}) * T2+{menc}) * T2{aibm})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,n) = INTERMVO(e,n) +
!      & (1.0/6.0)*FOCK(o,f)*T2(n,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,c) = INTERMOV(m,c) +
!      &  INTERMVO(e,n)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(m,c)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/12*( ( (F{oc}*T2{enfo}) * T2+{menf}) * T2{aibm})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOV(c,e,n,f) = INTERMVVOV(c,e,n,f) -
!     & (1.0/12.0)*FOCK(o,c)*T2(n,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(c,m) = INTERMVO(c,m) +
!      &  INTERMVVOV(c,e,n,f)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,m)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (F{of}*T2{eifo}) * T2+{mcne}) * T2{bman})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,i) = INTERMVO(e,i) +
!      & (1.0/12.0)*FOCK(o,f)*T2(i,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVO(i,m,c,n) = INTERMOOVO(i,m,c,n) +
!      &  INTERMVO(e,i)*T2(m,n,c,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOVO(i,m,c,n)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (T2+{menf}*T2{eifn}) * T2{boam}) * F{oc})
! c      Commutator = 3, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/4, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) +
!      & (1.0/12.0)*T2(m,n,e,f)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVOV(i,b,o,a) = INTERMOVOV(i,b,o,a) +
!      &  INTERMOO(m,i)*T2(o,m,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVOV(i,b,o,a)*FOCK(o,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!         ENDIF
!
c     2/3*( ( (T2+{menc}*T2{bneo}) * F{of}) * T2{aifm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,c,b,o) = INTERMOVVO(m,c,b,o) +
     & (2.0/3.0)*T2(m,n,e,c)*T2(n,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(m,c,b,f) = INTERMOVVV(m,c,b,f) +
     &  INTERMOVVO(m,c,b,o)*FOCK(o,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(m,c,b,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T2+{menc}*T2{aneo}) * F{of}) * T2{bifm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,c,a,o) = INTERMOVVO(m,c,a,o) -
     &          (2.0/3.0)*T2(m,n,e,c)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(m,c,a,f) = INTERMOVVV(m,c,a,f) +
     &  INTERMOVVO(m,c,a,o)*FOCK(o,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(m,c,a,f)*T2(i,m,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (F{of}*T2{bmfn}) * T2+{menc}) * T2{aieo})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(o,b,m,n) = INTERMOVOO(o,b,m,n) -
     &          (1.0/3.0)*FOCK(o,f)*T2(m,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVV(o,b,e,c) = INTERMOVVV(o,b,e,c) +
     &  INTERMOVOO(o,b,m,n)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(o,b,e,c)*T2(i,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     1/3*( ( (F{of}*T2{amfn}) * T2+{menc}) * T2{bieo})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(o,a,m,n) = INTERMOVOO(o,a,m,n) +
     & (1.0/3.0)*FOCK(o,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVV(o,a,e,c) = INTERMOVVV(o,a,e,c) +
     &  INTERMOVOO(o,a,m,n)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMOVVV(o,a,e,c)*T2(i,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{of}*T2{bmfo}) * T2+{menc}) * T2{aien})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(b,m) = INTERMVO(b,m) +
     & (1.0/6.0)*FOCK(o,f)*T2(m,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0

        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(b,e,n,c) = INTERMVVOV(b,e,n,c) +
     &  INTERMVO(b,m)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV(b,e,n,c)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{of}*T2{amfo}) * T2+{menc}) * T2{bien})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(a,m) = INTERMVO(a,m) -
     &          (1.0/6.0)*FOCK(o,f)*T2(m,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(a,e,n,c) = INTERMVVOV(a,e,n,c) +
     &  INTERMVO(a,m)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV(a,e,n,c)*T2(i,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{oc}*T2{bnfo}) * T2+{menf}) * T2{aiem})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(c,b,n,f) = INTERMVVOV(c,b,n,f) +
     & (1.0/6.0)*FOCK(o,c)*T2(n,o,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV2 = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV2(c,b,m,e) = INTERMVVOV2(c,b,m,e) +
     &  INTERMVVOV(c,b,n,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV2(c,b,m,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV2)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{oc}*T2{anfo}) * T2+{menf}) * T2{biem})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(c,a,n,f) = INTERMVVOV(c,a,n,f) -
     &          (1.0/6.0)*FOCK(o,c)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV2(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV2(c,a,m,e) = INTERMVVOV2(c,a,m,e) +
     &  INTERMVVOV(c,a,n,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV2(c,a,m,e)*T2(i,m,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV2)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{menf}*T2{bmfn}) * T2{aieo}) * F{oc})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,i,o) = INTERMVVOO(b,a,i,o) +
     &  INTERMVV(e,b)*T2(i,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOO(b,a,i,o)*FOCK(o,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{amfn}) * T2{bieo}) * F{oc})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,i,o) = INTERMVVOO(a,b,i,o) +
     &  INTERMVV(e,a)*T2(i,o,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOO(a,b,i,o)*FOCK(o,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/12*( ( (F{oi}*T1{eo}) * T2+{menc}) * T2{bnam})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(i,e) = INTERMOV(i,e) -
!     & (1.0/12.0)*FOCK(o,i)*T1(o,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOV(i,m,n,c) = INTERMOOOV(i,m,n,c) +
!      &  INTERMOV(i,e)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOOV(i,m,n,c)*T2(n,m,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{om}*T1{eo}) * T2+{menc}) * T2{aibn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,e) = INTERMOV(m,e) +
!      & (1.0/2.0)*FOCK(o,m)*T1(o,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOV2(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV2 = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV2(n,c) = INTERMOV2(n,c) +
!      &  INTERMOV(m,e)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV2(n,c)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV2)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (F{ef}*T1{fm}) * T2+{menc}) * T2{aibn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,m) = INTERMVO(e,m) -
!     & (1.0/2.0)*FOCK(e,f)*T1(m,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) +
!      &  INTERMVO(e,m)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(n,c)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{menc}*T1{em}) * T2{bnao}) * F{oi})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) +
!      & (1.0/6.0)*T2(m,n,e,c)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVVVO(c,b,a,o) = INTERMVVVO(c,b,a,o) +
!      &  INTERMOV(n,c)*T2(n,o,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(c,b,a,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T2+{menc}*T1{ei}) * F{om}) * T2{bnao})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVO(m,n,c,i) = INTERMOOVO(m,n,c,i) -
!     & (1.0/2.0)*T2(m,n,e,c)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVOO(n,c,i,o) = INTERMOVOO(n,c,i,o) +
!      &  INTERMOOVO(m,n,c,i)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVOO(n,c,i,o)*T2(n,o,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T2+{menc}*T1{en}) * F{om}) * T2{boai})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,c) = INTERMOV(m,c) -
!     & (1.0/2.0)*T2(m,n,e,c)*T1(n,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVO(c,o) = INTERMVO(c,o) +
!      &  INTERMOV(m,c)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,o)*T2(o,i,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*T1{fm}) * F{ec}) * T2{bnai})
! c      Commutator = 3, PT Order = 4
! !
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
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) +
!      &  INTERMVO(e,n)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOV(n,c)*T2(n,i,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{menf}*F{ec}) * T1{fi}) * T2{bman})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,f,c) = INTERMOOVV(m,n,f,c) +
!      & (1.0/6.0)*T2(m,n,e,f)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO(m,n,c,i) = INTERMOOVO(m,n,c,i) +
!      &  INTERMOOVV(m,n,f,c)*T1(i,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOVO(m,n,c,i)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( ( (F{ef}*T1{fi}) * T2+{mcne}) * T2{bman})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,i) = INTERMVO(e,i) +
!      & (1.0/4.0)*FOCK(e,f)*T1(i,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVO(i,m,c,n) = INTERMOOVO(i,m,c,n) +
!      &  INTERMVO(e,i)*T2(m,n,c,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOOVO(i,m,c,n)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO)
!         ENDIF
!
c     1/12*( ( (F{bf}*T2{amfn}) * T2+{menc}) * T1{ei})
c      Commutator = 3, PT Order = 4
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
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(b,a,m,n) = INTERMVVOO(b,a,m,n) +
     & (1.0/12.0)*FOCK(b,f)*T2(m,n,a,f)
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVVV(b,a,e,c) = INTERMVVVV(b,a,e,c) +
     &  INTERMVVOO(b,a,m,n)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVVV(b,a,e,c)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/12*( ( (F{af}*T2{bmfn}) * T2+{menc}) * T1{ei})
c      Commutator = 3, PT Order = 4
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
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,b,m,n) = INTERMVVOO(a,b,m,n) -
     &          (1.0/12.0)*FOCK(a,f)*T2(m,n,b,f)
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVVV(a,b,e,c) = INTERMVVVV(a,b,e,c) +
     &  INTERMVVOO(a,b,m,n)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVVV(a,b,e,c)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (F{om}*T1{bo}) * T2+{menc}) * T2{aien})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOV(m,b) = INTERMOV(m,b) -
!     & (1.0/2.0)*FOCK(o,m)*T1(o,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,e,n,c) = INTERMVVOV(b,e,n,c) +
!      &  INTERMOV(m,b)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOV(b,e,n,c)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{om}*T1{ao}) * T2+{menc}) * T2{bien})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOV(m,a) = INTERMOV(m,a) +
!      & (1.0/2.0)*FOCK(o,m)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(a,e,n,c) = INTERMVVOV(a,e,n,c) +
!      &  INTERMOV(m,a)*T2(m,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOV(a,e,n,c)*T2(i,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T2+{menc}*T1{bm}) * F{ef}) * T2{aifn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,b) = INTERMVOVV(e,n,c,b) +
!      & (1.0/2.0)*T2(m,n,e,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE

!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVV(n,c,b,f) = INTERMOVVV(n,c,b,f) +
!      &  INTERMVOVV(e,n,c,b)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVVV(n,c,b,f)*T2(i,n,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T2+{menc}*T1{am}) * F{ef}) * T2{bifn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,a) = INTERMVOVV(e,n,c,a) -
!     & (1.0/2.0)*T2(m,n,e,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVV(n,c,a,f) = INTERMOVVV(n,c,a,f) +
!      &  INTERMVOVV(e,n,c,a)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVVV(n,c,a,f)*T2(i,n,b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T2+{menc}*T1{bn}) * F{om}) * T2{aieo})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,c,b) = INTERMOVVV(m,e,c,b) +
!      & (1.0/2.0)*T2(m,n,e,c)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(e,c,b,o) = INTERMVVVO(e,c,b,o) +
!      &  INTERMOVVV(m,e,c,b)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(e,c,b,o)*T2(i,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T2+{menc}*T1{an}) * F{om}) * T2{bieo})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,c,a) = INTERMOVVV(m,e,c,a) -
!     & (1.0/2.0)*T2(m,n,e,c)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(e,c,a,o) = INTERMVVVO(e,c,a,o) +
!      &  INTERMOVVV(m,e,c,a)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(e,c,a,o)*T2(i,o,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
c     -1/6*( ( (T2+{menc}*T1{em}) * T2{aifn}) * F{bf})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(n,c) = INTERMOV(n,c) -
     &          (1.0/6.0)*T2(m,n,e,c)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(c,a,i,f) = INTERMVVOV(c,a,i,f) +
     &  INTERMOV(n,c)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV(c,a,i,f)*FOCK(b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
!
c     1/6*( ( (T2+{menc}*T1{em}) * T2{bifn}) * F{af})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(n,c) = INTERMOV(n,c) +
     & (1.0/6.0)*T2(m,n,e,c)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(c,b,i,f) = INTERMVVOV(c,b,i,f) +
     &  INTERMOV(n,c)*T2(i,n,b,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV(c,b,i,f)*FOCK(a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{menc}*T1{bm}) * T2{aneo}) * F{oi})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,b) = INTERMVOVV(e,n,c,b) +
!      & (1.0/6.0)*T2(m,n,e,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(c,b,a,o) = INTERMVVVO(c,b,a,o) +
!      &  INTERMVOVV(e,n,c,b)*T2(n,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(c,b,a,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T2+{menc}*T1{am}) * T2{bneo}) * F{oi})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,a) = INTERMVOVV(e,n,c,a) -
!     & (1.0/6.0)*T2(m,n,e,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(c,a,b,o) = INTERMVVVO(c,a,b,o) +
!      &  INTERMVOVV(e,n,c,a)*T2(n,o,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVO(c,a,b,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menf}*F{ec}) * T1{bm}) * T2{aifn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,f,c) = INTERMOOVV(m,n,f,c) -
!     & (1.0/3.0)*T2(m,n,e,f)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVV(n,f,c,b) = INTERMOVVV(n,f,c,b) +
!      &  INTERMOOVV(m,n,f,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVVV(n,f,c,b)*T2(i,n,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menf}*F{ec}) * T1{am}) * T2{bifn})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,f,c) = INTERMOOVV(m,n,f,c) +
!      & (1.0/3.0)*T2(m,n,e,f)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVV(n,f,c,a) = INTERMOVVV(n,f,c,a) +
!      &  INTERMOOVV(m,n,f,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVVV(n,f,c,a)*T2(i,n,b,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (F{af}*T2{fien}) * T2+{mcne}) * T1{bm})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(a,i,e,n) = INTERMVOVO(a,i,e,n) -
!     & (1.0/6.0)*FOCK(a,f)*T2(i,n,f,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(a,i,m,c) = INTERMVOOV(a,i,m,c) +
!      &  INTERMVOVO(a,i,e,n)*T2(m,n,c,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVOOV(a,i,m,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (F{bf}*T2{fien}) * T2+{mcne}) * T1{am})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(b,i,e,n) = INTERMVOVO(b,i,e,n) +
!      & (1.0/6.0)*FOCK(b,f)*T2(i,n,f,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(b,i,m,c) = INTERMVOOV(b,i,m,c) +
!      &  INTERMVOVO(b,i,e,n)*T2(m,n,c,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVOOV(b,i,m,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
c     1/6*( ( (F{bf}*T1{fm}) * T2+{menc}) * T2{aien})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(b,m) = INTERMVO(b,m) +
     & (1.0/6.0)*FOCK(b,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(b,e,n,c) = INTERMVVOV(b,e,n,c) +
     &  INTERMVO(b,m)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV(b,e,n,c)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
!
c     -1/6*( ( (F{af}*T1{fm}) * T2+{menc}) * T2{bien})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(a,m) = INTERMVO(a,m) -
     &          (1.0/6.0)*FOCK(a,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(a,e,n,c) = INTERMVVOV(a,e,n,c) +
     &  INTERMVO(a,m)*T2(m,n,e,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,i,c) = VT(a,b,i,c) +
     &  INTERMVVOV(a,e,n,c)*T2(i,n,b,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menc}*T1{bn}) * T1{am}) * F{ei})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,c,b) = INTERMOVVV(m,e,c,b) +
!      & (1.0/3.0)*T2(m,n,e,c)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVV(e,c,b,a) = INTERMVVVV(e,c,b,a) +
!      &  INTERMOVVV(m,e,c,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(e,c,b,a)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{menc}*F{bm}) * T1{an}) * T1{ei})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,b) = INTERMVOVV(e,n,c,b) -
!     & (1.0/3.0)*T2(m,n,e,c)*FOCK(b,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVV(e,c,b,a) = INTERMVVVV(e,c,b,a) +
!      &  INTERMVOVV(e,n,c,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(e,c,b,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{menc}*F{am}) * T1{bn}) * T1{ei})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,c,a) = INTERMVOVV(e,n,c,a) +
!      & (1.0/3.0)*T2(m,n,e,c)*FOCK(a,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVV(e,c,a,b) = INTERMVVVV(e,c,a,b) +
!      &  INTERMVOVV(e,n,c,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVVV(e,c,a,b)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{ne}*T1{ei}) * T2{bman}) * T1+{mc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) +
!      & (1.0/2.0)*FOCK(n,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVOV(i,b,m,a) = INTERMOVOV(i,b,m,a) +
!      &  INTERMOO(n,i)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVOV(i,b,m,a)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{mc}*F{em}) * T1+{ne}) * T2{aibn})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(c,e) = INTERMVV(c,e) +
!      & (1.0/2.0)*T1(m,c)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(c,n) = INTERMVO(c,n) +
!      &  INTERMVV(c,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,n)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{mc}*T1{em}) * F{ne}) * T2{aibn})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(c,e) = INTERMVV(c,e) +
!      & (1.0/2.0)*T1(m,c)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(c,n) = INTERMVO(c,n) +
!      &  INTERMVV(c,e)*FOCK(n,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE

!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,n)*T2(i,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{ne}*T1{bn}) * T2{aiem}) * T1+{mc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/2.0)*FOCK(n,e)*T1(n,b)
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
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,i,m) = INTERMVVOO(b,a,i,m) +
!      &  INTERMVV(e,b)*T2(i,m,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(b,a,i,m)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (F{ne}*T1{an}) * T2{biem}) * T1+{mc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/2.0)*FOCK(n,e)*T1(n,a)
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
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,i,m) = INTERMVVOO(a,b,i,m) +
!      &  INTERMVV(e,a)*T2(i,m,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(a,b,i,m)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{ne}*F{ei}) * T2{bman}) * T1+{mc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(n,i) = INTERMOO(n,i) +
!      & (1.0/6.0)*T1(n,e)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVOV(i,b,m,a) = INTERMOVOV(i,b,m,a) +
!      &  INTERMOO(n,i)*T2(m,n,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVOV(i,b,m,a)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{ei}) * T2{bnam}) * F{nc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) +
!      & (1.0/6.0)*T1(m,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVOV(i,b,n,a) = INTERMOVOV(i,b,n,a) +
!      &  INTERMOO(m,i)*T2(n,m,b,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMOVOV(i,b,n,a)*FOCK(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (F{nc}*T1{en}) * T1+{me}) * T2{aibm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(c,e) = INTERMVV(c,e) +
!      & (1.0/6.0)*FOCK(n,c)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(c,m) = INTERMVO(c,m) +
!      &  INTERMVV(c,e)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVO(c,m)*T2(i,m,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{bm}) * T2{aien}) * F{nc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(b,a,i,n)*FOCK(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*T1{am}) * T2{bien}) * F{nc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(a,b,i,n)*FOCK(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*F{bm}) * T2{aien}) * T1+{nc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/6.0)*T1(m,e)*FOCK(b,m)
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(b,a,i,n)*T1(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T1+{me}*F{am}) * T2{bien}) * T1+{nc})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/6.0)*T1(m,e)*FOCK(a,m)
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,i,c) = VT(a,b,i,c) +
!      &  INTERMVVOO(a,b,i,n)*T1(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF

        DO  a  = NOS+1, NOS+NVS
        DO  b  = NOS+1, NOS+NVS
        DO  c  = NOS+1, NOS+NVS
        DO  i  = 1, NOS
        VT(b,a,c,i) = VT(a,b,i,c)
        VT(i,c,a,b) = VT(a,b,i,c)
        VT(c,i,b,a) = VT(a,b,i,c)
        ENDDO; ENDDO; ENDDO; ENDDO

      END
