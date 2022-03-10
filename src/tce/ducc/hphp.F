      SUBROUTINE HPHP(NOS,NVS,NVAS,NACTV,T1,T2,
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
      ! DOUBLE PRECISION, ALLOCATABLE :: INTERMOV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO2(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVV(:,:,:,:)
      ! DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV(:,:,:,:)
      ! DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV2(:,:,:,:)
      ! DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVV(:,:,:,:)

!
c     1*V{aijb}
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
        VT(a,j,i,b) = VT(a,j,i,b) + V(a,j,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{mejb}*T2{aiem}
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     & V(m,j,e,b)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*T2+{mejb}*V{aiem}
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     & T2(m,j,e,b)*V(a,e,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{aejb}*T1{ei}
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     & V(a,j,e,b)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*T1+{mb}*V{aijm}
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
!         VT(a,j,i,b) = VT(a,j,i,b) -
!     & (1.0)*T1(m,b)*V(a,j,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     1*T1+{je}*V{aieb}
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     & T1(j,e)*V(a,e,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*V{mijb}*T1{am}
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
!         VT(a,j,i,b) = VT(a,j,i,b) -
!     & (1.0)*V(m,j,i,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     -1*( (T2+{mejb}*F{nm}) * T2{aien})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVOVO(e,j,b,n) = INTERMVOVO(e,j,b,n) -
     &          (1.0)*T2(m,j,e,b)*FOCK(n,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(e,j,b,n)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1*( (T2+{mejb}*F{ef}) * T2{aifm})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,b,f) = INTERMOOVV(m,j,b,f) +
     & T2(m,j,e,b)*FOCK(e,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVV(m,j,b,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/2*( (T2+{mejb}*T2{amen}) * F{ni})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(j,b,a,n) = INTERMOVVO(j,b,a,n) +
     & (1.0/2.0)*T2(m,j,e,b)*T2(m,n,a,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVO(j,b,a,n)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1/2*( (T2+{menb}*T2{aien}) * F{jm})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,b,a,i) = INTERMOVVO(m,b,a,i) +
     & (1.0/2.0)*T2(m,n,e,b)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVO(m,b,a,i)*FOCK(j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1/2*( (T2+{mejf}*F{eb}) * T2{aifm})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,f,b) = INTERMOOVV(m,j,f,b) -
     &          (1.0/2.0)*T2(m,j,e,f)*FOCK(e,b)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVV(m,j,f,b)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/2*( (F{af}*T2{fiem}) * T2+{mejb})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,i,e,m) = INTERMVOVO(a,i,e,m) +
     & (1.0/2.0)*FOCK(a,f)*T2(i,m,f,e)
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
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(a,i,e,m)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1*( (T2+{menb}*V{emjf}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(n,b,j,f) = INTERMOVOV(n,b,j,f) -
     &          (1.0)*T2(m,n,e,b)*V(e,j,m,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(n,b,j,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1*( (V{einf}*T2{amfn}) * T2+{mejb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(e,i,a,m) = INTERMVOVO(e,i,a,m) -
     &          (1.0)*V(e,n,i,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(e,i,a,m)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (T2+{menb}*T2{amfn}) * V{eijf})
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
        INTERMVVVV(e,b,a,f) = INTERMVVVV(e,b,a,f) +
     & (1.0/2.0)*T2(m,n,e,b)*T2(m,n,a,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVVV(e,b,a,f)*V(e,j,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     1*( (T2+{mejb}*V{emnf}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(j,b,n,f) = INTERMOVOV(j,b,n,f) +
     & T2(m,j,e,b)*V(e,n,m,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(j,b,n,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1*( (T2+{mejf}*V{emnb}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(j,f,n,b) = INTERMOVOV(j,f,n,b) -
     &          (1.0)*T2(m,j,e,f)*V(e,n,m,b)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(j,f,n,b)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1*( (V{amnf}*T2{eifn}) * T2+{mejb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,m,e,i) = INTERMVOVO(a,m,e,i) -
     &          (1.0)*V(a,n,m,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(a,m,e,i)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (T2+{mejf}*T2{eifn}) * V{amnb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,j,i,n) = INTERMOOOO(m,j,i,n) +
     & (1.0/2.0)*T2(m,j,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOOO(m,j,i,n)*V(a,n,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1*( (T2+{menb}*T2{aoen}) * V{oijm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,b,a,o) = INTERMOVVO(m,b,a,o) -
     &          (1.0)*T2(m,n,e,b)*T2(o,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVO(m,b,a,o)*V(o,j,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1*( (T2+{mejf}*T2{gifm}) * V{ageb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(e,j,g,i) = INTERMVOVO(e,j,g,i) -
     &          (1.0)*T2(m,j,e,f)*T2(i,m,g,f)
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
        DO  g  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(e,j,g,i)*V(a,e,g,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (V{niom}*T2{aneo}) * T2+{mejb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVV(i,m,a,e) = INTERMOOVV(i,m,a,e) +
     & (1.0/2.0)*V(n,o,i,m)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVV(i,m,a,e)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/2*( (T2+{menb}*V{omjn}) * T2{aieo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOO(e,b,o,j) = INTERMVVOO(e,b,o,j) +
     & (1.0/2.0)*T2(m,n,e,b)*V(o,j,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOO(e,b,o,j)*T2(i,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     1/2*( (T2+{mejf}*V{egfb}) * T2{aigm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,g,b) = INTERMOOVV(m,j,g,b) +
     & (1.0/2.0)*T2(m,j,e,f)*V(e,f,g,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVV(m,j,g,b)*T2(i,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/2*( (T2+{mejb}*T2{emfn}) * V{ainf})
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
        INTERMOVVO(j,b,f,n) = INTERMOVVO(j,b,f,n) +
     & (1.0/2.0)*T2(m,j,e,b)*T2(m,n,e,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVO(j,b,f,n)*V(a,n,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1/2*( (V{afeg}*T2{figm}) * T2+{mejb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,e,i,m) = INTERMVVOO(a,e,i,m) +
     & (1.0/2.0)*V(a,e,f,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOO(a,e,i,m)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     1*( (T2+{menb}*T2{eifn}) * V{amjf})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,b,i,f) = INTERMOVOV(m,b,i,f) +
     & T2(m,n,e,b)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(m,b,i,f)*V(a,j,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1*( (T2+{mejf}*T2{amfn}) * V{einb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(e,j,a,n) = INTERMVOVO(e,j,a,n) +
     & T2(m,j,e,f)*T2(m,n,a,f)
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
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(e,j,a,n)*V(e,n,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (T2+{menf}*V{emjb}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(n,f,j,b) = INTERMOVOV(n,f,j,b) +
     & (1.0/2.0)*T2(m,n,e,f)*V(e,j,m,b)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(n,f,j,b)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{eifn}) * V{amjb})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOO(m,i)*V(a,j,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{mejf}*T2{emfn}) * V{ainb})
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
        INTERMOO(j,n) = INTERMOO(j,n) -
     &          (1.0/4.0)*T2(m,j,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOO(j,n)*V(a,n,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{menb}*T2{emfn}) * V{aijf})
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
        INTERMVV(b,f) = INTERMVV(b,f) -
     &          (1.0/4.0)*T2(m,n,e,b)*T2(m,n,e,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVV(b,f)*V(a,j,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{amfn}) * V{eijb})
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
        INTERMVV(e,a) = INTERMVV(e,a) -
     &          (1.0/4.0)*T2(m,n,e,f)*T2(m,n,a,f)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVV(e,a)*V(e,j,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{mb}*T2{aiem}) * F{je})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) -
!     & (1.0/2.0)*T1(m,b)*T2(i,m,a,e)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOV(b,a,i,e)*FOCK(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
c     -1/2*( (F{mb}*T2{aiem}) * T1+{je})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
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
        INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) -
     &          (1.0/2.0)*FOCK(m,b)*T2(i,m,a,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(b,a,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
!
c     -1/2*( (T2+{mejb}*F{am}) * T1{ei})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVOVV(e,j,b,a) = INTERMVOVV(e,j,b,a) -
     &          (1.0/2.0)*T2(m,j,e,b)*FOCK(a,m)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVV(e,j,b,a)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T2+{mejb}*T1{am}) * F{ei})
! c      Commutator = 2, PT Order = 3
! !
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
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
!         INTERMVOVV(e,j,b,a) = INTERMVOVV(e,j,b,a) -
!     & (1.0/2.0)*T2(m,j,e,b)*T1(m,a)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOVV(e,j,b,a)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
c     1/4*( (T1+{je}*T2{amen}) * V{minb})
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
        INTERMOVOO(j,a,m,n) = INTERMOVOO(j,a,m,n) +
     & (1.0/4.0)*T1(j,e)*T2(m,n,a,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOO(j,a,m,n)*V(m,n,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1*( (T1+{me}*T2{aien}) * V{nmjb})
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOO(m,a,i,n)*V(n,j,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T1+{me}*T2{amen}) * V{nijb})
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
     & (1.0/2.0)*T1(m,e)*T2(m,n,a,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVO(a,n)*V(n,j,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1/4*( (T2+{menb}*T1{ei}) * V{amjn})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,b,i) = INTERMOOVO(m,n,b,i) +
     & (1.0/4.0)*T2(m,n,e,b)*T1(i,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVO(m,n,b,i)*V(a,j,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/2*( (T2+{menb}*T1{en}) * V{aijm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(m,b) = INTERMOV(m,b) +
     & (1.0/2.0)*T2(m,n,e,b)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOV(m,b)*V(a,j,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1*( (T1+{me}*V{efjb}) * T2{aifm})
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
        INTERMOVOV(m,f,j,b) = INTERMOVOV(m,f,j,b) +
     & T1(m,e)*V(e,j,f,b)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(m,f,j,b)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1/2*( (T1+{me}*T2{fiem}) * V{afjb})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,i) = INTERMVO(f,i) +
     & (1.0/2.0)*T1(m,e)*T2(i,m,f,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVO(f,i)*V(a,j,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1*( (T2+{mejb}*T1{en}) * V{ainm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,j,b,n) = INTERMOOVO(m,j,b,n) -
     &          (1.0)*T2(m,j,e,b)*T1(n,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVO(m,j,b,n)*V(a,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/2*( (T1+{je}*T2{fiem}) * V{afmb})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(j,f,i,m) = INTERMOVOO(j,f,i,m) -
     &          (1.0/2.0)*T1(j,e)*T2(i,m,f,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOO(j,f,i,m)*V(a,m,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( (T1+{mb}*T2{eifm}) * V{aejf})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOV(b,e,i,f) = INTERMVVOV(b,e,i,f) -
!     & (1.0/4.0)*T1(m,b)*T2(i,m,e,f)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOV(b,e,i,f)*V(a,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
c     -1*( (T1+{je}*V{efmb}) * T2{aifm})
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
        INTERMOVOV(j,f,m,b) = INTERMOVOV(j,f,m,b) -
     &          (1.0)*T1(j,e)*V(e,m,f,b)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(j,f,m,b)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T2+{mbne}*T1{an}) * V{eijm})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,b,e,a) = INTERMOVVV(m,b,e,a) +
!      & (1.0/2.0)*T2(m,n,b,e)*T1(n,a)
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
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVVV(m,b,e,a)*V(e,j,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{mb}*V{nmje}) * T2{aien})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(b,n,j,e) = INTERMVOOV(b,n,j,e) +
!      & T1(m,b)*V(n,j,m,e)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOOV(b,n,j,e)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( (T1+{mb}*T2{anem}) * V{nije})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,a,n,e) = INTERMVVOV(b,a,n,e) +
!      & (1.0/2.0)*T1(m,b)*T2(n,m,a,e)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOV(b,a,n,e)*V(n,j,i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
c     -1/2*( (T2+{mejf}*T1{fm}) * V{aieb})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) -
     &          (1.0/2.0)*T2(m,j,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVO(e,j)*V(a,e,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{einm}*T1{an}) * T2+{mejb})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(e,i,m,a) = INTERMVOOV(e,i,m,a) +
!      & V(e,n,i,m)*T1(n,a)
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
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOOV(e,i,m,a)*T2(m,j,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
c     1*( (V{aief}*T1{fm}) * T2+{mejb})
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
        INTERMVOVO(a,i,e,m) = INTERMVOVO(a,i,e,m) +
     & V(a,e,i,f)*T1(m,f)
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
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(a,i,e,m)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1*( (V{amef}*T1{fi}) * T2+{mejb})
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
        INTERMVOVO(a,m,e,i) = INTERMVOVO(a,m,e,i) -
     &          (1.0)*V(a,e,m,f)*T1(i,f)
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
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(a,m,e,i)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1/2*( (T2+{mejf}*T1{ei}) * V{amfb})
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
        INTERMOOVO(m,j,f,i) = INTERMOOVO(m,j,f,i) -
     &          (1.0/2.0)*T2(m,j,e,f)*T1(i,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVO(m,j,f,i)*V(a,f,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( (T2+{mejf}*T1{am}) * V{eifb})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,j,f,a) = INTERMVOVV(e,j,f,a) -
!     & (1.0/4.0)*T2(m,j,e,f)*T1(m,a)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOVV(e,j,f,a)*V(e,f,i,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (T1+{mb}*V{aiem}) * T1+{je})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) -
!     & (1.0)*T1(m,b)*V(a,e,i,m)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOV(b,a,i,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
c     1*( (T1+{je}*V{afeb}) * T1{fi})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(j,a,f,b) = INTERMOVVV(j,a,f,b) +
     & T1(j,e)*V(a,e,f,b)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVV(j,a,f,b)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (T1+{mb}*V{nijm}) * T1{an})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOO(NOS+1:NOS+NVS,1:NOS,1:NOS,1:NOS))
!         INTERMVOOO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOOO(b,n,i,j) = INTERMVOOO(b,n,i,j) +
!      & T1(m,b)*V(n,j,i,m)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOOO(b,n,i,j)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{mejb}*T1{am}) * T1{ei})
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
!         INTERMVOVV(e,j,b,a) = INTERMVOVV(e,j,b,a) -
!     & (1.0)*V(m,j,e,b)*T1(m,a)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOVV(e,j,b,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{me}*T1{am}) * V{eijb})
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
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/2.0)*T1(m,e)*T1(m,a)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVV(e,a)*V(e,j,i,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{mb}*T1{em}) * V{aije})
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
!         INTERMVV(b,e) = INTERMVV(b,e) -
!     & (1.0/2.0)*T1(m,b)*T1(m,e)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVV(b,e)*V(a,j,i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
c     -1/2*( (T1+{me}*T1{ei}) * V{amjb})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &          (1.0/2.0)*T1(m,e)*T1(i,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOO(m,i)*V(a,j,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (T1+{mb}*V{aejm}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(b,a,e,j) = INTERMVVVO(b,a,e,j) -
!     & (1.0)*T1(m,b)*V(a,j,e,m)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVVO(b,a,e,j)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{mieb}*T1{am}) * T1+{je})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(i,e,b,a) = INTERMOVVV(i,e,b,a) -
!     & (1.0)*V(m,e,i,b)*T1(m,a)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVVV(i,e,b,a)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
c     -1/2*( (T1+{je}*T1{em}) * V{aimb})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(j,m) = INTERMOO(j,m) -
     &          (1.0/2.0)*T1(j,e)*T1(m,e)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOO(j,m)*V(a,m,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/6*( ( (F{nf}*T2{amfn}) * T2+{mejb}) * T1{ei})
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
        INTERMVO(a,m) = INTERMVO(a,m) -
     &          (1.0/6.0)*FOCK(n,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
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
        INTERMVVOV(a,e,j,b) = INTERMVVOV(a,e,j,b) +
     &  INTERMVO(a,m)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(a,e,j,b)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T1+{me}*T2{aien}) * F{fm}) * T2+{nfjb})
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
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(a,i,n,f) = INTERMVOOV(a,i,n,f) +
     &  INTERMOVOO(m,a,i,n)*FOCK(f,m)
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
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOOV(a,i,n,f)*T2(n,j,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (F{nf}*T1{fm}) * T2+{mejb}) * T2{aien})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,m) = INTERMOO(n,m) -
     &          (2.0/3.0)*FOCK(n,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOV(n,e,j,b) = INTERMOVOV(n,e,j,b) +
     &  INTERMOO(n,m)*T2(m,j,e,b)
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
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(n,e,j,b)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T1+{me}*F{en}) * T2+{nfjb}) * T2{aifm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,n) = INTERMOO(m,n) -
     &          (2.0/3.0)*T1(m,e)*FOCK(e,n)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,f,j,b) = INTERMOVOV(m,f,j,b) +
     &  INTERMOO(m,n)*T2(n,j,f,b)
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
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOV(m,f,j,b)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{menb}*T1{em}) * T2{aifn}) * F{jf})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(n,b) = INTERMOV(n,b) -
     &          (1.0/6.0)*T2(m,n,e,b)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
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
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(b,a,i,f) = INTERMVVOV(b,a,i,f) +
     &  INTERMOV(n,b)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(b,a,i,f)*FOCK(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{menb}*T1{ei}) * T2{amfn}) * F{jf})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,b,i) = INTERMOOVO(m,n,b,i) +
     & (1.0/12.0)*T2(m,n,e,b)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(b,i,a,f) = INTERMVOVV(b,i,a,f) +
     &  INTERMOOVO(m,n,b,i)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVV(b,i,a,f)*FOCK(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*T2{anem}) * T2+{nfjb}) * F{fi})
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
     &          (1.0/6.0)*T1(m,e)*T2(n,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
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
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(a,f,j,b) = INTERMVVOV(a,f,j,b) +
     &  INTERMVO(a,n)*T2(n,j,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(a,f,j,b)*FOCK(f,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (F{nf}*T2{eifn}) * T2+{mejb}) * T1{am})
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
!         INTERMVO(e,i) = INTERMVO(e,i) -
!     & (1.0/6.0)*FOCK(n,f)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOV(i,m,j,b) = INTERMOOOV(i,m,j,b) +
!      &  INTERMVO(e,i)*T2(m,j,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOOOV(i,m,j,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
!         ENDIF
!
c     -2/3*( ( (T2+{mejb}*T1{en}) * F{nf}) * T2{aifm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,j,b,n) = INTERMOOVO(m,j,b,n) -
     &          (2.0/3.0)*T2(m,j,e,b)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,b,f) = INTERMOOVV(m,j,b,f) +
     &  INTERMOOVO(m,j,b,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVV(m,j,b,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{mejf}*T1{em}) * T2{aifn}) * F{nb})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) -
     &          (1.0/6.0)*T2(m,j,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(j,a,i,n) = INTERMOVOO(j,a,i,n) +
     &  INTERMOV(j,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVOO(j,a,i,n)*FOCK(n,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/6*( ( (T2+{nejf}*F{en}) * T2{aifm}) * T1+{mb})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(j,f) = INTERMOV(j,f) -
!     & (1.0/6.0)*T2(n,j,e,f)*FOCK(e,n)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(j,a,i,m) = INTERMOVOO(j,a,i,m) +
!      &  INTERMOV(j,f)*T2(i,m,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVOO(j,a,i,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
!         ENDIF
!
c     1/6*( ( (T1+{me}*T2{eifm}) * T2+{nfjb}) * F{an})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) +
     & (1.0/6.0)*T1(m,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(i,n,j,b) = INTERMOOOV(i,n,j,b) +
     &  INTERMOV(i,f)*T2(n,j,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOOV(i,n,j,b)*FOCK(a,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{mfnb}*F{fm}) * T2{aien}) * T1+{je})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,b) = INTERMOV(n,b) -
     &          (1.0/6.0)*T2(m,n,f,b)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO; ENDDO
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
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) +
     &  INTERMOV(n,b)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(b,a,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (T2+{mejf}*T1{am}) * T2{eifn}) * F{nb})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,j,f,a) = INTERMVOVV(e,j,f,a) +
!      & (1.0/12.0)*T2(m,j,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(j,a,i,n) = INTERMOVOO(j,a,i,n) +
!      &  INTERMVOVV(e,j,f,a)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVOO(j,a,i,n)*FOCK(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (T1+{mb}*T2{eifm}) * T2+{nejf}) * F{an})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOV(b,e,i,f) = INTERMVVOV(b,e,i,f) +
!      & (1.0/12.0)*T1(m,b)*T2(i,m,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOOO(NOS+1:NOS+NVS,1:NOS,1:NOS,1:NOS))
!         INTERMVOOO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOOO(b,i,n,j) = INTERMVOOO(b,i,n,j) +
!      &  INTERMVVOV(b,e,i,f)*T2(n,j,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOOO(b,i,n,j)*FOCK(a,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOO)
!         ENDIF
!         ENDIF
!
c     1/12*( ( (T1+{je}*T2{amen}) * T2+{mfnb}) * F{fi})
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
        INTERMOVOO(j,a,m,n) = INTERMOVOO(j,a,m,n) +
     & (1.0/12.0)*T1(j,e)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(j,a,f,b) = INTERMOVVV(j,a,f,b) +
     &  INTERMOVOO(j,a,m,n)*T2(m,n,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVV(j,a,f,b)*FOCK(f,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (F{jf}*T2{eifn}) * T2+{menb}) * T1{am})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(j,e,i,n) = INTERMOVOO(j,e,i,n) +
!      & (1.0/6.0)*FOCK(j,f)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOV(j,i,m,b) = INTERMOOOV(j,i,m,b) +
!      &  INTERMOVOO(j,e,i,n)*T2(m,n,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOOOV(j,i,m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
!         ENDIF
!
c     1/6*( ( (T2+{menb}*T2{aien}) * T1{fm}) * F{jf})
c      Commutator = 3, PT Order = 4

!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,b,a,i) = INTERMOVVO(m,b,a,i) +
     & (1.0/6.0)*T2(m,n,e,b)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(b,a,i,f) = INTERMVVOV(b,a,i,f) +
     &  INTERMOVVO(m,b,a,i)*T1(m,f)
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
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(b,a,i,f)*FOCK(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{nf}*T1{an}) * T2{eifm}) * T2+{mejb})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(f,a) = INTERMVV(f,a) +
!      & (1.0/2.0)*FOCK(n,f)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOO(a,e,i,m) = INTERMVVOO(a,e,i,m) +
!      &  INTERMVV(f,a)*T2(i,m,e,f)
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
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOO(a,e,i,m)*T2(m,j,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
!
c     1/2*( ( (T2+{mejb}*T2{amen}) * F{nf}) * T1{fi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(j,b,a,n) = INTERMOVVO(j,b,a,n) +
     & (1.0/2.0)*T2(m,j,e,b)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(j,b,a,f) = INTERMOVVV(j,b,a,f) +
     &  INTERMOVVO(j,b,a,n)*FOCK(n,f)
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
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVV(j,b,a,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{nb}*T1{en}) * T2+{mejf}) * T2{aifm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(b,e) = INTERMVV(b,e) +
     & (1.0/6.0)*FOCK(n,b)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(b,m,j,f) = INTERMVOOV(b,m,j,f) +
     &  INTERMVV(b,e)*T2(m,j,e,f)
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
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOOV(b,m,j,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
        ENDIF
!
c     1/2*( ( (T2+{mfnb}*T2{aifn}) * F{em}) * T1+{je})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(m,b,a,i) = INTERMOVVO(m,b,a,i) +
     & (1.0/2.0)*T2(m,n,f,b)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
     &  INTERMOVVO(m,b,a,i)*FOCK(e,m)
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
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(b,a,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{nejf}*F{ei}) * T2{anfm}) * T1+{mb})
! c      Commutator = 3, PT Order = 4
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO(n,j,f,i) = INTERMOOVO(n,j,f,i) +
!      & (1.0/6.0)*T2(n,j,e,f)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOVO2(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMOOVO2 = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVO2(j,i,a,m) = INTERMOOVO2(j,i,a,m) +
!      &  INTERMOOVO(n,j,f,i)*T2(n,m,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOOVO2(j,i,a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVO2)
!         ENDIF
!         ENDIF
!
c     1/6*( ( (T2+{mejf}*T1{ei}) * T2{amfn}) * F{nb})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,j,f,i) = INTERMOOVO(m,j,f,i) +
     & (1.0/6.0)*T2(m,j,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO2(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO2(j,i,a,n) = INTERMOOVO2(j,i,a,n) +
     &  INTERMOOVO(m,j,f,i)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOVO2(j,i,a,n)*FOCK(n,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO2)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{mb}*F{em}) * T2+{nejf}) * T2{aifn})
! c      Commutator = 3, PT Order = 4
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
!         INTERMVV(b,e) = INTERMVV(b,e) +
!      & (1.0/2.0)*T1(m,b)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOOV(b,n,j,f) = INTERMVOOV(b,n,j,f) +
!      &  INTERMVV(b,e)*T2(n,j,e,f)
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
!         DO  f  = NOS+1, NOS+NVS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOOV(b,n,j,f)*T2(i,n,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!         ENDIF
!
c     -1/6*( ( (T2+{nfjb}*T2{amfn}) * T1+{me}) * F{ei})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVO(j,b,a,m) = INTERMOVVO(j,b,a,m) -
     &          (1.0/6.0)*T2(n,j,f,b)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVV(j,b,a,e) = INTERMOVVV(j,b,a,e) +
     &  INTERMOVVO(j,b,a,m)*T1(m,e)
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
        DO  e  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOVVV(j,b,a,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*F{am}) * T2{eifn}) * T2+{nfjb})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) -
     &          (1.0/6.0)*T1(m,e)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,i,f,n) = INTERMVOVO(a,i,f,n) +
     &  INTERMVV(e,a)*T2(i,n,e,f)
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
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVOVO(a,i,f,n)*T2(n,j,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{je}*T2{fien}) * T2+{mfnb}) * F{am})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(j,f,i,n) = INTERMOVOO(j,f,i,n) +
     & (1.0/6.0)*T1(j,e)*T2(i,n,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(j,i,m,b) = INTERMOOOV(j,i,m,b) +
     &  INTERMOVOO(j,f,i,n)*T2(m,n,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMOOOV(j,i,m,b)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{mbne}*T1{an}) * T1{ei}) * F{jm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,b,e,a) = INTERMOVVV(m,b,e,a) +
!      & (1.0/3.0)*T2(m,n,b,e)*T1(n,a)
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
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVO(m,b,a,i) = INTERMOVVO(m,b,a,i) +
!      &  INTERMOVVV(m,b,e,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVVO(m,b,a,i)*FOCK(j,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
!
c     -1/6*( ( (T1+{me}*F{eb}) * T2{aifm}) * T1+{jf})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(m,b) = INTERMOV(m,b) -
     &          (1.0/6.0)*T1(m,e)*FOCK(e,b)
        ENDDO; ENDDO; ENDDO
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
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(b,a,i,f) = INTERMVVOV(b,a,i,f) +
     &  INTERMOV(m,b)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(b,a,i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{nb}*T2{aien}) * T1+{me}) * F{jm})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) +
!      & (1.0/6.0)*T1(n,b)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!      &  INTERMVVOV(b,a,i,e)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOO(b,a,i,m)*FOCK(j,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T1+{mb}*T2{anem}) * T1+{je}) * F{ni})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,a,n,e) = INTERMVVOV(b,a,n,e) +
!      & (1.0/3.0)*T1(m,b)*T2(n,m,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         INTERMVVOO(b,a,n,j) = INTERMVVOO(b,a,n,j) +
!      &  INTERMVVOV(b,a,n,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOO(b,a,n,j)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T1+{je}*F{ef}) * T2{aifm}) * T1+{mb})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(j,f) = INTERMOV(j,f) -
!     & (1.0/2.0)*T1(j,e)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(j,a,i,m) = INTERMOVOO(j,a,i,m) +
!      &  INTERMOV(j,f)*T2(i,m,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVOO(j,a,i,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (T1+{mb}*F{nm}) * T2{aien}) * T1+{je})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVO(b,n) = INTERMVO(b,n) +
!      & (1.0/2.0)*T1(m,b)*FOCK(n,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(b,a,i,e) = INTERMVVOV(b,a,i,e) +
!      &  INTERMVO(b,n)*T2(i,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOV(b,a,i,e)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (F{af}*T2{fiem}) * T1+{mb}) * T1+{je})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(a,i,e,m) = INTERMVOVO(a,i,e,m) -
!     & (1.0/3.0)*FOCK(a,f)*T2(i,m,f,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         INTERMVOVV(a,i,e,b) = INTERMVOVV(a,i,e,b) +
!      &  INTERMVOVO(a,i,e,m)*T1(m,b)
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
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOVV(a,i,e,b)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{mejb}*T1{am}) * T1{en}) * F{ni})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
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
!         INTERMVOVV(e,j,b,a) = INTERMVOVV(e,j,b,a) +
!      & (1.0/6.0)*T2(m,j,e,b)*T1(m,a)
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
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVO(j,b,a,n) = INTERMOVVO(j,b,a,n) +
!      &  INTERMVOVV(e,j,b,a)*T1(n,e)
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOVVO(j,b,a,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{mejf}*F{fb}) * T1{am}) * T1{ei})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOV(m,e,j,b) = INTERMOVOV(m,e,j,b) -
!     & (1.0/3.0)*T2(m,j,e,f)*FOCK(f,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         INTERMVOVV(e,j,b,a) = INTERMVOVV(e,j,b,a) +
!      &  INTERMOVOV(m,e,j,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVOVV(e,j,b,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/2*( ( (F{nm}*T1{an}) * T2+{mejb}) * T1{ei})
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
!         INTERMOV(m,a) = INTERMOV(m,a) +
!      & (1.0/2.0)*FOCK(n,m)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(a,e,j,b) = INTERMVVOV(a,e,j,b) +
!      &  INTERMOV(m,a)*T2(m,j,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMVVOV(a,e,j,b)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (F{ef}*T1{fi}) * T2+{mejb}) * T1{am})
! c      Commutator = 3, PT Order = 5
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,i) = INTERMVO(e,i) -
!     & (1.0/2.0)*FOCK(e,f)*T1(i,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOOOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOOV(i,m,j,b) = INTERMOOOV(i,m,j,b) +
!      &  INTERMVO(e,i)*T2(m,j,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
!         VT(a,j,i,b) = VT(a,j,i,b) +
!      &  INTERMOOOV(i,m,j,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOOV)
!         ENDIF
!
c     -1/6*( ( (F{af}*T1{fm}) * T2+{mejb}) * T1{ei})
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
        INTERMVO(a,m) = INTERMVO(a,m) -
     &          (1.0/6.0)*FOCK(a,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
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
        INTERMVVOV(a,e,j,b) = INTERMVVOV(a,e,j,b) +
     &  INTERMVO(a,m)*T2(m,j,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
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
        VT(a,j,i,b) = VT(a,j,i,b) +
     &  INTERMVVOV(a,e,j,b)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF

        DO  a  = NOS+1, NOS+NVS
        DO  b  = NOS+1, NOS+NVS
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        VT(a,j,i,b) = VT(a,j,i,b)
        VT(i,b,a,j) = VT(a,j,i,b)
        VT(a,j,b,i) = -VT(a,j,i,b)
        VT(j,a,i,b) = -VT(a,j,i,b)
        ENDDO; ENDDO; ENDDO; ENDDO

      END
