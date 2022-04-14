      SUBROUTINE HHHP(NOS,NVS,NVAS,NACTV,T1,T2,
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
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO2(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO2(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVO(:,:,:,:)

!
c     1*V{kjai}
c      Commutator = 0, PT Order = 1
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) + V(k,a,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*F{ke}*T2{ejai}
c      Commutator = 1, PT Order = 1
!
        IF(NONHF)THEN
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     & FOCK(k,e)*T2(j,i,e,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
        ENDIF
!
c     -1*V{mjke}*T2{aiem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) -
     &          (1.0)*V(m,k,j,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{mike}*T2{ajem}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     & V(m,k,i,e)*T2(j,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*V{aekf}*T2{fjei}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     & (1.0/2.0)*V(a,k,e,f)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*V{kjmi}*T1{am}
! c      Commutator = 1, PT Order = 3
! !
!         IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         VT(k,a,j,i) = VT(k,a,j,i) -
!     & (1.0)*V(k,m,j,i)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     1*T1+{ke}*V{ejai}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     & T1(k,e)*V(e,a,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*V{ajke}*T1{ei}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) -
     &          (1.0)*V(a,k,j,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{aike}*T1{ej}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     & V(a,k,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*( (T2+{mekf}*F{em}) * T2{fjai})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     & (1.0/2.0)*T2(m,k,e,f)*FOCK(e,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,f)*T2(j,i,f,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/4*( (T2+{mekf}*F{am}) * T2{fjei})
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
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
     &          (1.0/4.0)*T2(m,k,e,f)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOVV(e,k,f,a)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     -1/2*( (T2+{mekf}*F{ej}) * T2{aifm})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,j) = INTERMOOVO(m,k,f,j) -
     &          (1.0/2.0)*T2(m,k,e,f)*FOCK(e,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(m,k,f,j)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/2*( (T2+{mekf}*F{ei}) * T2{ajfm})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,i) = INTERMOOVO(m,k,f,i) +
     & (1.0/2.0)*T2(m,k,e,f)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(m,k,f,i)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1/2*( (T1+{ke}*T2{aiem}) * F{mj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(k,a,i,m) = INTERMOVOO(k,a,i,m) -
     &          (1.0/2.0)*T1(k,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,a,i,m)*FOCK(m,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T1+{ke}*T2{ajem}) * F{mi})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(k,a,j,m) = INTERMOVOO(k,a,j,m) +
     & (1.0/2.0)*T1(k,e)*T2(j,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,a,j,m)*FOCK(m,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (F{af}*T2{ejfi}) * T1+{ke})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,e,j,i) = INTERMVVOO(a,e,j,i) +
     & (1.0/2.0)*FOCK(a,f)*T2(j,i,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVOO(a,e,j,i)*T1(k,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/2*( (T1+{me}*F{km}) * T2{aiej})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(e,k) = INTERMVO(e,k) -
     &          (1.0/2.0)*T1(m,e)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVO(e,k)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1*( (T1+{ke}*F{ef}) * T2{ajfi})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) -
     &          (1.0)*T1(k,e)*FOCK(e,f)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,f)*T2(j,i,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T2+{mekf}*V{emfg}) * T2{aigj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(k,g) = INTERMOV(k,g) +
     & (1.0/2.0)*T2(m,k,e,f)*V(e,f,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,g)*T2(i,j,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1*( (T2+{mekf}*V{ameg}) * T2{fjgi})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVVV(k,f,a,g) = INTERMOVVV(k,f,a,g) -
     &          (1.0)*T2(m,k,e,f)*V(a,e,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVVV(k,f,a,g)*T2(j,i,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     1/4*( (T2+{menf}*V{fmkn}) * T2{aiej})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,k) = INTERMVO(e,k) +
     & (1.0/4.0)*T2(m,n,e,f)*V(f,k,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVO(e,k)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1/8*( (T2+{menf}*V{amkn}) * T2{fjei})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,1:NOS))
        INTERMVVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVO(e,f,a,k) = INTERMVVVO(e,f,a,k) +
     & (1.0/8.0)*T2(m,n,e,f)*V(a,k,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVVO(e,f,a,k)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
!
c     1*( (T2+{mekf}*T2{amfn}) * V{njei})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(e,k,a,n) = INTERMVOVO(e,k,a,n) +
     & T2(m,k,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOVO(e,k,a,n)*V(n,e,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1/4*( (T2+{mekf}*T2{emfn}) * V{ainj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(k,n) = INTERMOO(k,n) -
     &          (1.0/4.0)*T2(m,k,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOO(k,n)*V(a,n,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{amfn}) * V{kjei})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVV(e,a)*V(k,e,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1*( (T2+{mekf}*V{ejnm}) * T2{aifn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(k,f,j,n) = INTERMOVOO(k,f,j,n) +
     & T2(m,k,e,f)*V(e,n,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,f,j,n)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1*( (T2+{mekf}*V{einm}) * T2{ajfn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(k,f,i,n) = INTERMOVOO(k,f,i,n) -
     &          (1.0)*T2(m,k,e,f)*V(e,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,f,i,n)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( (T2+{menf}*T2{aifn}) * V{ejkm})
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
        INTERMOVVO(m,e,a,i) = INTERMOVVO(m,e,a,i) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVVO(m,e,a,i)*V(e,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1/2*( (T2+{menf}*T2{ajfn}) * V{eikm})
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
        INTERMOVVO(m,e,a,j) = INTERMOVVO(m,e,a,j) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVVO(m,e,a,j)*V(e,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1/2*( (T2+{mekf}*V{ejfg}) * T2{aigm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,k,j,g) = INTERMOOOV(m,k,j,g) -
     &          (1.0/2.0)*T2(m,k,e,f)*V(e,f,j,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,k,j,g)*T2(i,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     1/2*( (T2+{mekf}*V{eifg}) * T2{ajgm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,k,i,g) = INTERMOOOV(m,k,i,g) +
     & (1.0/2.0)*T2(m,k,e,f)*V(e,f,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,k,i,g)*T2(j,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{ejfn}) * V{aikm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOO(m,j)*V(a,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/4*( (T2+{menf}*T2{eifn}) * V{ajkm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOO(m,i)*V(a,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1*( (T2+{mekf}*T2{fjgm}) * V{aieg})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOOV(e,k,j,g) = INTERMVOOV(e,k,j,g) +
     & T2(m,k,e,f)*T2(j,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOOV(e,k,j,g)*V(a,e,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     -1*( (T2+{mekf}*T2{figm}) * V{ajeg})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOOV(e,k,i,g) = INTERMVOOV(e,k,i,g) -
     &          (1.0)*T2(m,k,e,f)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOOV(e,k,i,g)*V(a,e,j,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     1/2*( (T2+{mekf}*T2{ejfn}) * V{ainm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,j,n) = INTERMOOOO(m,k,j,n) +
     & (1.0/2.0)*T2(m,k,e,f)*T2(j,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOO(m,k,j,n)*V(a,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/2*( (T2+{mekf}*T2{eifn}) * V{ajnm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,i,n) = INTERMOOOO(m,k,i,n) -
     &          (1.0/2.0)*T2(m,k,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOO(m,k,i,n)*V(a,n,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1*( (T1+{me}*T2{fjei}) * V{amkf})
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
        INTERMOVOO(m,f,j,i) = INTERMOVOO(m,f,j,i) -
     &          (1.0)*T1(m,e)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(m,f,j,i)*V(a,k,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/4*( (T2+{mekf}*T1{am}) * V{fjei})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
!     & (1.0/4.0)*T2(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVOVV(e,k,f,a)*V(f,e,j,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (V{mekf}*T1{am}) * T2{fjei})
! c      Commutator = 2, PT Order = 4
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
!     & (1.0/2.0)*V(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVOVV(e,k,f,a)*T2(j,i,f,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!
c     -1/2*( (T2+{mekf}*T1{fm}) * V{aiej})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,k) = INTERMVO(e,k) -
     &          (1.0/2.0)*T2(m,k,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVO(e,k)*V(a,e,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1/2*( (T1+{ke}*V{afeg}) * T2{gjfi})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVVV(k,a,f,g) = INTERMOVVV(k,a,f,g) +
     & (1.0/2.0)*T1(k,e)*V(a,e,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVVV(k,a,f,g)*T2(j,i,g,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     1*( (T1+{me}*V{emkf}) * T2{aifj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN

        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     & T1(m,e)*V(e,k,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,f)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1*( (V{mekf}*T1{em}) * T2{aifj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     & V(m,k,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,f)*T2(i,j,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/4*( (T1+{ke}*T2{anem}) * V{mjni})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(k,a,n,m) = INTERMOVOO(k,a,n,m) +
     & (1.0/4.0)*T1(k,e)*T2(n,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,a,n,m)*V(m,n,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T1+{me}*T2{amen}) * V{kjni})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVO(a,n)*V(k,n,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/2*( (T2+{mekf}*T1{fi}) * V{ajem})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,k,i) = INTERMOVOO(m,e,k,i) -
     &          (1.0/2.0)*T2(m,k,e,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(m,e,k,i)*V(a,e,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T2+{mekf}*T1{fj}) * V{aiem})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,k,j) = INTERMOVOO(m,e,k,j) +
     & (1.0/2.0)*T2(m,k,e,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(m,e,k,j)*V(a,e,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T1+{me}*T2{eifm}) * V{ajkf})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(i,f)*V(a,k,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( (T1+{me}*T2{ejfm}) * V{aikf})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(j,f)*V(a,k,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1*( (V{mekf}*T1{ei}) * T2{ajfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,i) = INTERMOOVO(m,k,f,i) +
     & V(m,k,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(m,k,f,i)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1*( (V{mekf}*T1{ej}) * T2{aifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,j) = INTERMOOVO(m,k,f,j) -
     &          (1.0)*V(m,k,e,f)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(m,k,f,j)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1*( (T1+{me}*V{eikf}) * T2{ajfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(m,i,k,f) = INTERMOOOV(m,i,k,f) +
     & T1(m,e)*V(e,k,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,i,k,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1*( (T1+{me}*V{ejkf}) * T2{aifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(m,j,k,f) = INTERMOOOV(m,j,k,f) -
     &          (1.0)*T1(m,e)*V(e,k,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,j,k,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1*( (T1+{me}*T2{ajen}) * V{nikm})
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
        INTERMOVOO(m,a,j,n) = INTERMOVOO(m,a,j,n) -
     &          (1.0)*T1(m,e)*T2(j,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(m,a,j,n)*V(n,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1*( (T1+{me}*T2{aien}) * V{njkm})
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
        INTERMOVOO(m,a,i,n) = INTERMOVOO(m,a,i,n) +
     & T1(m,e)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(m,a,i,n)*V(n,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T1+{ke}*T2{ejfm}) * V{aimf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(k,j,f,m) = INTERMOOVO(k,j,f,m) +
     & (1.0/2.0)*T1(k,e)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(k,j,f,m)*V(a,m,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/2*( (T1+{ke}*T2{eifm}) * V{ajmf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(k,i,f,m) = INTERMOOVO(k,i,f,m) -
     &          (1.0/2.0)*T1(k,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(k,i,f,m)*V(a,m,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1*( (T1+{ke}*V{ejmf}) * T2{aifm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(k,j,m,f) = INTERMOOOV(k,j,m,f) +
     & T1(k,e)*V(e,m,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(k,j,m,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1*( (T1+{ke}*V{eimf}) * T2{ajfm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(k,i,m,f) = INTERMOOOV(k,i,m,f) -
     &          (1.0)*T1(k,e)*V(e,m,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(k,i,m,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{ejmi}*T1{am}) * T1+{ke})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(e,j,i,a) = INTERMVOOV(e,j,i,a) -
!     & (1.0)*V(e,m,j,i)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVOOV(e,j,i,a)*T1(k,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
c     1*( (V{aekf}*T1{fj}) * T1{ei})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,e,k,j) = INTERMVVOO(a,e,k,j) +
     & V(a,k,e,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVOO(a,e,k,j)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1*( (T1+{ke}*V{ajef}) * T1{fi})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(k,a,j,f) = INTERMOVOV(k,a,j,f) -
     &          (1.0)*T1(k,e)*V(a,e,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(k,a,j,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1*( (T1+{ke}*V{aief}) * T1{fj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(k,a,i,f) = INTERMOVOV(k,a,i,f) +
     & T1(k,e)*V(a,e,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(k,a,i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1*( (V{mike}*T1{am}) * T1{ej})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(i,k,e,a) = INTERMOOVV(i,k,e,a) -
!     & (1.0)*V(m,k,i,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMOOVV(i,k,e,a)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1*( (V{mjke}*T1{am}) * T1{ei})
! c      Commutator = 2, PT Order = 5
! !
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(j,k,e,a) = INTERMOOVV(j,k,e,a) +
!      & V(m,k,j,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMOOVV(j,k,e,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( (T1+{me}*T1{am}) * V{kjei})
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
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVV(e,a)*V(k,e,j,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
c     -1/2*( (T1+{ke}*T1{em}) * V{aimj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(k,m) = INTERMOO(k,m) -
     &          (1.0/2.0)*T1(k,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOO(k,m)*V(a,m,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/2*( (T1+{me}*T1{ej}) * V{aikm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOO(m,j)*V(a,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/2*( (T1+{me}*T1{ei}) * V{ajkm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOO(m,i)*V(a,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/24*( ( (F{kg}*T2{amgn}) * T2+{menf}) * T2{fjei})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(k,a,m,n) = INTERMOVOO(k,a,m,n) +
     & (1.0/24.0)*FOCK(k,g)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(k,a,e,f) = INTERMOVVV(k,a,e,f) +
     &  INTERMOVOO(k,a,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVVV(k,a,e,f)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     1/12*( ( (F{kg}*T2{fmgn}) * T2+{menf}) * T2{aiej})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(k,f,m,n) = INTERMOVOO(k,f,m,n) +
     & (1.0/12.0)*FOCK(k,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,e) = INTERMOV(k,e) +
     &  INTERMOVOO(k,f,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,e)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (F{ng}*T2{amgn}) * T2+{mekf}) * T2{fjei})
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
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(a,m) = INTERMVO(a,m) -
     &          (1.0/12.0)*FOCK(n,g)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(a,e,k,f) = INTERMVVOV(a,e,k,f) +
     &  INTERMVO(a,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVOV(a,e,k,f)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{amfn}) * T2{gjei}) * F{kg})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,g,j,i) = INTERMVVOO(a,g,j,i) +
     &  INTERMVV(e,a)*T2(j,i,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVOO(a,g,j,i)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{menf}*T2{ejfm}) * T2{aign}) * F{kg})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(j,a,i,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{eifm}) * T2{ajgn}) * F{kg})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(i,a,j,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     1/3*( ( (T2+{mekf}*T2{ejfn}) * F{ng}) * T2{aigm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,j,n) = INTERMOOOO(m,k,j,n) +
     & (1.0/3.0)*T2(m,k,e,f)*T2(j,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,k,j,g) = INTERMOOOV(m,k,j,g) +
     &  INTERMOOOO(m,k,j,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,k,j,g)*T2(i,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (T2+{mekf}*T2{eifn}) * F{ng}) * T2{ajgm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,i,n) = INTERMOOOO(m,k,i,n) -
     &          (1.0/3.0)*T2(m,k,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,k,i,g) = INTERMOOOV(m,k,i,g) +
     &  INTERMOOOO(m,k,i,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,k,i,g)*T2(j,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/4*( ( (T2+{mekf}*T2{emfn}) * F{ng}) * T2{aigj})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(k,n) = INTERMOO(k,n) -
     &          (1.0/4.0)*T2(m,k,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(k,g) = INTERMOV(k,g) +
     &  INTERMOO(k,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,g)*T2(i,j,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{ng}*T2{fmgn}) * T2+{mekf}) * T2{aiej})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(f,m) = INTERMVO(f,m) -
     &          (1.0/6.0)*FOCK(n,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO2(NOS+1:NOS+NVS,1:NOS))
        INTERMVO2 = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO2(e,k) = INTERMVO2(e,k) +
     &  INTERMVO(f,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVO2(e,k)*T2(i,j,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO2)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{ng}*T2{fjgn}) * T2+{mekf}) * T2{aiem})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(f,j) = INTERMVO(f,j) +
     & (1.0/6.0)*FOCK(n,g)*T2(j,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(j,m,e,k) = INTERMOOVO(j,m,e,k) +
     &  INTERMVO(f,j)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(j,m,e,k)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{ng}*T2{fign}) * T2+{mekf}) * T2{ajem})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(f,i) = INTERMVO(f,i) -
     &          (1.0/6.0)*FOCK(n,g)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(i,m,e,k) = INTERMOOVO(i,m,e,k) +
     &  INTERMVO(f,i)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(i,m,e,k)*T2(j,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{menf}*T2{aifn}) * T2{ejgm}) * F{kg})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        INTERMOVVO(m,e,a,i) = INTERMOVVO(m,e,a,i) -
     &          (1.0/6.0)*T2(m,n,e,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOOV(a,i,j,g) = INTERMVOOV(a,i,j,g) +
     &  INTERMOVVO(m,e,a,i)*T2(j,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOOV(a,i,j,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T2+{menf}*T2{ajfn}) * T2{eigm}) * F{kg})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        INTERMOVVO(m,e,a,j) = INTERMOVVO(m,e,a,j) +
     & (1.0/6.0)*T2(m,n,e,f)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOOV(a,j,i,g) = INTERMVOOV(a,j,i,g) +
     &  INTERMOVVO(m,e,a,j)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOOV(a,j,i,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
        ENDIF
!
c     2/3*( ( (T2+{mekf}*T2{amfn}) * F{ng}) * T2{gjei})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(e,k,a,n) = INTERMVOVO(e,k,a,n) +
     & (2.0/3.0)*T2(m,k,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVV(e,k,a,g) = INTERMVOVV(e,k,a,g) +
     &  INTERMVOVO(e,k,a,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOVV(e,k,a,g)*T2(j,i,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     2/3*( ( (F{ng}*T2{ejgm}) * T2+{mekf}) * T2{aifn})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(n,e,j,m) = INTERMOVOO(n,e,j,m) +
     & (2.0/3.0)*FOCK(n,g)*T2(j,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(n,j,k,f) = INTERMOOOV(n,j,k,f) +
     &  INTERMOVOO(n,e,j,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(n,j,k,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (F{ng}*T2{eigm}) * T2+{mekf}) * T2{ajfn})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(n,e,i,m) = INTERMOVOO(n,e,i,m) -
     &          (2.0/3.0)*FOCK(n,g)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(n,i,k,f) = INTERMOOOV(n,i,k,f) +
     &  INTERMOVOO(n,e,i,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(n,i,k,f)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/3*( ( (T2+{menf}*T1{en}) * F{km}) * T2{fjai})
c      Commutator = 3, PT Order = 4
!
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
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,k) = INTERMVO(f,k) +
     &  INTERMOV(m,f)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVO(f,k)*T2(j,i,f,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/3*( ( (T2+{menf}*T1{ej}) * T2{aifn}) * F{km})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,f,j) = INTERMOOVO(m,n,f,j) -
     &          (1.0/3.0)*T2(m,n,e,f)*T1(j,e)
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
        INTERMOOVO2(m,j,a,i) = INTERMOOVO2(m,j,a,i) +
     &  INTERMOOVO(m,n,f,j)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO2(m,j,a,i)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO2)
        ENDIF
!
c     1/3*( ( (T2+{menf}*T1{ei}) * T2{ajfn}) * F{km})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,f,i) = INTERMOOVO(m,n,f,i) +
     & (1.0/3.0)*T2(m,n,e,f)*T1(i,e)
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
        INTERMOOVO2(m,i,a,j) = INTERMOOVO2(m,i,a,j) +
     &  INTERMOOVO(m,n,f,i)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO2(m,i,a,j)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO2)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/12*( ( (T2+{mekf}*T1{am}) * T2{eifn}) * F{nj})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) +
!      & (1.0/12.0)*T2(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  k  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(k,a,i,n) = INTERMOVOO(k,a,i,n) +
!      &  INTERMVOVV(e,k,f,a)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  n  = 1, NOS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMOVOO(k,a,i,n)*FOCK(n,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/12*( ( (T2+{mekf}*T1{am}) * T2{ejfn}) * F{ni})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
!     & (1.0/12.0)*T2(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(k,a,j,n) = INTERMOVOO(k,a,j,n) +
!      &  INTERMVOVV(e,k,f,a)*T2(j,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  n  = 1, NOS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMOVOO(k,a,j,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
!
c     -1/12*( ( (F{ag}*T1{gm}) * T2+{mekf}) * T2{fjei})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(a,m) = INTERMVO(a,m) -
     &          (1.0/12.0)*FOCK(a,g)*T1(m,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(a,e,k,f) = INTERMVVOV(a,e,k,f) +
     &  INTERMVO(a,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVOV(a,e,k,f)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (T2+{mekf}*T1{am}) * F{eg}) * T2{fjgi})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
!     & (1.0/2.0)*T2(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMOVVV(k,f,a,g) = INTERMOVVV(k,f,a,g) +
!      &  INTERMVOVV(e,k,f,a)*FOCK(e,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMOVVV(k,f,a,g)*T2(j,i,f,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
c     -1/2*( ( (F{nm}*T1{en}) * T2+{mekf}) * T2{fjai})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) -
     &          (1.0/2.0)*FOCK(n,m)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV2(1:NOS,NOS+1:NOS+NVS))
        INTERMOV2 = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV2(k,f) = INTERMOV2(k,f) +
     &  INTERMOV(m,e)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV2(k,f)*T2(j,i,f,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV2)
        ENDIF
!
c     -1/2*( ( (T2+{mekf}*T1{fm}) * F{eg}) * T2{gjai})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,k) = INTERMVO(e,k) -
     &          (1.0/2.0)*T2(m,k,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(k,g) = INTERMOV(k,g) +
     &  INTERMVO(e,k)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,g)*T2(j,i,g,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( ( (F{eg}*T1{gm}) * T2+{mekf}) * T2{ajfi})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,m) = INTERMVO(e,m) -
     &          (1.0/2.0)*FOCK(e,g)*T1(m,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     &  INTERMVO(e,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,f)*T2(j,i,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( ( (T2+{mekf}*T1{fj}) * F{eg}) * T2{aigm})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,k,j) = INTERMOVOO(m,e,k,j) +
     & (1.0/2.0)*T2(m,k,e,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,k,j,g) = INTERMOOOV(m,k,j,g) +
     &  INTERMOVOO(m,e,k,j)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,k,j,g)*T2(i,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1/2*( ( (T2+{mekf}*T1{fi}) * F{eg}) * T2{ajgm})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,k,i) = INTERMOVOO(m,e,k,i) -
     &          (1.0/2.0)*T2(m,k,e,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,k,i,g) = INTERMOOOV(m,k,i,g) +
     &  INTERMOVOO(m,e,k,i)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(m,k,i,g)*T2(j,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     1/2*( ( (T2+{mekf}*T1{ej}) * F{nm}) * T2{aifn})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,j) = INTERMOOVO(m,k,f,j) +
     & (1.0/2.0)*T2(m,k,e,f)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(k,f,j,n) = INTERMOVOO(k,f,j,n) +
     &  INTERMOOVO(m,k,f,j)*FOCK(n,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,f,j,n)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( ( (T2+{mekf}*T1{ei}) * F{nm}) * T2{ajfn})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,i) = INTERMOOVO(m,k,f,i) -
     &          (1.0/2.0)*T2(m,k,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(k,f,i,n) = INTERMOVOO(k,f,i,n) +
     &  INTERMOOVO(m,k,f,i)*FOCK(n,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOO(k,f,i,n)*T2(j,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( ( (F{eg}*T1{gj}) * T2+{mekf}) * T2{aifm})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) -
     &          (1.0/2.0)*FOCK(e,g)*T1(j,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(j,m,k,f) = INTERMOOOV(j,m,k,f) +
     &  INTERMVO(e,j)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(j,m,k,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     1/2*( ( (F{eg}*T1{gi}) * T2+{mekf}) * T2{ajfm})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) +
     & (1.0/2.0)*FOCK(e,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(i,m,k,f) = INTERMOOOV(i,m,k,f) +
     &  INTERMVO(e,i)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(i,m,k,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T2+{menf}*T1{an}) * T2{ejfi}) * F{km})
! c      Commutator = 3, PT Order = 4
! !
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,f,a) = INTERMOVVV(m,e,f,a) +
!      & (1.0/6.0)*T2(m,n,e,f)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMOVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOVOO(m,a,j,i) = INTERMOVOO(m,a,j,i) +
!      &  INTERMOVVV(m,e,f,a)*T2(j,i,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMOVOO(m,a,j,i)*FOCK(k,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOO)
!         ENDIF
!
c     1/6*( ( (T2+{mekf}*T1{em}) * T2{fjgi}) * F{ag})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     & (1.0/6.0)*T2(m,k,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(k,j,g,i) = INTERMOOVO(k,j,g,i) +
     &  INTERMOV(k,f)*T2(j,i,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(k,j,g,i)*FOCK(a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/6*( ( (T2+{mekf}*T1{ei}) * T2{amfn}) * F{nj})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,i) = INTERMOOVO(m,k,f,i) +
     & (1.0/6.0)*T2(m,k,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO2(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO2(k,i,a,n) = INTERMOOVO2(k,i,a,n) +
     &  INTERMOOVO(m,k,f,i)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO2(k,i,a,n)*FOCK(n,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO2)
        ENDIF
!
c     -1/6*( ( (T2+{mekf}*T1{ej}) * T2{amfn}) * F{ni})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,k,f,j) = INTERMOOVO(m,k,f,j) -
     &          (1.0/6.0)*T2(m,k,e,f)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO2(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO2(k,j,a,n) = INTERMOOVO2(k,j,a,n) +
     &  INTERMOOVO(m,k,f,j)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO2(k,j,a,n)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO2)
        ENDIF
!
c     -1/6*( ( (T2+{mekf}*T1{em}) * T2{fjan}) * F{ni})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) -
     &          (1.0/6.0)*T2(m,k,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(k,j,a,n) = INTERMOOVO(k,j,a,n) +
     &  INTERMOV(k,f)*T2(j,n,f,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(k,j,a,n)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF

!
c     1/6*( ( (T2+{mekf}*T1{em}) * T2{fian}) * F{nj})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     & (1.0/6.0)*T2(m,k,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(k,i,a,n) = INTERMOOVO(k,i,a,n) +
     &  INTERMOV(k,f)*T2(i,n,f,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  n  = 1, NOS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOVO(k,i,a,n)*FOCK(n,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/6*( ( (F{ag}*T2{gifm}) * T2+{mekf}) * T1{ej})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(a,i,f,m) = INTERMVOVO(a,i,f,m) -
     &          (1.0/6.0)*FOCK(a,g)*T2(i,m,g,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO2(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO2(a,i,e,k) = INTERMVOVO2(a,i,e,k) +
     &  INTERMVOVO(a,i,f,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOVO2(a,i,e,k)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO2)
        ENDIF
!
c     1/6*( ( (F{ag}*T2{gjfm}) * T2+{mekf}) * T1{ei})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(a,j,f,m) = INTERMVOVO(a,j,f,m) +
     & (1.0/6.0)*FOCK(a,g)*T2(j,m,g,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO2(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO2 = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO2(a,j,e,k) = INTERMVOVO2(a,j,e,k) +
     &  INTERMVOVO(a,j,f,m)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOVO2(a,j,e,k)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO2)
        ENDIF
!
c     1/6*( ( (F{nj}*T1{en}) * T2+{mekf}) * T2{aifm})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(j,e) = INTERMOV(j,e) +
     & (1.0/6.0)*FOCK(n,j)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(j,m,k,f) = INTERMOOOV(j,m,k,f) +
     &  INTERMOV(j,e)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(j,m,k,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1/6*( ( (F{ni}*T1{en}) * T2+{mekf}) * T2{ajfm})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(i,e) = INTERMOV(i,e) -
     &          (1.0/6.0)*FOCK(n,i)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(i,m,k,f) = INTERMOOOV(i,m,k,f) +
     &  INTERMOV(i,e)*T2(m,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOOOV(i,m,k,f)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/4*( ( (F{nm}*T1{an}) * T2+{mekf}) * T2{fjei})
! c      Commutator = 3, PT Order = 4
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
!      & (1.0/4.0)*FOCK(n,m)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOV(a,e,k,f) = INTERMVVOV(a,e,k,f) +
!      &  INTERMOV(m,a)*T2(m,k,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVVOV(a,e,k,f)*T2(j,i,f,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
c     -1/6*( ( (T1+{me}*F{am}) * T2{fjei}) * T1+{kf})
c      Commutator = 3, PT Order = 5
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
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOO(a,f,j,i) = INTERMVVOO(a,f,j,i) +
     &  INTERMVV(e,a)*T2(j,i,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVVOO(a,f,j,i)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/2*( ( (F{mf}*T1{am}) * T2{ejfi}) * T1+{ke})
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
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(f,a) = INTERMVV(f,a) -
!     & (1.0/2.0)*FOCK(m,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOO(a,e,j,i) = INTERMVVOO(a,e,j,i) +
!      &  INTERMVV(f,a)*T2(j,i,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVVOO(a,e,j,i)*T1(k,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/6*( ( (T1+{me}*T1{am}) * T2{ejfi}) * F{kf})
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
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      & (1.0/6.0)*T1(m,e)*T1(m,a)
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
!         INTERMVOVO(a,j,f,i) = INTERMVOVO(a,j,f,i) +
!      &  INTERMVV(e,a)*T2(j,i,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVOVO(a,j,f,i)*FOCK(k,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!         ENDIF
!
c     -1/6*( ( (F{kf}*T1{fm}) * T1+{me}) * T2{ejai})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(k,m) = INTERMOO(k,m) -
     &          (1.0/6.0)*FOCK(k,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(k,e) = INTERMOV(k,e) +
     &  INTERMOO(k,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,e)*T2(j,i,e,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/2*( ( (T1+{ke}*T1{em}) * F{mf}) * T2{fjai})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(k,m) = INTERMOO(k,m) -
     &          (1.0/2.0)*T1(k,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(k,f) = INTERMOV(k,f) +
     &  INTERMOO(k,m)*FOCK(m,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,f)*T2(j,i,f,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/2*( ( (T1+{kf}*F{fm}) * T1+{me}) * T2{ejai})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(k,m) = INTERMOO(k,m) -
     &          (1.0/2.0)*T1(k,f)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(k,e) = INTERMOV(k,e) +
     &  INTERMOO(k,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOV(k,e)*T2(j,i,e,a)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (T2+{mekf}*F{am}) * T1{fj}) * T1{ei})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
     &          (1.0/3.0)*T2(m,k,e,f)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(e,k,a,j) = INTERMVOVO(e,k,a,j) +
     &  INTERMVOVV(e,k,f,a)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMVOVO(e,k,a,j)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*T1{ej}) * T2{aifm}) * F{kf})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,j) = INTERMOO(m,j) -
     &          (1.0/6.0)*T1(m,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(j,a,i,f) = INTERMOVOV(j,a,i,f) +
     &  INTERMOO(m,j)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(j,a,i,f)*FOCK(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{me}*T1{ei}) * T2{ajfm}) * F{kf})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) +
     & (1.0/6.0)*T1(m,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(i,a,j,f) = INTERMOVOV(i,a,j,f) +
     &  INTERMOO(m,i)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(i,a,j,f)*FOCK(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/2*( ( (F{mf}*T1{fj}) * T2{aiem}) * T1+{ke})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,j) = INTERMOO(m,j) -
     &          (1.0/2.0)*FOCK(m,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOV(j,a,i,e) = INTERMOVOV(j,a,i,e) +
     &  INTERMOO(m,j)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(j,a,i,e)*T1(k,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     1/2*( ( (F{mf}*T1{fi}) * T2{ajem}) * T1+{ke})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) +
     & (1.0/2.0)*FOCK(m,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOV(i,a,j,e) = INTERMOVOV(i,a,j,e) +
     &  INTERMOO(m,i)*T2(j,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(i,a,j,e)*T1(k,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     -1/3*( ( (T2+{mekf}*T1{am}) * F{fj}) * T1{ei})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) -
!     & (1.0/3.0)*T2(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(e,k,a,j) = INTERMVOVO(e,k,a,j) +
!      &  INTERMVOVV(e,k,f,a)*FOCK(f,j)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVOVO(e,k,a,j)*T1(i,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!         ENDIF
! !
! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
! c     1/3*( ( (T2+{mekf}*T1{am}) * F{fi}) * T1{ej})
! c      Commutator = 3, PT Order = 5
! !
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  k  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,k,f,a) = INTERMVOVV(e,k,f,a) +
!      & (1.0/3.0)*T2(m,k,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVO(e,k,a,i) = INTERMVOVO(e,k,a,i) +
!      &  INTERMVOVV(e,k,f,a)*FOCK(f,i)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  j  = 1, NOS
!         DO  k  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(k,a,j,i) = VT(k,a,j,i) +
!      &  INTERMVOVO(e,k,a,i)*T1(j,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!         ENDIF
!
c     -1/6*( ( (T1+{me}*F{ej}) * T2{aifm}) * T1+{kf})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,j) = INTERMOO(m,j) -
     &          (1.0/6.0)*T1(m,e)*FOCK(e,j)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(j,a,i,f) = INTERMOVOV(j,a,i,f) +
     &  INTERMOO(m,j)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(j,a,i,f)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{me}*F{ei}) * T2{ajfm}) * T1+{kf})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) +
     & (1.0/6.0)*T1(m,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(i,a,j,f) = INTERMOVOV(i,a,j,f) +
     &  INTERMOO(m,i)*T2(j,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,a,j,i) = VT(k,a,j,i) +
     &  INTERMOVOV(i,a,j,f)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF

        DO  a  = NOS+1, NOS+NVS
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        VT(a,k,i,j) = VT(k,a,j,i)
        VT(j,i,k,a) = VT(k,a,j,i)
        VT(i,j,a,k) = VT(k,a,j,i)
        ENDDO; ENDDO; ENDDO; ENDDO

      END
