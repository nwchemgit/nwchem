      SUBROUTINE HP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,HT)
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
      DOUBLE PRECISION :: HT(NOS+NVS,NOS+NVS)                 ! TRANSFORMED 1E- INTEGRALS
C
C     CLUSTER OPERATORS: OCCUPIED IS NUMBERED 1--NOS AND VIRTUAL IS NUMBERED (NOS+1)--(NOS+NVS)
      DOUBLE PRECISION :: T1(NOS,NOS+1:NOS+NVS)                   ! (OCCUPIED,VIRTUAL)
      DOUBLE PRECISION :: T2(NOS,NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS) ! (OCCUPIED,OCCUPIED,VIRTUAL,VIRTUAL)
C
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVO(:,:,:,:)


!
c     1*F{ai}
c      Commutator = 0, PT Order = 0
!
        IF(NONHF)THEN
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        HT(a,i) = HT(a,i) + FOCK(a,i)
        ENDDO; ENDDO
        ENDIF
        ENDIF
!
c     1*F{me}*T2{aiem}
c      Commutator = 1, PT Order = 1
!
        IF(NONHF)THEN
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     & FOCK(m,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
        ENDIF
!
c     1*F{ae}*T1{ei}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     & FOCK(a,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*F{mi}*T1{am}
c      Commutator = 1, PT Order = 2
!
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         HT(a,i) = HT(a,i) -
!     & (1.0)*FOCK(m,i)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ENDIF
!
c     -1/2*V{mine}*T2{amen}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) -
     &          (1.0/2.0)*V(m,n,i,e)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*V{aemf}*T2{eifm}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     & (1.0/2.0)*V(a,m,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*T1+{me}*V{aiem}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     & T1(m,e)*V(a,e,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{aime}*T1{em}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     & V(a,m,i,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*( (T2+{menf}*F{em}) * T2{aifn})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,f) = INTERMOV(n,f) +
     & (1.0/2.0)*T2(m,n,e,f)*FOCK(e,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/4*( (T2+{menf}*F{am}) * T2{eifn})
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
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) -
     &          (1.0/4.0)*T2(m,n,e,f)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVOVV(e,n,f,a)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{amfn}) * F{ei})
c      Commutator = 2, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
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
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(e,a)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
c     -1*( (T1+{me}*F{nm}) * T2{aien})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(e,n) = INTERMVO(e,n) -
     &          (1.0)*T1(m,e)*FOCK(n,m)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(e,n)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1*( (T1+{me}*F{ef}) * T2{aifm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,f) = INTERMOV(m,f) +
     & T1(m,e)*FOCK(e,f)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(m,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T1+{me}*T2{amen}) * F{ni})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
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
        DO  n  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(a,n)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/2*( (T1+{me}*T2{eifm}) * F{af})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
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
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(i,f)*FOCK(a,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( (T2+{menf}*T2{fmgn}) * V{aieg})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVV(e,g) = INTERMVV(e,g) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(e,g)*V(a,e,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1*( (T2+{menf}*V{ameg}) * T2{fign})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVVV(n,f,a,g) = INTERMOVVV(n,f,a,g) +
     & T2(m,n,e,f)*V(a,e,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOVVV(n,f,a,g)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     1/2*( (T2+{menf}*V{emfg}) * T2{aign})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(n,g) = INTERMOV(n,g) +
     & (1.0/2.0)*T2(m,n,e,f)*V(e,f,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,g)*T2(i,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1*( (T2+{menf}*T2{anfo}) * V{eiom})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
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
        INTERMOVVO(m,e,a,o) = INTERMOVVO(m,e,a,o) -
     &          (1.0)*T2(m,n,e,f)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOVVO(m,e,a,o)*V(e,o,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1/2*( (T2+{menf}*V{emon}) * T2{aifo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,o) = INTERMVO(f,o) -
     &          (1.0/2.0)*T2(m,n,e,f)*V(e,o,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(f,o)*T2(i,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/4*( (T2+{menf}*V{eifg}) * T2{amgn}) ! not the most efficient but avoind VVVV intermediate
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(m,n,i,g) = INTERMOOOV(m,n,i,g) -
     &          (1.0/4.0)*T2(m,n,e,f)*V(e,f,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOOOV(m,n,i,g)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     1/2*( (T2+{menf}*T2{enfo}) * V{aiom})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,o) = INTERMOO(m,o) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMOO(m,o)*V(a,o,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/4*( (T2+{menf}*V{amon}) * T2{eifo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,1:NOS))
        INTERMVVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVO(e,f,a,o) = INTERMVVVO(e,f,a,o) +
     & (1.0/4.0)*T2(m,n,e,f)*V(a,o,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVVVO(e,f,a,o)*T2(i,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{me}*T1{am}) * F{ei})
c      Commutator = 2, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
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
!         DO  e  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVV(e,a)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
!
c     -1/2*( (T1+{me}*F{am}) * T1{ei})
c      Commutator = 2, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) -
     &          (1.0/2.0)*T1(m,e)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(e,a)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (F{me}*T1{am}) * T1{ei})
c      Commutator = 2, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0)*FOCK(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVV(e,a)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
!
c     -1*( (T1+{me}*V{amnf}) * T2{eifn})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,
     &  NOS+1:NOS+NVS))
        INTERMVVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(e,a,n,f) = INTERMVVOV(e,a,n,f) -
     &          (1.0)*T1(m,e)*V(a,n,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVVOV(e,a,n,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
!
c     -1*( (T1+{me}*T2{amfn}) * V{einf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,1:NOS))
        INTERMVVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVO(e,a,f,n) = INTERMVVVO(e,a,f,n) -
     &          (1.0)*T1(m,e)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVVVO(e,a,f,n)*V(e,n,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
!
c     1*( (T1+{me}*V{emnf}) * T2{aifn})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,f) = INTERMOV(n,f) +
     & T1(m,e)*V(e,n,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T1+{me}*T2{aneo}) * V{niom})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,a,n,o) = INTERMOVOO(m,a,n,o) +
     & (1.0/2.0)*T1(m,e)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMOVOO(m,a,n,o)*V(n,o,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T2+{menf}*T1{fn}) * V{aiem})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) +
     & (1.0/2.0)*T2(m,n,e,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(m,e)*V(a,e,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T1+{me}*V{afeg}) * T2{figm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVVV(m,a,f,g) = INTERMOVVV(m,a,f,g) +
     & (1.0/2.0)*T1(m,e)*V(a,e,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOVVV(m,a,f,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
!
c     1/2*( (T1+{me}*T2{emfn}) * V{ainf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,n) = INTERMVO(f,n) +
     & (1.0/2.0)*T1(m,e)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(f,n)*V(a,n,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1*( (V{menf}*T1{em}) * T2{aifn})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,f) = INTERMOV(n,f) +
     & V(m,n,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/4*( (T2+{menf}*V{amen}) * T1{fi})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
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
     & (1.0/4.0)*T2(m,n,e,f)*V(a,e,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(f,a)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( (T2+{menf}*T1{an}) * V{eifm})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
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
!      & (1.0/4.0)*T2(m,n,e,f)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMOVVV(m,e,f,a)*V(e,f,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
c     -1/2*( (V{menf}*T2{amfn}) * T1{ei})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
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
     &          (1.0/2.0)*V(m,n,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(e,a)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (V{menf}*T1{am}) * T2{eifn})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
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
!     & (1.0/2.0)*V(m,n,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVOVV(e,n,f,a)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!
c     1*( (V{aemf}*T1{fm}) * T1{ei})
c      Commutator = 2, PT Order = 5
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
        INTERMVV(a,e) = INTERMVV(a,e) +
     & V(a,m,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1*( (T1+{me}*T1{fm}) * V{aief})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,f) = INTERMVV(e,f) +
     & T1(m,e)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(e,f)*V(a,e,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1*( (T1+{me}*V{amef}) * T1{fi})
c      Commutator = 2, PT Order = 5
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
        INTERMVV(a,f) = INTERMVV(a,f) -
     &          (1.0)*T1(m,e)*V(a,e,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (V{mine}*T1{en}) * T1{am})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,i) = INTERMOO(m,i) -
!     & (1.0)*V(m,n,i,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         HT(a,i) = HT(a,i) +
!      &  INTERMOO(m,i)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -1*( (T1+{me}*T1{en}) * V{ainm})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,n) = INTERMOO(m,n) -
     &          (1.0)*T1(m,e)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMOO(m,n)*V(a,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{me}*V{einm}) * T1{an})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(i,n) = INTERMOO(i,n) +
!      & T1(m,e)*V(e,n,i,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         HT(a,i) = HT(a,i) +
!      &  INTERMOO(i,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!
c     -2/3*( ( (T2+{menf}*T2{anfo}) * F{og}) * T2{eigm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
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
        INTERMOVVO(m,e,a,o) = INTERMOVVO(m,e,a,o) -
     &          (2.0/3.0)*T2(m,n,e,f)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVVV(m,e,a,g) = INTERMOVVV(m,e,a,g) +
     &  INTERMOVVO(m,e,a,o)*FOCK(o,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOVVV(m,e,a,g)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{og}*T2{fngo}) * T2+{menf}) * T2{aiem})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(f,n) = INTERMVO(f,n) +
     & (1.0/6.0)*FOCK(o,g)*T2(n,o,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) +
     &  INTERMVO(f,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(m,e)*T2(i,m,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     1/12*( ( (F{og}*T2{ango}) * T2+{menf}) * T2{eifm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(a,n) = INTERMVO(a,n) +
     & (1.0/12.0)*FOCK(o,g)*T2(n,o,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO
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
        INTERMVOVV(a,m,e,f) = INTERMVOVV(a,m,e,f) +
     &  INTERMVO(a,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVOVV(a,m,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     1/3*( ( (T2+{menf}*T2{enfo}) * F{og}) * T2{aigm})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,o) = INTERMOO(m,o) +
     & (1.0/3.0)*T2(m,n,e,f)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(m,g) = INTERMOV(m,g) +
     &  INTERMOO(m,o)*FOCK(o,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(m,g)*T2(i,m,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     1/3*( ( (F{og}*T2{fmgn}) * T2+{menf}) * T2{aieo})
c      Commutator = 3, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/4, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS

        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(o,f,m,n) = INTERMOVOO(o,f,m,n) +
     & (1.0/3.0)*FOCK(o,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(o,e) = INTERMOV(o,e) +
     &  INTERMOVOO(o,f,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(o,e)*T2(i,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{amfn}) * T2{eigo}) * F{og})
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
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(a,i,g,o) = INTERMVOVO(a,i,g,o) +
     &  INTERMVV(e,a)*T2(i,o,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVOVO(a,i,g,o)*FOCK(o,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
        ENDIF
!
c     1/6*( ( (F{og}*T2{amgn}) * T2+{menf}) * T2{eifo})
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
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(o,a,m,n) = INTERMOVOO(o,a,m,n) +
     & (1.0/6.0)*FOCK(o,g)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMOVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVVV(o,a,e,f) = INTERMOVVV(o,a,e,f) +
     &  INTERMOVOO(o,a,m,n)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOVVV(o,a,e,f)*T2(i,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVV)
        ENDIF
        ENDIF
!
c     -1/2*( ( (T2+{menf}*T1{fm}) * F{eg}) * T2{aign})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,n) = INTERMVO(e,n) -
     &          (1.0/2.0)*T2(m,n,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(n,g) = INTERMOV(n,g) +
     &  INTERMVO(e,n)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,g)*T2(i,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( ( (T2+{menf}*T1{en}) * F{om}) * T2{aifo})
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
     & (1.0/2.0)*T2(m,n,e,f)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,o) = INTERMVO(f,o) +
     &  INTERMOV(m,f)*FOCK(o,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(f,o)*T2(i,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( ( (T2+{menf}*T1{an}) * F{om}) * T2{eifo})
c      Commutator = 3, PT Order = 4
!
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
!         INTERMOVVV(m,e,f,a) = INTERMOVVV(m,e,f,a) -
!     & (1.0/4.0)*T2(m,n,e,f)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVO(e,f,a,o) = INTERMVVVO(e,f,a,o) +
!      &  INTERMOVVV(m,e,f,a)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVVVO(e,f,a,o)*T2(i,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
c     1/4*( ( (F{eg}*T2{amgn}) * T2+{menf}) * T1{fi})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(e,a,m,n) = INTERMVVOO(e,a,m,n) +
     & (1.0/4.0)*FOCK(e,g)*T2(m,n,a,g)
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
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/6*( ( (T2+{menf}*T1{em}) * T2{fign}) * F{ag})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,f) = INTERMOV(n,f) -
     &          (1.0/6.0)*T2(m,n,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV2(1:NOS,NOS+1:NOS+NVS))
        INTERMOV2 = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV2(i,g) = INTERMOV2(i,g) +
     &  INTERMOV(n,f)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV2(i,g)*FOCK(a,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV2)
        ENDIF
!
c     1/6*( ( (T2+{menf}*T1{em}) * T2{anfo}) * F{oi})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,f) = INTERMOV(n,f) +
     & (1.0/6.0)*T2(m,n,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(a,o) = INTERMVO(a,o) +
     &  INTERMOV(n,f)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(a,o)*FOCK(o,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T2+{menf}*T1{am}) * F{eg}) * T2{fign})
c      Commutator = 3, PT Order = 4
!
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
!      & (1.0/2.0)*T2(m,n,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         INTERMOVVV(n,f,a,g) = INTERMOVVV(n,f,a,g) +
!      &  INTERMVOVV(e,n,f,a)*FOCK(e,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         DO  g  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMOVVV(n,f,a,g)*T2(i,n,f,g)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
c     -1/2*( ( (F{om}*T2{anfo}) * T2+{menf}) * T1{ei})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,a,n,f) = INTERMOVOV(m,a,n,f) -
     &          (1.0/2.0)*FOCK(o,m)*T2(n,o,a,f)
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
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/2*( ( (F{eg}*T1{gm}) * T2+{menf}) * T2{aifn})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,m) = INTERMVO(e,m) +
     & (1.0/2.0)*FOCK(e,g)*T1(m,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,f) = INTERMOV(n,f) +
     &  INTERMVO(e,m)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( ( (F{om}*T1{eo}) * T2+{menf}) * T2{aifn})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) -
     &          (1.0/2.0)*FOCK(o,m)*T1(o,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV2(1:NOS,NOS+1:NOS+NVS))
        INTERMOV2 = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV2(n,f) = INTERMOV2(n,f) +
     &  INTERMOV(m,e)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV2(n,f)*T2(i,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV2)
        ENDIF
!
c     -1/12*( ( (F{ag}*T1{gm}) * T2+{menf}) * T2{eifn})
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
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVOV(a,e,n,f) = INTERMVVOV(a,e,n,f) +
     &  INTERMVO(a,m)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVVOV(a,e,n,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOV)
        ENDIF
!
c     1/12*( ( (F{ag}*T2{fmgn}) * T2+{menf}) * T1{ei})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,f,m,n) = INTERMVVOO(a,f,m,n) +
     & (1.0/12.0)*FOCK(a,g)*T2(m,n,f,g)
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
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/12*( ( (T2+{menf}*T1{am}) * T2{enfo}) * F{oi})
c      Commutator = 3, PT Order = 4
!
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
!     & (1.0/12.0)*T2(m,n,e,f)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(a,o) = INTERMVO(a,o) +
!      &  INTERMVOVV(e,n,f,a)*T2(n,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  o  = 1, NOS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVO(a,o)*FOCK(o,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
c     1/12*( ( (T2+{menf}*T2{amfn}) * T1{eo}) * F{oi})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
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
     & (1.0/12.0)*T2(m,n,e,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(a,o) = INTERMVO(a,o) +
     &  INTERMVV(e,a)*T1(o,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  o  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMVO(a,o)*FOCK(o,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( ( (F{om}*T1{ao}) * T2+{menf}) * T2{eifn})
c      Commutator = 3, PT Order = 4
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOV(m,a) = INTERMOV(m,a) +
!      & (1.0/4.0)*FOCK(o,m)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVOV(a,e,n,f) = INTERMVVOV(a,e,n,f) +
!      &  INTERMOV(m,a)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVVOV(a,e,n,f)*T2(i,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
c     -1/4*( ( (T2+{menf}*T2{amfn}) * F{eg}) * T1{gi})
c      Commutator = 3, PT Order = 4
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
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
        DO  i  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV2(a,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV2)
        ENDIF
!
c     -2/3*( ( (F{nf}*T1{fm}) * T1+{me}) * T2{aien})
c      Commutator = 3, PT Order = 5
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
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(n,e) = INTERMOV(n,e) +
     &  INTERMOO(n,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,e)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T1+{me}*T2{eifm}) * F{nf}) * T1{an})
c      Commutator = 3, PT Order = 5
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(i,f) = INTERMOV(i,f) +
!      & (1.0/2.0)*T1(m,e)*T2(i,m,e,f)
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
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         HT(a,i) = HT(a,i) +
!      &  INTERMOO(i,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
c     -2/3*( ( (T1+{nf}*F{fm}) * T1+{me}) * T2{aien})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,m) = INTERMOO(n,m) -
     &          (2.0/3.0)*T1(n,f)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(n,e) = INTERMOV(n,e) +
     &  INTERMOO(n,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(n,e)*T2(i,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T1+{me}*T1{en}) * F{nf}) * T2{aifm})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,n) = INTERMOO(m,n) -
     &          (2.0/3.0)*T1(m,e)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,f) = INTERMOV(m,f) +
     &  INTERMOO(m,n)*FOCK(n,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMOV(m,f)*T2(i,m,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{nf}*T2{amfn}) * T1+{me}) * T1{ei})
c      Commutator = 3, PT Order = 5
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
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{nf}*T2{amfn}) * T1+{me}) * F{ei})
c      Commutator = 3, PT Order = 5
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
     &          (1.0/6.0)*T1(n,f)*T2(m,n,a,f)
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
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{nf}*T2{eifn}) * T1+{me}) * F{am})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) -
     &          (1.0/6.0)*T1(n,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(i,m) = INTERMOO(i,m) +
     &  INTERMVO(e,i)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        HT(a,i) = HT(a,i) +
     &  INTERMOO(i,m)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (F{nf}*T2{eifn}) * T1+{me}) * T1{am})
c      Commutator = 3, PT Order = 5
!
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
!         DO  i  = 1, NOS
!         DO  m  = 1, NOS
!         HT(a,i) = HT(a,i) +
!      &  INTERMOO(i,m)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/3*( ( (T2+{menf}*T1{fn}) * T1{am}) * F{ei})
c      Commutator = 3, PT Order = 5
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(m,e) = INTERMOV(m,e) -
!     & (1.0/3.0)*T2(m,n,e,f)*T1(n,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      &  INTERMOV(m,e)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVV(e,a)*FOCK(e,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
!
c     1/2*( ( (T1+{me}*T2{amen}) * F{nf}) * T1{fi})
c      Commutator = 3, PT Order = 5
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
        INTERMVO(a,n) = INTERMVO(a,n) +
     & (1.0/2.0)*T1(m,e)*T2(m,n,a,e)
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
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (T2+{menf}*T1{fn}) * F{am}) * T1{ei})
c      Commutator = 3, PT Order = 5
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) -
     &          (1.0/3.0)*T2(m,n,e,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVV(e,a) = INTERMVV(e,a) +
     &  INTERMOV(m,e)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(e,a)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/3*( ( (T2+{menf}*F{em}) * T1{an}) * T1{fi})
c      Commutator = 3, PT Order = 5
!
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
!         DO  i  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVV(f,a)*T1(i,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!         ENDIF
!
c     -1/6*( ( (F{af}*T1{fm}) * T1+{me}) * T1{ei})
c      Commutator = 3, PT Order = 6
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
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,i) = HT(a,i) +
     &  INTERMVV(a,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (F{nm}*T1{an}) * T1+{me}) * T1{ei})
c      Commutator = 3, PT Order = 6
!
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
!         DO  i  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVV(a,e)*T1(i,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( ( (T1+{me}*T1{am}) * F{ef}) * T1{fi})
c      Commutator = 3, PT Order = 6
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
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
!         DO  i  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVV2(a,f)*T1(i,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV2)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{me}*T1{am}) * T1{en}) * F{ni})
c      Commutator = 3, PT Order = 6
!
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
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,n) = INTERMVO(a,n) +
!      &  INTERMVV(e,a)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  i  = 1, NOS
!         DO  n  = 1, NOS
!         HT(a,i) = HT(a,i) +
!      &  INTERMVO(a,n)*FOCK(n,i)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF

      DO  a  = NOS+1, NOS+NVS
        DO  i  = 1, NOS
          HT(i,a) = HT(a,i)
      ENDDO; ENDDO

      END
