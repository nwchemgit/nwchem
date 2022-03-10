      SUBROUTINE PP(NOS,NVS,NVAS,NACTV,T1,T2,FOCK,V,DUCC_MODEL,NONHF,HT)
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
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVO(:,:,:,:)

!
c     1*F{ab}
c      Commutator = 0, PT Order = 0
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        HT(a,b) = HT(a,b) + FOCK(a,b)
        ENDDO; ENDDO
        ENDIF
!
c     1/2*V{menb}*T2{amen}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     & (1.0/2.0)*V(m,n,e,b)*T2(m,n,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*T2+{menb}*V{amen}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     & (1.0/2.0)*T2(m,n,e,b)*V(a,e,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*T1+{mb}*F{am}
c      Commutator = 1, PT Order = 2
!
!         IF(NONHF)THEN
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      & -(1.0)*T1(m,b)*FOCK(a,m)
!         ENDDO; ENDDO; ENDDO
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*F{mb}*T1{am}
c      Commutator = 1, PT Order = 2
!
!         IF(NONHF)THEN
!         IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      & -(1.0)*FOCK(m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         ENDIF
!         ENDIF
!
c     -1*V{aemb}*T1{em}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) -
     &          (1.0)*V(a,m,e,b)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*T1+{me}*V{ameb}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) -
     &          (1.0)*T1(m,e)*V(a,e,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*( (T2+{menb}*F{ef}) * T2{amfn})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVV(m,n,b,f) = INTERMOOVV(m,n,b,f) +
     & (1.0/2.0)*T2(m,n,e,b)*FOCK(e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOOVV(m,n,b,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1*( (T2+{menb}*F{om}) * T2{aneo})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVOVO(e,n,b,o) = INTERMVOVO(e,n,b,o) +
     & T2(m,n,e,b)*FOCK(o,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVOVO(e,n,b,o)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1/4*( (T2+{menb}*T2{emfn}) * F{af})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
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
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVV(b,f)*FOCK(a,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{amfn}) * F{eb})
c      Commutator = 2, PT Order = 2
!
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
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVV(e,a)*FOCK(e,b)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1*( (T2+{menb}*V{emof}) * T2{anfo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(n,b,o,f) = INTERMOVOV(n,b,o,f) -
     &          (1.0)*T2(m,n,e,b)*V(e,o,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOVOV(n,b,o,f)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1/2*( (T2+{menf}*T2{fmgn}) * V{ageb})
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
        INTERMVV(e,g) = INTERMVV(e,g) +
     & (1.0/2.0)*T2(m,n,e,f)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVV(e,g)*V(a,e,g,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1*( (T2+{menf}*V{emob}) * T2{anfo})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(n,f,o,b) = INTERMOVOV(n,f,o,b) +
     & T2(m,n,e,f)*V(e,o,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOVOV(n,f,o,b)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1*( (T2+{menb}*T2{enfo}) * V{amof})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
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
        INTERMOVVO(m,b,f,o) = INTERMOVVO(m,b,f,o) +
     & T2(m,n,e,b)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOVVO(m,b,f,o)*V(a,o,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1/4*( (T2+{menb}*V{ompn}) * T2{aoep})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  p  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVVOO(e,b,o,p) = INTERMVVOO(e,b,o,p) +
     & (1.0/4.0)*T2(m,n,e,b)*V(o,p,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  o  = 1, NOS
        DO  p  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVVOO(e,b,o,p)*T2(o,p,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     1/4*( (V{afeg}*T2{fmgn}) * T2+{menb})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMVVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVVOO(a,e,m,n) = INTERMVVOO(a,e,m,n) +
     & (1.0/4.0)*V(a,e,f,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVVOO(a,e,m,n)*T2(m,n,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVOO)
        ENDIF
!
c     -1/2*( (T2+{menf}*T2{enfo}) * V{amob})
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
        INTERMOO(m,o) = INTERMOO(m,o) -
     &          (1.0/2.0)*T2(m,n,e,f)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMOO(m,o)*V(a,o,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/4*( (T2+{menf}*V{egfb}) * T2{amgn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVV(m,n,g,b) = INTERMOOVV(m,n,g,b) +
     & (1.0/4.0)*T2(m,n,e,f)*V(e,f,g,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOOVV(m,n,g,b)*T2(m,n,a,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T2+{menb}*F{em}) * T1{an})
c      Commutator = 2, PT Order = 3
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,b) = INTERMOV(n,b) +
!      & -(1.0/2.0)*T2(m,n,e,b)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(n,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!         ENDIF
!
c     1/2*( (T2+{menb}*T1{en}) * F{am})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
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
        DO  m  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMOV(m,b)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (F{ne}*T2{amen}) * T1+{mb})
c      Commutator = 2, PT Order = 3
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,m) = INTERMVO(a,m) +
!      & -(1.0/2.0)*FOCK(n,e)*T2(m,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
c     1/2*( (T1+{me}*T2{amen}) * F{nb})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
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
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMVO(a,n)*FOCK(n,b)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{mb}*F{nm}) * T1{an})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVO(b,n) = INTERMVO(b,n) +
!      & T1(m,b)*FOCK(n,m)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{me}*F{eb}) * T1{am})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,b) = INTERMOV(m,b) +
!      & -(1.0/2.0)*T1(m,e)*FOCK(e,b)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (F{ae}*T1{em}) * T1+{mb})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,m) = INTERMVO(a,m) +
!      & -(1.0/2.0)*FOCK(a,e)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
c     -1/2*( (T2+{menf}*T1{fn}) * V{ameb})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) -
     &          (1.0/2.0)*T2(m,n,e,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOV(m,e)*V(a,e,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( (T2+{menb}*T1{eo}) * V{amon})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,b,o) = INTERMOOVO(m,n,b,o) -
     &          (1.0/2.0)*T2(m,n,e,b)*T1(o,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMOOVO(m,n,b,o)*V(a,o,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1*( (V{amef}*T1{fn}) * T2+{menb})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(a,m,e,n) = INTERMVOVO(a,m,e,n) +
     & V(a,e,m,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVOVO(a,m,e,n)*T2(m,n,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1*( (T1+{me}*V{efnb}) * T2{amfn})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,f,n,b) = INTERMOVOV(m,f,n,b) +
     & T1(m,e)*V(e,n,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOVOV(m,f,n,b)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/2*( (T1+{me}*T2{emfn}) * V{afnb})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,n) = INTERMVO(f,n) -
     &          (1.0/2.0)*T1(m,e)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVO(f,n)*V(a,n,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T2+{menb}*V{emon}) * T1{ao})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(b,o) = INTERMVO(b,o) +
!      & (1.0/2.0)*T2(m,n,e,b)*V(e,o,m,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,o)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( (T2+{menf}*V{emfb}) * T1{an})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,b) = INTERMOV(n,b) +
!      & -(1.0/4.0)*T2(m,n,e,f)*V(e,f,m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(n,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (V{nmoe}*T2{aneo}) * T1+{mb})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,a) = INTERMOV(m,a) +
!      & (1.0/2.0)*V(n,o,m,e)*T2(n,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(m,a)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
c     -1/2*( (T1+{me}*T2{aneo}) * V{nmob})
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
        INTERMOVOO(m,a,n,o) = INTERMOVOO(m,a,n,o) -
     &          (1.0/2.0)*T1(m,e)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMOVOO(m,a,n,o)*V(n,o,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( (V{aenf}*T2{emfn}) * T1+{mb})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(a,m) = INTERMVO(a,m) +
!      & -(1.0/4.0)*V(a,n,e,f)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
c     -1*( (T1+{me}*T1{fm}) * V{afeb})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,f) = INTERMVV(e,f) -
     &          (1.0)*T1(m,e)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVV(e,f)*V(a,e,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (V{menb}*T1{en}) * T1{am})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,b) = INTERMOV(m,b) +
!      & V(m,n,e,b)*T1(n,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{me}*V{amen}) * T1+{nb})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,n) = INTERMVO(a,n) +
!      & T1(m,e)*V(a,e,m,n)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
c     1*( (T1+{me}*T1{en}) * V{amnb})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,n) = INTERMOO(m,n) +
     & T1(m,e)*T1(n,e)
        ENDDO; ENDDO; ENDDO
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMOO(m,n)*V(a,n,m,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (V{amne}*T1{en}) * T1+{mb})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,m) = INTERMVO(a,m) +
!      & -(1.0)*V(a,n,m,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE

!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (T1+{me}*V{emnb}) * T1{an})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,b) = INTERMOV(n,b) +
!      & -(1.0)*T1(m,e)*V(e,n,m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(n,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
c     -1/6*( ( (T2+{menb}*T1{em}) * T2{anfo}) * F{of})
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
        ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,1:NOS))
        INTERMVVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVO(b,a,f,o) = INTERMVVVO(b,a,f,o) +
     &  INTERMOV(n,b)*T2(n,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVVVO(b,a,f,o)*FOCK(o,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVO)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{me}*T2{ameo}) * T2+{nfob}) * F{fn})
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
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(a,o) = INTERMVO(a,o) +
     & (1.0/6.0)*T1(m,e)*T2(m,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS))
        INTERMVOVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVV(a,n,f,b) = INTERMVOVV(a,n,f,b) +
     &  INTERMVO(a,o)*T2(n,o,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVOVV(a,n,f,b)*FOCK(f,n)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (T2+{menb}*T1{eo}) * F{of}) * T2{amfn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,b,o) = INTERMOOVO(m,n,b,o) -
     &          (1.0/3.0)*T2(m,n,e,b)*T1(o,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVV(m,n,b,f) = INTERMOOVV(m,n,b,f) +
     &  INTERMOOVO(m,n,b,o)*FOCK(o,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOOVV(m,n,b,f)*T2(m,n,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
        ENDIF
!
c     -1/3*( ( (T1+{me}*T2{aneo}) * F{fm}) * T2+{nfob})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOO(m,a,n,o) = INTERMOVOO(m,a,n,o) -
     &          (1.0/3.0)*T1(m,e)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOOV(a,n,o,f) = INTERMVOOV(a,n,o,f) +
     &  INTERMOVOO(m,a,n,o)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVOOV(a,n,o,f)*T2(n,o,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T2+{neof}*F{en}) * T2{amfo}) * T1+{mb})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(o,f) = INTERMOV(o,f) +
!      & -(1.0/6.0)*T2(n,o,e,f)*FOCK(e,n)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(a,m) = INTERMVO(a,m) +
!      &  INTERMOV(o,f)*T2(m,o,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
c     1/6*( ( (T1+{me}*T2{emfo}) * T2+{nfob}) * F{an})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,o) = INTERMVO(f,o) +
     & (1.0/6.0)*T1(m,e)*T2(m,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(n,b) = INTERMOV(n,b) +
     &  INTERMVO(f,o)*T2(n,o,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMOV(n,b)*FOCK(a,n)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T1+{me}*F{en}) * T2+{nfob}) * T2{amfo})
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
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(m,f,o,b) = INTERMOVOV(m,f,o,b) +
     &  INTERMOO(m,n)*T2(n,o,f,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOVOV(m,f,o,b)*T2(m,o,a,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     2/3*( ( (F{of}*T1{fm}) * T2+{menb}) * T2{aneo})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(o,m) = INTERMOO(o,m) +
     & (2.0/3.0)*FOCK(o,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVOV(o,e,n,b) = INTERMOVOV(o,e,n,b) +
     &  INTERMOO(o,m)*T2(m,n,e,b)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMOVOV(o,e,n,b)*T2(n,o,a,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{nfob}*T2{enfo}) * T1+{me}) * F{am})
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
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,e) = INTERMVV(b,e) -
     &          (1.0/12.0)*T2(n,o,f,b)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(b,m) = INTERMVO(b,m) +
     &  INTERMVV(b,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMVO(b,m)*FOCK(a,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( ( (T2+{neof}*T2{anfo}) * F{em}) * T1+{mb})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) +
!      & (1.0/4.0)*T2(n,o,e,f)*T2(n,o,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,m) = INTERMVO(a,m) +
!      &  INTERMVV(e,a)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
c     1/6*( ( (T2+{menf}*T1{em}) * T2{anfo}) * F{ob})
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
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  o  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMVO(a,o)*FOCK(o,b)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (F{of}*T2{enfo}) * T2+{menb}) * T1{am})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
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
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,b) = INTERMOV(m,b) +
!      &  INTERMVO(e,n)*T2(m,n,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/12*( ( (T1+{mb}*T2{emfo}) * T2+{neof}) * F{an})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVVVO(b,e,f,o) = INTERMVVVO(b,e,f,o) +
!      & (1.0/12.0)*T1(m,b)*T2(m,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(b,n) = INTERMVO(b,n) +
!      &  INTERMVVVO(b,e,f,o)*T2(n,o,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,n)*FOCK(a,n)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
c     1/12*( ( (T2+{menf}*T2{amfn}) * T1{eo}) * F{ob})
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
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  o  = 1, NOS
        HT(a,b) = HT(a,b) +
     &  INTERMVO(a,o)*FOCK(o,b)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( ( (T2+{menb}*T2{emfn}) * F{of}) * T1{ao})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVV(b,f) = INTERMVV(b,f) +
!      & (1.0/4.0)*T2(m,n,e,b)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(b,o) = INTERMVO(b,o) +
!      &  INTERMVV(b,f)*FOCK(o,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,o)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/12*( ( (T2+{menf}*T1{am}) * T2{enfo}) * F{ob})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,
!      &  NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,f,a) = INTERMVOVV(e,n,f,a) +
!      & -(1.0/12.0)*T2(m,n,e,f)*T1(m,a)
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
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,o)*FOCK(o,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T1+{me}*F{ef}) * T2{amfn}) * T1+{nb})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(m,f) = INTERMOV(m,f) +
!      & (1.0/2.0)*T1(m,e)*FOCK(e,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(a,n) = INTERMVO(a,n) +
!      &  INTERMOV(m,f)*T2(m,n,a,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (F{ef}*T1{fn}) * T2+{menb}) * T1{am})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(e,n) = INTERMVO(e,n) +
!      & (1.0/2.0)*FOCK(e,f)*T1(n,f)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,b) = INTERMOV(m,b) +
!      &  INTERMVO(e,n)*T2(m,n,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (F{om}*T1{eo}) * T2+{menb}) * T1{an})
c      Commutator = 3, PT Order = 5
!
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
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV2(n,b) = INTERMOV2(n,b) +
!      &  INTERMOV(m,e)*T2(m,n,e,b)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV2(n,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV2)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( ( (T2+{menb}*T1{en}) * F{om}) * T1{ao})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(m,b) = INTERMOV(m,b) +
!      & -(1.0/2.0)*T2(m,n,e,b)*T1(n,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVO(b,o) = INTERMVO(b,o) +
!      &  INTERMOV(m,b)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,o)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( ( (T1+{ne}*T2{aneo}) * F{om}) * T1+{mb})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,o) = INTERMVO(a,o) +
!      & -(1.0/2.0)*T1(n,e)*T2(n,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO2(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO2 = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVO2(a,m) = INTERMVO2(a,m) +
!      &  INTERMVO(a,o)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO2(a,m)*T1(m,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO2)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/3*( ( (T2+{menf}*T1{fn}) * F{eb}) * T1{am})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(m,e) = INTERMOV(m,e) +
!      & -(1.0/3.0)*T2(m,n,e,f)*T1(n,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOV2(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV2 = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV2(m,b) = INTERMOV2(m,b) +
!      &  INTERMOV(m,e)*FOCK(e,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV2(m,b)*T1(m,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV2)
!         ENDIF
!
c     -1/6*( ( (T2+{menb}*T1{em}) * T1{fn}) * F{af})
c      Commutator = 3, PT Order = 5
!
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
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(b,f) = INTERMVV(b,f) +
     &  INTERMOV(n,b)*T1(n,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  f  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVV(b,f)*FOCK(a,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T1+{me}*F{om}) * T2{aneo}) * T1+{nb})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(e,o) = INTERMVO(e,o) +
!      & (1.0/2.0)*T1(m,e)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO2(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO2 = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO2(a,n) = INTERMVO2(a,n) +
!      &  INTERMVO(e,o)*T2(n,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO2(a,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO2)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/3*( ( (T1+{me}*T2{emfn}) * F{af}) * T1+{nb})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMVO(f,n) = INTERMVO(f,n) +
!      & -(1.0/3.0)*T1(m,e)*T2(m,n,e,f)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOV(n,a) = INTERMOV(n,a) +
!      &  INTERMVO(f,n)*FOCK(a,f)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMOV(n,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOV)
!         ENDIF
!
c     -1/6*( ( (T1+{nf}*T2{amfn}) * T1+{me}) * F{eb})
c      Commutator = 3, PT Order = 5
!
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
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        HT(a,b) = HT(a,b) +
     &  INTERMVV(a,e)*FOCK(e,b)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T1+{mb}*F{em}) * T1+{ne}) * T1{an})
c      Commutator = 3, PT Order = 6
!
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
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(b,n) = INTERMVO(b,n) +
!      &  INTERMVV(b,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{me}*F{am}) * T1{en}) * T1+{nb})
c      Commutator = 3, PT Order = 6
!
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
!      & (1.0/6.0)*T1(m,e)*FOCK(a,m)
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
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{me}*T1{am}) * T1{en}) * F{nb})
c      Commutator = 3, PT Order = 6
!
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
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(a,n)*FOCK(n,b)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T1+{mb}*T1{em}) * F{ne}) * T1{an})
c      Commutator = 3, PT Order = 6
!
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
!      & (1.0/2.0)*T1(m,b)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(b,n) = INTERMVO(b,n) +
!      &  INTERMVV(b,e)*FOCK(n,e)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         HT(a,b) = HT(a,b) +
!      &  INTERMVO(b,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVO)
!         ENDIF
!         ENDIF


      END
