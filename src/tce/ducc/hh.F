      SUBROUTINE HH(NOS,NVS,T1,T2,FOCK,V,DUCC_MODEL,NONHF,HT)
#ifndef OPTIMIZE_DUCC
CDIR$ OPTIMIZE: 2
#endif
      IMPLICIT NONE
      INTEGER          :: A, B, C, D   ! FIXED VIRTUAL SPIN ORBITAL INDICIES
      INTEGER          :: E, F, G, H   ! FREE VIRTUAL SPIN ORBITAL INDICIES
      INTEGER          :: I, J, K, L   ! FIXED OCCUPIED SPIN ORBITAL INDICIES
      INTEGER          :: M, N, O, P   ! FREE OCCUPIED SPIN ORBITAL INDICIES
      INTEGER          :: NOS, NVS     ! # OF OCCUPIED/VIRTUAL SPINORBITALS
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
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
C
!
c     1*F{ji}
c      Commutator = 0, PT Order = 0
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        HT(j,i) = HT(j,i) + FOCK(j,i)
        ENDDO; ENDDO
        ENDIF
!
c     -1/2*V{mejf}*T2{eifm}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) -
     &  (1.0/2.0)*V(m,j,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1/2*T2+{mejf}*V{eifm}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) -
     &  (1.0/2.0)*T2(m,j,e,f)*V(e,f,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*T1+{je}*F{ei}
c      Commutator = 1, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     & T1(j,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        ENDIF
        ENDIF
!
c     1*F{je}*T1{ei}
c      Commutator = 1, PT Order = 2
!
        IF(NONHF)THEN
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     & FOCK(j,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        ENDIF
        ENDIF
!
c     -1*V{mije}*T1{em}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) -
     &  (1.0)*V(m,j,i,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*T1+{me}*V{eijm}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) -
     &  (1.0)*T1(m,e)*V(e,j,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*( (T2+{mejf}*F{nm}) * T2{eifn})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVOVO(e,j,f,n) = INTERMVOVO(e,j,f,n) +
     & (1.0/2.0)*T2(m,j,e,f)*FOCK(n,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVOVO(e,j,f,n)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1*( (T2+{mejf}*F{eg}) * T2{figm})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,f,g) = INTERMOOVV(m,j,f,g) +
     & T2(m,j,e,f)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOOVV(m,j,f,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{eifn}) * F{jm})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/4.0)*T2(m,n,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(m,i)*FOCK(j,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{mejf}*T2{emfn}) * F{ni})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(j,n) = INTERMOO(j,n) -
     &  (1.0/4.0)*T2(m,j,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(j,n)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1*( (T2+{mejf}*V{emng}) * T2{fign})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(j,f,n,g) = INTERMOVOV(j,f,n,g) +
     & T2(m,j,e,f)*V(e,n,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOV(j,f,n,g)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1*( (T2+{mejf}*T2{fmgn}) * V{eing})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(e,j,g,n) = INTERMVOVO(e,j,g,n) -
     &  (1.0)*T2(m,j,e,f)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVOVO(e,j,g,n)*V(e,n,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1*( (T2+{menf}*V{emjg}) * T2{fign})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(n,f,j,g) = INTERMOVOV(n,f,j,g) -
     &  (1.0)*T2(m,n,e,f)*V(e,j,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOV(n,f,j,g)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/2*( (T2+{menf}*T2{enfo}) * V{oijm})
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
     &  (1.0/2.0)*T2(m,n,e,f)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(m,o)*V(o,j,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{mejf}*T2{enfo}) * V{niom})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,j,n,o) = INTERMOOOO(m,j,n,o) -
     &  (1.0/4.0)*T2(m,j,e,f)*T2(n,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOOOO(m,j,n,o)*V(n,o,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{eifo}) * V{omjn})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,n,i,o) = INTERMOOOO(m,n,i,o) -
     &  (1.0/4.0)*T2(m,n,e,f)*T2(i,o,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOOOO(m,n,i,o)*V(o,j,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     1/2*( (T2+{menf}*T2{fmgn}) * V{eijg})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVV(e,g)*V(e,j,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/4*( (T2+{mejf}*V{egfh}) * T2{gihm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        DO  h  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,g,h) = INTERMOOVV(m,j,g,h) -
     &  (1.0/4.0)*T2(m,j,e,f)*V(e,f,g,h)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        DO  h  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOOVV(m,j,g,h)*T2(i,m,g,h)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/2*( (F{mf}*T2{eifm}) * T1+{je})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) +
     & (1.0/2.0)*FOCK(m,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,i)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
c     -1/2*( (T1+{me}*T2{eifm}) * F{jf})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) -
     &  (1.0/2.0)*T1(m,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,f)*FOCK(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/2*( (T2+{mejf}*T1{fm}) * F{ei})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) -
     &  (1.0/2.0)*T2(m,j,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,j)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
c     1/2*( (T2+{mejf}*F{em}) * T1{fi})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     & (1.0/2.0)*T2(m,j,e,f)*FOCK(e,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/2*( (T1+{me}*T1{ei}) * F{jm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/2.0)*T1(m,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(m,i)*FOCK(j,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/2*( (T1+{je}*T1{em}) * F{mi})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(j,m) = INTERMOO(j,m) -
     &  (1.0/2.0)*T1(j,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(j,m)*FOCK(m,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1*( (T1+{je}*F{ef}) * T1{fi})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/3, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     & T1(j,e)*FOCK(e,f)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( (T1+{me}*T2{emfn}) * V{nijf})
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
     &  (1.0/2.0)*T1(m,e)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(f,n)*V(n,j,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1*( (T1+{me}*T2{eifn}) * V{nmjf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,i,f,n) = INTERMOOVO(m,i,f,n) +
     & T1(m,e)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOOVO(m,i,f,n)*V(n,j,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/2*( (V{eifg}*T1{gm}) * T2+{mejf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(e,i,f,m) = INTERMVOVO(e,i,f,m) -
     &  (1.0/2.0)*V(e,f,i,g)*T1(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVOVO(e,i,f,m)*T2(m,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (V{efmg}*T2{figm}) * T1+{je})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) +
     & (1.0/2.0)*V(e,m,f,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,i)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/2*( (T1+{me}*V{efjg}) * T2{figm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(m,f,j,g) = INTERMOVOV(m,f,j,g) -
     &  (1.0/2.0)*T1(m,e)*V(e,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOV(m,f,j,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/4*( (V{minf}*T2{emfn}) * T1+{je})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,e) = INTERMOV(i,e) -
     &  (1.0/4.0)*V(m,n,i,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T2+{mejf}*V{emfg}) * T1{gi})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(j,g) = INTERMOV(j,g) +
     & (1.0/2.0)*T2(m,j,e,f)*V(e,f,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1*( (T2+{mejf}*T1{fn}) * V{einm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,j,n) = INTERMOVOO(m,e,j,n) +
     & T2(m,j,e,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOO(m,e,j,n)*V(e,n,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( (T2+{menf}*T1{fn}) * V{eijm})
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
     &  (1.0/2.0)*T2(m,n,e,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(m,e)*V(e,j,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/4*( (T2+{menf}*V{emjn}) * T1{fi})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,j) = INTERMVO(f,j) -
     &  (1.0/4.0)*T2(m,n,e,f)*V(e,j,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(f,j)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1*( (T1+{me}*V{eifm}) * T1+{jf})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) -
     &  (1.0)*T1(m,e)*V(e,f,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1*( (T1+{me}*T1{en}) * V{nijm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(m,n)*V(n,j,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1*( (V{mejf}*T1{fm}) * T1{ei})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) -
     &  (1.0)*V(m,j,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,j)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1*( (T1+{me}*T1{fm}) * V{eijf})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(e,f) = INTERMVV(e,f) -
     &  (1.0)*T1(m,e)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVV(e,f)*V(e,j,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1*( (T1+{me}*V{emjf}) * T1{fi})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     & T1(m,e)*V(e,j,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1*( (V{eimf}*T1{fm}) * T1+{je})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) +
     & V(e,m,i,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,i)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1/6*( ( (T2+{mejf}*T1{em}) * T2{fign}) * F{ng})
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
        INTERMOV(j,f) = INTERMOV(j,f) +
     & (1.0/6.0)*T2(m,j,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(j,i,g,n) = INTERMOOVO(j,i,g,n) +
     &  INTERMOV(j,f)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOOVO(j,i,g,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/3*( ( (F{ng}*T1{gm}) * T2+{mejf}) * T2{eifn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOO(n,m) = INTERMOO(n,m) +
     & (1.0/3.0)*FOCK(n,g)*T1(m,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOV(n,e,j,f) = INTERMOVOV(n,e,j,f) +
     &  INTERMOO(n,m)*T2(m,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOV(n,e,j,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*T2{eigm}) * T2+{nfjg}) * F{fn})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(i,g) = INTERMOV(i,g) -
     &  (1.0/6.0)*T1(m,e)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(i,n,f,j) = INTERMOOVO(i,n,f,j) +
     &  INTERMOV(i,g)*T2(n,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOOVO(i,n,f,j)*FOCK(f,n)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/3*( ( (T1+{me}*F{en}) * T2+{nfjg}) * T2{figm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,n) = INTERMOO(m,n) +
     & (1.0/3.0)*T1(m,e)*FOCK(e,n)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(m,f,j,g) = INTERMOVOV(m,f,j,g) +
     &  INTERMOO(m,n)*T2(n,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOV(m,f,j,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*T2{emgn}) * T2+{nfjg}) * F{fi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(g,n) = INTERMVO(g,n) -
     &  (1.0/6.0)*T1(m,e)*T2(m,n,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO2(NOS+1:NOS+NVS,1:NOS))
        INTERMVO2 = 0.0d0
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO2(f,j) = INTERMVO2(f,j) +
     &  INTERMVO(g,n)*T2(n,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO2(f,j)*FOCK(f,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO2)
        ENDIF
        ENDIF
!
c     1/6*( ( (T2+{mfng}*F{fm}) * T2{eign}) * T1+{je})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(n,g) = INTERMOV(n,g) +
     & (1.0/6.0)*T2(m,n,f,g)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) +
     &  INTERMOV(n,g)*T2(i,n,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,i)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
        ENDIF
!
c     2/3*( ( (T1+{me}*T2{eign}) * F{fm}) * T2+{nfjg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(m,i,g,n) = INTERMOOVO(m,i,g,n) +
     & (2.0/3.0)*T1(m,e)*T2(i,n,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(i,g,n,f) = INTERMOVOV(i,g,n,f) +
     &  INTERMOOVO(m,i,g,n)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOVOV(i,g,n,f)*T2(n,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T2+{mejf}*T1{en}) * F{ng}) * T2{figm})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,j,f,n) = INTERMOOVO(m,j,f,n) -
     &  (2.0/3.0)*T2(m,j,e,f)*T1(n,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVV(m,j,f,g) = INTERMOOVV(m,j,f,g) +
     &  INTERMOOVO(m,j,f,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOOVV(m,j,f,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{ng}*T2{fmgn}) * T2+{mejf}) * T1{ei})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(f,m) = INTERMVO(f,m) -
     &  (1.0/6.0)*FOCK(n,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO2(NOS+1:NOS+NVS,1:NOS))
        INTERMVO2 = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO2(e,j) = INTERMVO2(e,j) +
     &  INTERMVO(f,m)*T2(m,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO2(e,j)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO2)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{menf}*T1{em}) * T2{fign}) * F{jg})
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
        INTERMOV(n,f) = INTERMOV(n,f) -
     &  (1.0/6.0)*T2(m,n,e,f)*T1(m,e)
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV2(i,g)*FOCK(j,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV2)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{nfjg}*T2{fmgn}) * T1+{me}) * F{ei})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOO(j,m) = INTERMOO(j,m) +
     & (1.0/12.0)*T2(n,j,f,g)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(j,e) = INTERMOV(j,e) +
     &  INTERMOO(j,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T1+{je}*T2{emgn}) * T2+{mfng}) * F{fi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(j,m,g,n) = INTERMOOVO(j,m,g,n) -
     &  (1.0/12.0)*T1(j,e)*T2(m,n,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     &  INTERMOOVO(j,m,g,n)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,f)*FOCK(f,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/4*( ( (T2+{mfng}*T2{fign}) * F{em}) * T1+{je})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/4.0)*T2(m,n,f,g)*T2(i,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(i,e) = INTERMOV(i,e) +
     &  INTERMOO(m,i)*FOCK(e,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{menf}*T1{ei}) * T2{fmgn}) * F{jg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,n,f,i) = INTERMOOVO(m,n,f,i) +
     & (1.0/12.0)*T2(m,n,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(i,g) = INTERMOV(i,g) +
     &  INTERMOOVO(m,n,f,i)*T2(m,n,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,g)*FOCK(j,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/4*( ( (T2+{mejf}*T2{emfn}) * F{ng}) * T1{gi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(j,n) = INTERMOO(j,n) -
     &  (1.0/4.0)*T2(m,j,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(j,g) = INTERMOV(j,g) +
     &  INTERMOO(j,n)*FOCK(n,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{menf}*T2{eifn}) * T1{gm}) * F{jg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/12.0)*T2(m,n,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(i,g) = INTERMOV(i,g) +
     &  INTERMOO(m,i)*T1(m,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,g)*FOCK(j,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     1/2*( ( (T1+{me}*F{nm}) * T2{eifn}) * T1+{jf})
c      Commutator = 3, PT Order = 5
!

        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(e,n) = INTERMVO(e,n) +
     & (1.0/2.0)*T1(m,e)*FOCK(n,m)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) +
     &  INTERMVO(e,n)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( ( (T1+{me}*F{eg}) * T2{figm}) * T1+{jf})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(m,g) = INTERMOV(m,g) +
     & (1.0/2.0)*T1(m,e)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(f,i) = INTERMVO(f,i) +
     &  INTERMOV(m,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(f,i)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/3*( ( (T1+{me}*T2{emfn}) * T1+{jf}) * F{ni})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(f,n) = INTERMVO(f,n) -
     &  (1.0/3.0)*T1(m,e)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(n,j) = INTERMOO(n,j) +
     &  INTERMVO(f,n)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(n,j)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/6*( ( (T1+{nf}*T2{eifn}) * T1+{me}) * F{jm})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) -
     &  (1.0/6.0)*T1(n,f)*T2(i,n,e,f)
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(i,m)*FOCK(j,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/2*( ( (T1+{me}*T2{eigm}) * F{fg}) * T1+{jf})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(i,g) = INTERMOV(i,g) -
     &  (1.0/2.0)*T1(m,e)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV2(1:NOS,NOS+1:NOS+NVS))
        INTERMOV2 = 0.0d0
        DO  i  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV2(i,f) = INTERMOV2(i,f) +
     &  INTERMOV(i,g)*FOCK(f,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV2(i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV2)
        ENDIF
!
c     1/2*( ( (F{nm}*T1{fn}) * T2+{mejf}) * T1{ei})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,f) = INTERMOV(m,f) +
     & (1.0/2.0)*FOCK(n,m)*T1(n,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) +
     &  INTERMOV(m,f)*T2(m,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMVO(e,j)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/3*( ( (T2+{menf}*T1{fn}) * T1{ei}) * F{jm})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(m,e) = INTERMOV(m,e) -
     &  (1.0/3.0)*T2(m,n,e,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) +
     &  INTERMOV(m,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(m,i)*FOCK(j,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/6*( ( (T2+{mejf}*T1{em}) * T1{fn}) * F{ni})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) -
     &  (1.0/6.0)*T2(m,j,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(j,n) = INTERMOO(j,n) +
     &  INTERMOV(j,f)*T1(n,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  n  = 1, NOS
        HT(j,i) = HT(j,i) +
     &  INTERMOO(j,n)*FOCK(n,i)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/2*( ( (F{eg}*T1{gm}) * T2+{mejf}) * T1{fi})
c      Commutator = 3, PT Order = 5
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
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     &  INTERMVO(e,m)*T2(m,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/2*( ( (T2+{mejf}*T1{fm}) * F{eg}) * T1{gi})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) -
     &  (1.0/2.0)*T2(m,j,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(j,g) = INTERMOV(j,g) +
     &  INTERMVO(e,j)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     -1/6*( ( (T1+{me}*F{ei}) * T1{fm}) * T1+{jf})
c      Commutator = 3, PT Order = 6
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/6.0)*T1(m,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) +
     &  INTERMOO(m,i)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/2*( ( (T1+{me}*T1{ei}) * F{fm}) * T1+{jf})
c      Commutator = 3, PT Order = 6
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/2.0)*T1(m,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) +
     &  INTERMOO(m,i)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*T1{ei}) * T1{fm}) * F{jf})
c      Commutator = 3, PT Order = 6
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,i) = INTERMOO(m,i) -
     &  (1.0/6.0)*T1(m,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) +
     &  INTERMOO(m,i)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(i,f)*FOCK(j,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF
!
c     -1/2*( ( (T1+{je}*T1{em}) * F{mf}) * T1{fi})
c      Commutator = 3, PT Order = 6
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(j,m) = INTERMOO(j,m) -
     &  (1.0/2.0)*T1(j,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     &  INTERMOO(j,m)*FOCK(m,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        HT(j,i) = HT(j,i) +
     &  INTERMOV(j,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
        ENDIF

      END
