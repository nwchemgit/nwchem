      SUBROUTINE HHHH(NOS,NVS,T1,T2,FOCK,V,DUCC_MODEL,NONHF,VT)
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
      DOUBLE PRECISION :: VT(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)     ! TRANSFORMED 2E- INTEGRALS
C
C     CLUSTER OPERATORS: OCCUPIED IS NUMBERED 1--NOS AND VIRTUAL IS NUMBERED (NOS+1)--(NOS+NVS)
      DOUBLE PRECISION :: T1(NOS,NOS+1:NOS+NVS)                   ! (OCCUPIED,VIRTUAL)
      DOUBLE PRECISION :: T2(NOS,NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS) ! (OCCUPIED,OCCUPIED,VIRTUAL,VIRTUAL)
C
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOOVV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOOV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVOVO(:,:,:,:)

!
c     1*V{kilj}
c      Commutator = 0, PT Order = 1
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) + V(k,l,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*T2+{lekf}*V{fiej}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     & (1.0/2.0)*T2(l,k,e,f)*V(f,e,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/2*V{lekf}*T2{fiej}
c      Commutator = 1, PT Order = 2
!
        IF(ANY((/1, 2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     & (1.0/2.0)*V(l,k,e,f)*T2(i,j,f,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*T1+{ke}*V{eilj}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     & T1(k,e)*V(e,l,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*T1+{le}*V{eikj}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) -
     &          (1.0)*T1(l,e)*V(e,k,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1*V{ljke}*T1{ei}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     & V(l,k,j,e)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     -1*V{like}*T1{ej}
c      Commutator = 1, PT Order = 3
!
        IF(ANY((/2, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) -
     &          (1.0)*V(l,k,i,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
c     1/4*( (T2+{mekf}*T2{eifj}) * F{lm})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,i,j) = INTERMOOOO(m,k,i,j) +
     & (1.0/4.0)*T2(m,k,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,k,i,j)*FOCK(l,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/4*( (T2+{melf}*T2{eifj}) * F{km})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,l,i,j) = INTERMOOOO(m,l,i,j) -
     &          (1.0/4.0)*T2(m,l,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,l,i,j)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/4*( (T2+{lekf}*T2{ejfm}) * F{mi})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,j,m) = INTERMOOOO(l,k,j,m) -
     &          (1.0/4.0)*T2(l,k,e,f)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(l,k,j,m)*FOCK(m,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     1*( (T2+{lekf}*F{eg}) * T2{figj})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVV(l,k,f,g) = INTERMOOVV(l,k,f,g) +
     & T2(l,k,e,f)*FOCK(e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVV(l,k,f,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/4*( (T2+{lekf}*T2{eifm}) * F{mj})
c      Commutator = 2, PT Order = 2
!
        IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,i,m) = INTERMOOOO(l,k,i,m) +
     & (1.0/4.0)*T2(l,k,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(l,k,i,m)*FOCK(m,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     1*( (V{ejmg}*T2{figm}) * T2+{lekf})
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
        INTERMVOVO(e,j,f,i) = INTERMVOVO(e,j,f,i) +
     & V(e,m,j,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOVO(e,j,f,i)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/4*( (T2+{melf}*T2{emfn}) * V{nikj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOO(l,n) = INTERMOO(l,n) +
     & (1.0/4.0)*T2(m,l,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(l,n)*V(n,k,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/4*( (T2+{menf}*T2{eifn}) * V{ljkm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(m,i)*V(l,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/4*( (T2+{menf}*T2{ejfn}) * V{likm})
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
        INTERMOO(m,j) = INTERMOO(m,j) +
     & (1.0/4.0)*T2(m,n,e,f)*T2(j,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(m,j)*V(l,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1*( (T2+{mekf}*V{emlg}) * T2{figj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(k,f,l,g) = INTERMOVOV(k,f,l,g) +
     & T2(m,k,e,f)*V(e,l,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOV(k,f,l,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/4*( (T2+{mekf}*T2{emfn}) * V{nilj})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(k,n)*V(n,l,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1*( (V{eimg}*T2{fjgm}) * T2+{lekf})
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
        INTERMVOVO(e,i,f,j) = INTERMVOVO(e,i,f,j) -
     &          (1.0)*V(e,m,i,g)*T2(j,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOVO(e,i,f,j)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     -1/8*( (T2+{lekf}*T2{emfn}) * V{minj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,m,n) = INTERMOOOO(l,k,m,n) -
     &          (1.0/8.0)*T2(l,k,e,f)*T2(m,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(l,k,m,n)*V(m,n,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/8*( (T2+{menf}*T2{eifj}) * V{lmkn})
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
        INTERMOOOO(m,n,i,j) = INTERMOOOO(m,n,i,j) -
     &          (1.0/8.0)*T2(m,n,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,n,i,j)*V(l,k,m,n)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1*( (T2+{melf}*V{emkg}) * T2{figj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(l,f,k,g) = INTERMOVOV(l,f,k,g) -
     &          (1.0)*T2(m,l,e,f)*V(e,k,m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOV(l,f,k,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     1/4*( (T2+{lekf}*V{egfh}) * T2{higj})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        DO  h  = NOS+1, NOS+NVS
        INTERMOOVV(l,k,g,h) = INTERMOOVV(l,k,g,h) +
     & (1.0/4.0)*T2(l,k,e,f)*V(e,f,g,h)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        DO  h  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVV(l,k,g,h)*T2(i,j,h,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
!
c     1/2*( (T2+{mekf}*T2{ejfn}) * V{nilm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,k,j,n)*V(n,l,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/2*( (T2+{mekf}*T2{eifn}) * V{njlm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,k,i,n)*V(n,l,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     1/2*( (T2+{melf}*T2{eifn}) * V{njkm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,l,i,n) = INTERMOOOO(m,l,i,n) +
     & (1.0/2.0)*T2(m,l,e,f)*T2(i,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,l,i,n)*V(n,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/2*( (T2+{melf}*T2{ejfn}) * V{nikm})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,l,j,n) = INTERMOOOO(m,l,j,n) -
     &          (1.0/2.0)*T2(m,l,e,f)*T2(j,n,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,l,j,n)*V(n,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1*( (T2+{melf}*T2{fjgm}) * V{eikg})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOOV(e,l,j,g) = INTERMVOOV(e,l,j,g) -
     &          (1.0)*T2(m,l,e,f)*T2(j,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOOV(e,l,j,g)*V(e,k,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     1*( (T2+{melf}*T2{figm}) * V{ejkg})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMVOOV = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOOV(e,l,i,g) = INTERMVOOV(e,l,i,g) +
     & T2(m,l,e,f)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOOV(e,l,i,g)*V(e,k,j,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     1*( (T2+{mekf}*T2{fjgm}) * V{eilg})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOOV(e,k,j,g)*V(e,l,i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     -1*( (T2+{mekf}*T2{figm}) * V{ejlg})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOOV(e,k,i,g)*V(e,l,j,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOOV)
        ENDIF
!
c     1/2*( (T2+{lekf}*F{ej}) * T1{fi})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,k,f,j) = INTERMOOVO(l,k,f,j) +
     & (1.0/2.0)*T2(l,k,e,f)*FOCK(e,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,k,f,j)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1/2*( (T2+{lekf}*F{ei}) * T1{fj})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,k,f,i) = INTERMOOVO(l,k,f,i) -
     &          (1.0/2.0)*T2(l,k,e,f)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,k,f,i)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1/2*( (T1+{le}*T2{eifj}) * F{kf})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,i,f,j) = INTERMOOVO(l,i,f,j) -
     &          (1.0/2.0)*T1(l,e)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,i,f,j)*FOCK(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/2*( (T1+{ke}*T2{eifj}) * F{lf})
c      Commutator = 2, PT Order = 3
!
        IF(NONHF)THEN
        IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(k,i,f,j) = INTERMOOVO(k,i,f,j) +
     & (1.0/2.0)*T1(k,e)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,i,f,j)*FOCK(l,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1*( (T1+{me}*T2{eifj}) * V{kmlf})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(m,i,f,j)*V(k,l,m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1*( (T2+{lekf}*T1{fm}) * V{miej})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(l,e,k,m) = INTERMOVOO(l,e,k,m) -
     &          (1.0)*T2(l,k,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(l,e,k,m)*V(m,e,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( (T2+{melf}*T1{fj}) * V{eikm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,l,j) = INTERMOVOO(m,e,l,j) -
     &          (1.0/2.0)*T2(m,l,e,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(m,e,l,j)*V(e,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T2+{mekf}*T1{fj}) * V{eilm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(m,e,k,j)*V(e,l,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1/2*( (T2+{melf}*T1{fm}) * V{eikj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVO(e,l) = INTERMVO(e,l) +
     & (1.0/2.0)*T2(m,l,e,f)*T1(m,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVO(e,l)*V(e,k,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/2*( (T2+{mekf}*T1{fi}) * V{ejlm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(m,e,k,i)*V(e,l,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( (T2+{mekf}*T1{fm}) * V{eilj})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVO(e,k)*V(e,l,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     1/2*( (T2+{melf}*T1{fi}) * V{ejkm})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(m,e,l,i) = INTERMOVOO(m,e,l,i) +
     & (1.0/2.0)*T2(m,l,e,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(m,e,l,i)*V(e,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( (T1+{me}*T2{eifm}) * V{ljkf})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOV(i,f)*V(l,k,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T1+{me}*T2{ejfm}) * V{likf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     & (1.0/2.0)*T1(m,e)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOV(j,f)*V(l,k,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T1+{ke}*T2{ejfm}) * V{milf})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,j,f,m)*V(m,l,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/2*( (T1+{ke}*T2{eifm}) * V{mjlf})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,i,f,m)*V(m,l,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/2*( (T1+{le}*T2{ejfm}) * V{mikf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,j,f,m) = INTERMOOVO(l,j,f,m) -
     &          (1.0/2.0)*T1(l,e)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,j,f,m)*V(m,k,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/2*( (T1+{le}*T2{eifm}) * V{mjkf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,i,f,m) = INTERMOOVO(l,i,f,m) +
     & (1.0/2.0)*T1(l,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,i,f,m)*V(m,k,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/2*( (T1+{ke}*V{eflg}) * T2{figj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(k,f,l,g) = INTERMOVOV(k,f,l,g) +
     & (1.0/2.0)*T1(k,e)*V(e,l,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOV(k,f,l,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/2*( (T1+{le}*V{efkg}) * T2{figj})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(l,f,k,g) = INTERMOVOV(l,f,k,g) -
     &          (1.0/2.0)*T1(l,e)*V(e,k,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOV(l,f,k,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
!
c     -1/2*( (V{eifg}*T1{gj}) * T2+{lekf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(e,i,f,j) = INTERMVOVO(e,i,f,j) -
     &          (1.0/2.0)*V(e,f,i,g)*T1(j,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOVO(e,i,f,j)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1/2*( (V{ejfg}*T1{gi}) * T2+{lekf})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMVOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVOVO(e,j,f,i) = INTERMVOVO(e,j,f,i) +
     & (1.0/2.0)*V(e,f,j,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMVOVO(e,j,f,i)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVOVO)
        ENDIF
!
c     1*( (T1+{le}*V{fiej}) * T1+{kf})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(l,f,i,j) = INTERMOVOO(l,f,i,j) +
     & T1(l,e)*V(f,e,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(l,f,i,j)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     1*( (V{lekf}*T1{fi}) * T1{ej})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOVOO(l,e,k,i) = INTERMOVOO(l,e,k,i) +
     & V(l,k,e,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(l,e,k,i)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
!
c     -1/2*( (T1+{le}*T1{em}) * V{kimj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(l,m) = INTERMOO(l,m) -
     &          (1.0/2.0)*T1(l,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS

        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(l,m)*V(k,m,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/2*( (T1+{ke}*T1{em}) * V{limj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(k,m) = INTERMOO(k,m) +
     & (1.0/2.0)*T1(k,e)*T1(m,e)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(k,m)*V(l,m,i,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1/2*( (T1+{me}*T1{ej}) * V{likm})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOO(1:NOS,1:NOS))
        INTERMOO = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOO(m,j) = INTERMOO(m,j) +
     & (1.0/2.0)*T1(m,e)*T1(j,e)
        ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(m,j)*V(l,k,i,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     -1/2*( (T1+{me}*T1{ei}) * V{ljkm})
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
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOO(m,i)*V(l,k,j,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOO)
        ENDIF
!
c     1*( (T1+{ke}*V{eilf}) * T1{fj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(k,i,l,f) = INTERMOOOV(k,i,l,f) +
     & T1(k,e)*V(e,l,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(k,i,l,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1*( (T1+{ke}*V{ejlf}) * T1{fi})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(k,j,l,f) = INTERMOOOV(k,j,l,f) -
     &          (1.0)*T1(k,e)*V(e,l,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(k,j,l,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1*( (T1+{le}*V{eikf}) * T1{fj})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(l,i,k,f) = INTERMOOOV(l,i,k,f) -
     &          (1.0)*T1(l,e)*V(e,k,i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,i,k,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     1*( (T1+{le}*V{ejkf}) * T1{fi})
c      Commutator = 2, PT Order = 5
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(l,j,k,f) = INTERMOOOV(l,j,k,f) +
     & T1(l,e)*V(e,k,j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,j,k,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1/4*( ( (T2+{lekf}*T2{ejfm}) * F{mg}) * T1{gi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,j,m) = INTERMOOOO(l,k,j,m) -
     &          (1.0/4.0)*T2(l,k,e,f)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(l,k,j,g) = INTERMOOOV(l,k,j,g) +
     &  INTERMOOOO(l,k,j,m)*FOCK(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,k,j,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{lekf}*T2{ejfm}) * F{mg}) * T1{gi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,j,m) = INTERMOOOO(l,k,j,m) -
     &          (1.0/6.0)*T2(l,k,e,f)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(l,k,j,g) = INTERMOOOV(l,k,j,g) +
     &  INTERMOOOO(l,k,j,m)*FOCK(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,k,j,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T2+{melf}*T1{em}) * T2{gifj}) * F{kg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(l,f) = INTERMOV(l,f) +
     & (1.0/6.0)*T2(m,l,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(l,g,i,j) = INTERMOVOO(l,g,i,j) +
     &  INTERMOV(l,f)*T2(i,j,g,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(l,g,i,j)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
        ENDIF
!
c     -1/6*( ( (F{mg}*T2{eigm}) * T2+{lekf}) * T1{fj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) -
     &          (1.0/6.0)*FOCK(m,g)*T2(i,m,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(i,l,k,f) = INTERMOOOV(i,l,k,f) +
     &  INTERMVO(e,i)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(i,l,k,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T2+{lekf}*T1{em}) * F{mg}) * T2{figj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,k,f,m) = INTERMOOVO(l,k,f,m) -
     &          (2.0/3.0)*T2(l,k,e,f)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMOOVV = 0.0d0
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVV(l,k,f,g) = INTERMOOVV(l,k,f,g) +
     &  INTERMOOVO(l,k,f,m)*FOCK(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVV(l,k,f,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVV)
        ENDIF
        ENDIF
!
c     -2/3*( ( (T1+{me}*T2{eifj}) * F{gm}) * T2+{lfkg})
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
        INTERMOOVO(m,i,f,j) = INTERMOOVO(m,i,f,j) -
     &          (2.0/3.0)*T1(m,e)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
        INTERMOVOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOV(i,f,j,g) = INTERMOVOV(i,f,j,g) +
     &  INTERMOOVO(m,i,f,j)*FOCK(g,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOV(i,f,j,g)*T2(l,k,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{me}*T2{ejfm}) * T2+{lfkg}) * F{gi})
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
     &          (1.0/6.0)*T1(m,e)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(j,l,k,g) = INTERMOOOV(j,l,k,g) +
     &  INTERMOV(j,f)*T2(l,k,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(j,l,k,g)*FOCK(g,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{mekf}*T1{em}) * T2{gifj}) * F{lg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
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
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(k,g,i,j) = INTERMOVOO(k,g,i,j) +
     &  INTERMOV(k,f)*T2(i,j,g,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(k,g,i,j)*FOCK(l,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{me}*T2{eifm}) * T2+{lfkg}) * F{gj})
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
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(i,l,k,g) = INTERMOOOV(i,l,k,g) +
     &  INTERMOV(i,f)*T2(l,k,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(i,l,k,g)*FOCK(g,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/4*( ( (T2+{lekf}*T2{eifm}) * F{mg}) * T1{gj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,i,m) = INTERMOOOO(l,k,i,m) +
     & (1.0/4.0)*T2(l,k,e,f)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(l,k,i,g) = INTERMOOOV(l,k,i,g) +
     &  INTERMOOOO(l,k,i,m)*FOCK(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,k,i,g)*T1(j,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T2+{mflg}*F{fm}) * T2{eigj}) * T1+{ke})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(l,g) = INTERMOV(l,g) +
     & (1.0/6.0)*T2(m,l,f,g)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(l,e,i,j) = INTERMOVOO(l,e,i,j) +
     &  INTERMOV(l,g)*T2(i,j,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(l,e,i,j)*T1(k,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{melf}*T2{eifj}) * T1{gm}) * F{kg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,l,i,j) = INTERMOOOO(m,l,i,j) -
     &          (1.0/12.0)*T2(m,l,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(l,i,j,g) = INTERMOOOV(l,i,j,g) +
     &  INTERMOOOO(m,l,i,j)*T1(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,i,j,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/12*( ( (T2+{lfkg}*T2{fjgm}) * T1+{me}) * F{ei})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,j,m) = INTERMOOOO(l,k,j,m) -
     &          (1.0/12.0)*T2(l,k,f,g)*T2(j,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOOV(l,k,j,e) = INTERMOOOV(l,k,j,e) +
     &  INTERMOOOO(l,k,j,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,k,j,e)*FOCK(e,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{mfkg}*F{fm}) * T2{eigj}) * T1+{le})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(k,g) = INTERMOV(k,g) -
     &          (1.0/6.0)*T2(m,k,f,g)*FOCK(f,m)
        ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(k,e,i,j) = INTERMOVOO(k,e,i,j) +
     &  INTERMOV(k,g)*T2(i,j,e,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOVOO(k,e,i,j)*T1(l,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        ENDIF
        ENDIF
!
c     -1/4*( ( (T2+{mflg}*T2{figj}) * F{em}) * T1+{ke})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOO(m,l,i,j) = INTERMOOOO(m,l,i,j) -
     &          (1.0/4.0)*T2(m,l,f,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOOV(l,i,j,e) = INTERMOOOV(l,i,j,e) +
     &  INTERMOOOO(m,l,i,j)*FOCK(e,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,i,j,e)*T1(k,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{mekf}*T2{eifj}) * T1{gm}) * F{lg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,i,j) = INTERMOOOO(m,k,i,j) +
     & (1.0/12.0)*T2(m,k,e,f)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(k,i,j,g) = INTERMOOOV(k,i,j,g) +
     &  INTERMOOOO(m,k,i,j)*T1(m,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(k,i,j,g)*FOCK(l,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/12*( ( (T2+{lfkg}*T2{figm}) * T1+{me}) * F{ej})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOO(l,k,i,m) = INTERMOOOO(l,k,i,m) +
     & (1.0/12.0)*T2(l,k,f,g)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOOV(l,k,i,e) = INTERMOOOV(l,k,i,e) +
     &  INTERMOOOO(l,k,i,m)*T1(m,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,k,i,e)*FOCK(e,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/4*( ( (T2+{mfkg}*T2{figj}) * F{em}) * T1+{le})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,i,j) = INTERMOOOO(m,k,i,j) +
     & (1.0/4.0)*T2(m,k,f,g)*T2(i,j,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOOOV(k,i,j,e) = INTERMOOOV(k,i,j,e) +
     &  INTERMOOOO(m,k,i,j)*FOCK(e,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(k,i,j,e)*T1(l,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{mekf}*T1{ej}) * T2{figm}) * F{lg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
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
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(k,j,i,g) = INTERMOOOV(k,j,i,g) +
     &  INTERMOOVO(m,k,f,j)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(k,j,i,g)*FOCK(l,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T2+{melf}*T1{ej}) * T2{figm}) * F{kg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,l,f,j) = INTERMOOVO(m,l,f,j) +
     & (1.0/6.0)*T2(m,l,e,f)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(l,j,i,g) = INTERMOOOV(l,j,i,g) +
     &  INTERMOOVO(m,l,f,j)*T2(i,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,j,i,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T2+{melf}*T1{ei}) * T2{fjgm}) * F{kg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,l,f,i) = INTERMOOVO(m,l,f,i) -
     &          (1.0/6.0)*T2(m,l,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(l,i,j,g) = INTERMOOOV(l,i,j,g) +
     &  INTERMOOVO(m,l,f,i)*T2(j,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(l,i,j,g)*FOCK(k,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     1/6*( ( (T2+{mekf}*T1{ei}) * T2{fjgm}) * F{lg})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
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
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOOV(k,i,j,g) = INTERMOOOV(k,i,j,g) +
     &  INTERMOOVO(m,k,f,i)*T2(j,m,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  g  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(k,i,j,g)*FOCK(l,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{ke}*T2{gjem}) * T2+{mflg}) * F{fi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(k,g,j,m) = INTERMOVOO(k,g,j,m) -
     &          (1.0/6.0)*T1(k,e)*T2(j,m,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(k,j,f,l) = INTERMOOVO(k,j,f,l) +
     &  INTERMOVOO(k,g,j,m)*T2(m,l,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,j,f,l)*FOCK(f,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{ke}*T2{giem}) * T2+{mflg}) * F{fj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(k,g,i,m) = INTERMOVOO(k,g,i,m) +
     & (1.0/6.0)*T1(k,e)*T2(i,m,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(k,i,f,l) = INTERMOOVO(k,i,f,l) +
     &  INTERMOVOO(k,g,i,m)*T2(m,l,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,i,f,l)*FOCK(f,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/6*( ( (T1+{le}*T2{gjem}) * T2+{mfkg}) * F{fi})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(l,g,j,m) = INTERMOVOO(l,g,j,m) +
     & (1.0/6.0)*T1(l,e)*T2(j,m,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(l,j,f,k) = INTERMOOVO(l,j,f,k) +
     &  INTERMOVOO(l,g,j,m)*T2(m,k,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,j,f,k)*FOCK(f,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     -1/6*( ( (T1+{le}*T2{giem}) * T2+{mfkg}) * F{fj})
c      Commutator = 3, PT Order = 4
!
        IF(NONHF)THEN
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVOO(1:NOS,NOS+1:NOS+NVS,1:NOS,1:NOS))
        INTERMOVOO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOVOO(l,g,i,m) = INTERMOVOO(l,g,i,m) -
     &          (1.0/6.0)*T1(l,e)*T2(i,m,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(l,i,f,k) = INTERMOOVO(l,i,f,k) +
     &  INTERMOVOO(l,g,i,m)*T2(m,k,f,g)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVOO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,i,f,k)*FOCK(f,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
        ENDIF
!
c     1/3*( ( (T2+{mekf}*T1{ei}) * T1{fj}) * F{lm})
c      Commutator = 3, PT Order = 5
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
     & (1.0/3.0)*T2(m,k,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,k,i,j) = INTERMOOOO(m,k,i,j) +
     &  INTERMOOVO(m,k,f,i)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)

        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,k,i,j)*FOCK(l,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/6*( ( (T1+{me}*F{km}) * T2{eifj}) * T1+{lf})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  k  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(e,k) = INTERMVO(e,k) -
     &          (1.0/6.0)*T1(m,e)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(k,i,f,j) = INTERMOOVO(k,i,f,j) +
     &  INTERMVO(e,k)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,i,f,j)*T1(l,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/6*( ( (F{mi}*T1{fm}) * T2+{lekf}) * T1{ej})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  i  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(i,f) = INTERMOV(i,f) -
     &          (1.0/6.0)*FOCK(m,i)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(i,l,e,k) = INTERMOOVO(i,l,e,k) +
     &  INTERMOV(i,f)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(i,l,e,k)*T1(j,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/6*( ( (T1+{me}*F{lm}) * T2{eifj}) * T1+{kf})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMVO(e,l) = INTERMVO(e,l) +
     & (1.0/6.0)*T1(m,e)*FOCK(l,m)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,i,f,j) = INTERMOOVO(l,i,f,j) +
     &  INTERMVO(e,l)*T2(i,j,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,i,f,j)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/3*( ( (T2+{melf}*T1{ei}) * T1{fj}) * F{km})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(m,l,f,i) = INTERMOOVO(m,l,f,i) -
     &          (1.0/3.0)*T2(m,l,e,f)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(m,l,i,j) = INTERMOOOO(m,l,i,j) +
     &  INTERMOOVO(m,l,f,i)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(m,l,i,j)*FOCK(k,m)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     -1/3*( ( (T1+{le}*T2{ejfm}) * T1+{kf}) * F{mi})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,j,f,m) = INTERMOOVO(l,j,f,m) -
     &          (1.0/3.0)*T1(l,e)*T2(j,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,j,m,k) = INTERMOOOO(l,j,m,k) +
     &  INTERMOOVO(l,j,f,m)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(l,j,m,k)*FOCK(m,i)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     1/3*( ( (T1+{le}*T2{eifm}) * T1+{kf}) * F{mj})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(l,i,f,m) = INTERMOOVO(l,i,f,m) +
     & (1.0/3.0)*T1(l,e)*T2(i,m,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOO(1:NOS,1:NOS,1:NOS,1:NOS))
        INTERMOOOO = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOO(l,i,m,k) = INTERMOOOO(l,i,m,k) +
     &  INTERMOOVO(l,i,f,m)*T1(k,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  m  = 1, NOS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOO(l,i,m,k)*FOCK(m,j)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOO)
        ENDIF
!
c     1/6*( ( (F{mj}*T1{fm}) * T2+{lekf}) * T1{ei})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  j  = 1, NOS
        DO  m  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        INTERMOV(j,f) = INTERMOV(j,f) +
     & (1.0/6.0)*FOCK(m,j)*T1(m,f)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOVO(j,l,e,k) = INTERMOOVO(j,l,e,k) +
     &  INTERMOV(j,f)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(j,l,e,k)*T1(i,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     -1/2*( ( (F{eg}*T1{gi}) * T2+{lekf}) * T1{fj})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  i  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,i) = INTERMVO(e,i) -
     &          (1.0/2.0)*FOCK(e,g)*T1(i,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  i  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(i,l,k,f) = INTERMOOOV(i,l,k,f) +
     &  INTERMVO(e,i)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(i,l,k,f)*T1(j,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF
!
c     -1/2*( ( (T1+{lf}*F{fg}) * T2{giej}) * T1+{ke})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(l,g) = INTERMOV(l,g) -
     &          (1.0/2.0)*T1(l,f)*FOCK(f,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(l,i,e,j) = INTERMOOVO(l,i,e,j) +
     &  INTERMOV(l,g)*T2(i,j,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(l,i,e,j)*T1(k,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/2*( ( (T1+{kf}*F{fg}) * T2{giej}) * T1+{le})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  k  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOV(k,g) = INTERMOV(k,g) +
     & (1.0/2.0)*T1(k,f)*FOCK(f,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOVO(1:NOS,1:NOS,NOS+1:NOS+NVS,1:NOS))
        INTERMOOVO = 0.0d0
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMOOVO(k,i,e,j) = INTERMOOVO(k,i,e,j) +
     &  INTERMOV(k,g)*T2(i,j,g,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOVO(k,i,e,j)*T1(l,e)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOVO)
        ENDIF
!
c     1/2*( ( (F{eg}*T1{gj}) * T2+{lekf}) * T1{fi})
c      Commutator = 3, PT Order = 5
!
        IF(ANY((/6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
        INTERMVO = 0.0d0
        DO  j  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  g  = NOS+1, NOS+NVS
        INTERMVO(e,j) = INTERMVO(e,j) +
     & (1.0/2.0)*FOCK(e,g)*T1(j,g)
        ENDDO; ENDDO; ENDDO
        ALLOCATE(INTERMOOOV(1:NOS,1:NOS,1:NOS,NOS+1:NOS+NVS))
        INTERMOOOV = 0.0d0
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMOOOV(j,l,k,f) = INTERMOOOV(j,l,k,f) +
     &  INTERMVO(e,j)*T2(l,k,e,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        DO  i  = 1, NOS
        DO  j  = 1, NOS
        DO  k  = 1, NOS
        DO  l  = 1, NOS
        DO  f  = NOS+1, NOS+NVS
        VT(k,l,i,j) = VT(k,l,i,j) +
     &  INTERMOOOV(j,l,k,f)*T1(i,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOOOV)
        ENDIF

      END
