      SUBROUTINE PPPP(NOS,NVS,NVAS,NACTV,T1,T2,
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
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVO(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMOVVO(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: INTERMVVVV(:,:,:,:)

!
c     1*V{acbd}
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        VT(a,b,c,d) = VT(a,b,c,d) + V(a,b,c,d)
        ENDDO; ENDDO; ENDDO; ENDDO
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*V{mcnd}*T2{ambn}
c      Commutator = 1, PT Order = 2
!
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      & (1.0/2.0)*V(m,n,c,d)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*T2+{mcnd}*V{ambn}
c      Commutator = 1, PT Order = 2
!
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      & (1.0/2.0)*T2(m,n,c,d)*V(a,b,m,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*V{acmd}*T1{bm}
c      Commutator = 1, PT Order = 3
!
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      & -(1.0)*V(a,m,c,d)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*V{bcmd}*T1{am}
c      Commutator = 1, PT Order = 3
!
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      & V(b,m,c,d)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*T1+{mc}*V{ambd}
c      Commutator = 1, PT Order = 3
!
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      & -(1.0)*T1(m,c)*V(a,b,m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*T1+{md}*V{ambc}
c      Commutator = 1, PT Order = 3
!
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      & T1(m,d)*V(a,b,m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (T2+{mdnc}*F{om}) * T2{anbo})
c      Commutator = 2, PT Order = 2
!
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOVO(d,n,c,o) = INTERMVOVO(d,n,c,o) +
!      & -(1.0)*T2(m,n,d,c)*FOCK(o,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVO(d,n,c,o)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( (F{ae}*T2{bmen}) * T2+{mdnc})
c      Commutator = 2, PT Order = 2
!
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,m,n) = INTERMVVOO(a,b,m,n) +
!      & (1.0/4.0)*FOCK(a,e)*T2(m,n,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(a,b,m,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( (F{be}*T2{amen}) * T2+{mdnc})
c      Commutator = 2, PT Order = 2
!
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,m,n) = INTERMVVOO(b,a,m,n) +
!      & -(1.0/4.0)*FOCK(b,e)*T2(m,n,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(b,a,m,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( (T2+{menc}*F{ed}) * T2{ambn})
c      Commutator = 2, PT Order = 2
!
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,c,d) = INTERMOOVV(m,n,c,d) +
!      & -(1.0/4.0)*T2(m,n,e,c)*FOCK(e,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOOVV(m,n,c,d)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( (T2+{mend}*F{ec}) * T2{ambn})
c      Commutator = 2, PT Order = 2
!
!         IF(ANY((/1, 3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,d,c) = INTERMOOVV(m,n,d,c) +
!      & (1.0/4.0)*T2(m,n,e,d)*FOCK(e,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOOVV(m,n,d,c)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T2+{menc}*V{emod}) * T2{anbo})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVOV(n,c,o,d) = INTERMOVOV(n,c,o,d) +
!      & T2(m,n,e,c)*V(e,o,m,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVOV(n,c,o,d)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (V{bmoe}*T2{aneo}) * T2+{mdnc})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVO(b,m,a,n) = INTERMVOVO(b,m,a,n) +
!      & V(b,o,m,e)*T2(n,o,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVO(b,m,a,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( (T2+{mdnc}*V{ompn}) * T2{apbo})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  p  = 1, NOS
!         INTERMVVOO(d,c,o,p) = INTERMVVOO(d,c,o,p) +
!      & (1.0/4.0)*T2(m,n,d,c)*V(o,p,m,n)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  p  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(d,c,o,p)*T2(p,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (T2+{mend}*V{emoc}) * T2{anbo})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVOV(n,d,o,c) = INTERMOVOV(n,d,o,c) +
!      & -(1.0)*T2(m,n,e,d)*V(e,o,m,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVOV(n,d,o,c)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/8*( (T2+{menf}*V{ecfd}) * T2{ambn})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         DO  f  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,c,d) = INTERMOOVV(m,n,c,d) +
!      & (1.0/8.0)*T2(m,n,e,f)*V(e,f,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOOVV(m,n,c,d)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOOVV)
!         ENDIF
!
c     -1/4*( (T2+{menf}*T2{amfn}) * V{bdec})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVV(e,a)*V(b,e,d,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     1/4*( (T2+{menf}*T2{bmfn}) * V{adec})
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
        INTERMVV(e,b) = INTERMVV(e,b) +
     & (1.0/4.0)*T2(m,n,e,f)*T2(m,n,b,f)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVV(e,b)*V(a,e,d,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
c     -1/4*( (T2+{mend}*T2{emfn}) * V{acbf})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVV = 0.0d0
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVV(d,f) = INTERMVV(d,f) -
     &          (1.0/4.0)*T2(m,n,e,d)*T2(m,n,e,f)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVV(d,f)*V(a,b,c,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (V{amoe}*T2{bneo}) * T2+{mdnc})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVO(a,m,b,n) = INTERMVOVO(a,m,b,n) +
!      & -(1.0)*V(a,o,m,e)*T2(n,o,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVO(a,m,b,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!
c     1/4*( (T2+{menc}*T2{emfn}) * V{adbf})
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
        INTERMVV(c,f) = INTERMVV(c,f) +
     & (1.0/4.0)*T2(m,n,e,c)*T2(m,n,e,f)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVV(c,f)*V(a,b,d,f)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/8*( (V{aebf}*T2{fmen}) * T2+{mdnc})
c      Commutator = 2, PT Order = 3
!
!         IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
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
!         INTERMVVOO(a,b,m,n) = INTERMVVOO(a,b,m,n) +
!      & (1.0/8.0)*V(a,b,e,f)*T2(m,n,f,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(a,b,m,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!
c     1/2*( (T2+{menc}*T2{bmfn}) * V{afed})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,NOS+1:NOS+NVS))
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVVVV(e,c,b,f)*V(a,e,f,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/2*( (T2+{menc}*T2{amfn}) * V{bfed})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,NOS+1:NOS+NVS))
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVVVV(e,c,a,f)*V(b,e,f,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     1/2*( (T2+{mend}*T2{amfn}) * V{bfec})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(e,d,a,f) = INTERMVVVV(e,d,a,f) +
     & (1.0/2.0)*T2(m,n,e,d)*T2(m,n,a,f)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVVVV(e,d,a,f)*V(b,e,f,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1/2*( (T2+{mend}*T2{bmfn}) * V{afec})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMVVVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,
     &  NOS+1:NOS+NVS,NOS+1:NOS+NVS))
        INTERMVVVV = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        INTERMVVVV(e,d,b,f) = INTERMVVVV(e,d,b,f) -
     &          (1.0/2.0)*T2(m,n,e,d)*T2(m,n,b,f)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  e  = NOS+1, NOS+NVS
        DO  f  = NOS+1, NOS+NVS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVVVV(e,d,b,f)*V(a,e,f,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVVVV)
        ENDIF
!
c     -1*( (T2+{mend}*T2{bneo}) * V{amoc})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  b  = NOS+1, NOS+NVS
        IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
        IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,d,b,o) = INTERMOVVO(m,d,b,o) -
     &          (1.0)*T2(m,n,e,d)*T2(n,o,b,e)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMOVVO(m,d,b,o)*V(a,o,m,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1*( (T2+{mend}*T2{aneo}) * V{bmoc})
c      Commutator = 2, PT Order = 3
!
        IF(ANY((/4, 5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOVVO(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
        INTERMOVVO = 0.0d0
        DO  a  = NOS+1, NOS+NVS
        IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
        IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  o  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOVVO(m,d,a,o) = INTERMOVVO(m,d,a,o) +
     & T2(m,n,e,d)*T2(n,o,a,e)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMOVVO(m,d,a,o)*V(b,o,m,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     1*( (T2+{menc}*T2{bneo}) * V{amod})
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMOVVO(m,c,b,o)*V(a,o,m,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
c     -1*( (T2+{menc}*T2{aneo}) * V{bmod})
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  o  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMOVVO(m,c,a,o)*V(b,o,m,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOVVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T2+{mdnc}*F{bm}) * T1{an})
c      Commutator = 2, PT Order = 3
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(d,n,c,b) = INTERMVOVV(d,n,c,b) +
!      & (1.0/2.0)*T2(m,n,d,c)*FOCK(b,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(d,n,c,b)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T2+{mdnc}*F{am}) * T1{bn})
c      Commutator = 2, PT Order = 3
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(d,n,c,a) = INTERMVOVV(d,n,c,a) -
!     & (1.0/2.0)*T2(m,n,d,c)*FOCK(a,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(d,n,c,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{md}*T2{ambn}) * F{nc})
c      Commutator = 2, PT Order = 3
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVVO(d,a,b,n) = INTERMVVVO(d,a,b,n) -
!     & (1.0/2.0)*T1(m,d)*T2(m,n,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(d,a,b,n)*FOCK(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T1+{mc}*T2{ambn}) * F{nd})
c      Commutator = 2, PT Order = 3
!
!         IF(NONHF)THEN
!         IF(ANY((/3, 4, 5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
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
!         INTERMVVVO(c,a,b,n) = INTERMVVVO(c,a,b,n) +
!      & (1.0/2.0)*T1(m,c)*T2(m,n,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(c,a,b,n)*FOCK(n,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{me}*V{ecnd}) * T2{ambn})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVOV(m,c,n,d) = INTERMOVOV(m,c,n,d) +
!      & T1(m,e)*V(e,n,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVOV(m,c,n,d)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T2+{mend}*T1{bn}) * V{amec})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,d,b) = INTERMOVVV(m,e,d,b) +
!      & (1.0/2.0)*T2(m,n,e,d)*T1(n,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVVV(m,e,d,b)*V(a,e,m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
c     1/2*( (T2+{menc}*T1{en}) * V{ambd})
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMOV(m,c)*V(a,b,m,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
c     1/2*( (T1+{me}*T2{bmen}) * V{acnd})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVO(b,n)*V(a,n,c,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (V{ambe}*T1{en}) * T2+{mcnd})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVO(a,m,b,n) = INTERMVOVO(a,m,b,n) +
!      & V(a,b,m,e)*T1(n,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVO(a,m,b,n)*T2(m,n,c,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T2+{menc}*T1{an}) * V{bmed})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVVV(m,e,c,a)*V(b,e,m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
c     -1/2*( (T1+{me}*T2{amen}) * V{bcnd})
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
        DO  c  = NOS+1, NOS+NVS
        IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
        IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  n  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMVO(a,n)*V(b,n,c,d)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMVO)
        ENDIF
!
c     -1/2*( (T2+{mend}*T1{en}) * V{ambc})
c      Commutator = 2, PT Order = 4
!
        IF(ANY((/5, 6/)==DUCC_MODEL))THEN
        ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
        INTERMOV = 0.0d0
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        DO  n  = 1, NOS
        DO  e  = NOS+1, NOS+NVS
        INTERMOV(m,d) = INTERMOV(m,d) -
     &          (1.0/2.0)*T2(m,n,e,d)*T1(n,e)
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
        DO  d  = NOS+1, NOS+NVS
        IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
        IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
        DO  m  = 1, NOS
        VT(a,b,c,d) = VT(a,b,c,d) +
     &  INTERMOV(m,d)*V(a,b,m,c)
        ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
        DEALLOCATE(INTERMOV)
        ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T2+{mend}*T1{an}) * V{bmec})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOVVV(m,e,d,a) = INTERMOVVV(m,e,d,a) -
!     & (1.0/2.0)*T2(m,n,e,d)*T1(n,a)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVVV(m,e,d,a)*V(b,e,m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T2+{menc}*T1{bn}) * V{amed})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVVV(m,e,c,b)*V(a,e,m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T1+{mc}*T2{bnem}) * V{aend})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
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
!         INTERMVVOV(c,b,n,e) = INTERMVVOV(c,b,n,e) +
!      & (1.0/2.0)*T1(m,c)*T2(n,m,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,b,n,e)*V(a,n,e,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{mc}*T2{anem}) * V{bend})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
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
!         INTERMVVOV(c,a,n,e) = INTERMVVOV(c,a,n,e) -
!     & (1.0/2.0)*T1(m,c)*T2(n,m,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,a,n,e)*V(b,n,e,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{md}*V{nmoc}) * T2{aobn})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOOV(d,n,o,c) = INTERMVOOV(d,n,o,c) -
!     & (1.0/2.0)*T1(m,d)*V(n,o,m,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOOV(d,n,o,c)*T2(o,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T1+{mc}*V{nmod}) * T2{aobn})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOOV(c,n,o,d) = INTERMVOOV(c,n,o,d) +
!      & (1.0/2.0)*T1(m,c)*V(n,o,m,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOOV(c,n,o,d)*T2(o,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{md}*T2{bnem}) * V{aenc})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(d,b,n,e) = INTERMVVOV(d,b,n,e) -
!     & (1.0/2.0)*T1(m,d)*T2(n,m,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,b,n,e)*V(a,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T1+{md}*T2{anem}) * V{benc})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(d,a,n,e) = INTERMVVOV(d,a,n,e) +
!      & (1.0/2.0)*T1(m,d)*T2(n,m,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,a,n,e)*V(b,n,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (V{bmon}*T1{ao}) * T2+{mdnc})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOOV(b,m,n,a) = INTERMVOOV(b,m,n,a) -
!     & (1.0/2.0)*V(b,o,m,n)*T1(o,a)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOOV(b,m,n,a)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (V{amon}*T1{bo}) * T2+{mdnc})
c      Commutator = 2, PT Order = 4
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOOV(a,m,n,b) = INTERMVOOV(a,m,n,b) +
!      & (1.0/2.0)*V(a,o,m,n)*T1(o,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOOV(a,m,n,b)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T1+{me}*T1{bm}) * V{adec})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/2.0)*T1(m,e)*T1(m,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVV(e,b)*V(a,e,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{mc}*V{ambn}) * T1+{nd})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
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
!         INTERMVVVO(c,a,b,n) = INTERMVVVO(c,a,b,n) +
!      & T1(m,c)*V(a,b,m,n)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(c,a,b,n)*T1(n,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (V{mcnd}*T1{am}) * T1{bn})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(c,n,d,a) = INTERMVOVV(c,n,d,a) +
!      & V(m,n,c,d)*T1(m,a)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(c,n,d,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{me}*T1{am}) * V{bdec})
c      Commutator = 2, PT Order = 5
!
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
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVV(e,a)*V(b,e,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( (T1+{md}*T1{em}) * V{acbe})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(d,e) = INTERMVV(d,e) -
!     & (1.0/2.0)*T1(m,d)*T1(m,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVV(d,e)*V(a,b,c,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( (T1+{mc}*T1{em}) * V{adbe})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
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
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  e  = NOS+1, NOS+NVS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVV(c,e)*V(a,b,d,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{mc}*V{bdnm}) * T1{an})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVVO(c,b,d,n) = INTERMVVVO(c,b,d,n) +
!      & T1(m,c)*V(b,n,d,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(c,b,d,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (T1+{mc}*V{adnm}) * T1{bn})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVVO(c,a,d,n) = INTERMVVVO(c,a,d,n) -
!     & (1.0)*T1(m,c)*V(a,n,d,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(c,a,d,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1*( (T1+{md}*V{bcnm}) * T1{an})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVVO(d,b,c,n) = INTERMVVVO(d,b,c,n) -
!     & (1.0)*T1(m,d)*V(b,n,c,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(d,b,c,n)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1*( (T1+{md}*V{acnm}) * T1{bn})
c      Commutator = 2, PT Order = 5
!
!         IF(ANY((/5, 6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVVO(d,a,c,n) = INTERMVVVO(d,a,c,n) +
!      & T1(m,d)*V(a,n,c,m)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(d,a,c,n)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (F{oe}*T2{bneo}) * T2+{mdnc}) * T1{am})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(b,n) = INTERMVO(b,n) -
!     & (1.0/6.0)*FOCK(o,e)*T2(n,o,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(b,m,d,c) = INTERMVOVV(b,m,d,c) +
!      &  INTERMVO(b,n)*T2(m,n,d,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(b,m,d,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -2/3*( ( (F{oe}*T1{em}) * T2+{mdnc}) * T2{anbo})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(o,m) = INTERMOO(o,m) -
!     & (2.0/3.0)*FOCK(o,e)*T1(m,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVOV(o,d,n,c) = INTERMOVOV(o,d,n,c) +
!      &  INTERMOO(o,m)*T2(m,n,d,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVOV(o,d,n,c)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (F{oe}*T2{aneo}) * T2+{mdnc}) * T1{bm})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,n) = INTERMVO(a,n) +
!      & (1.0/6.0)*FOCK(o,e)*T2(n,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(a,m,d,c) = INTERMVOVV(a,m,d,c) +
!      &  INTERMVO(a,n)*T2(m,n,d,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(a,m,d,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T2+{mend}*T1{em}) * T2{aobn}) * F{oc})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,d) = INTERMOV(n,d) +
!      & (1.0/6.0)*T2(m,n,e,d)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVVOV(d,a,o,b) = INTERMVVOV(d,a,o,b) +
!      &  INTERMOV(n,d)*T2(o,n,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,a,o,b)*FOCK(o,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( ( (F{oe}*T1{ao}) * T2{bmen}) * T2+{mdnc})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/4.0)*FOCK(o,e)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,m,n) = INTERMVVOO(a,b,m,n) +
!      &  INTERMVV(e,a)*T2(m,n,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(a,b,m,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T1+{me}*T2{bmeo}) * T2+{ncod}) * F{an})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(b,o) = INTERMVO(b,o) -
!     & (1.0/6.0)*T1(m,e)*T2(m,o,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOVV(b,n,c,d) = INTERMVOVV(b,n,c,d) +
!      &  INTERMVO(b,o)*T2(n,o,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(b,n,c,d)*FOCK(a,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/12*( ( (F{oc}*T1{eo}) * T2+{mend}) * T2{ambn})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(c,e) = INTERMVV(c,e) -
!     & (1.0/12.0)*FOCK(o,c)*T1(o,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(c,m,n,d) = INTERMVOOV(c,m,n,d) +
!      &  INTERMVV(c,e)*T2(m,n,e,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOOV(c,m,n,d)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -2/3*( ( (T1+{me}*F{en}) * T2+{ncod}) * T2{ambo})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOO(1:NOS,1:NOS))
!         INTERMOO = 0.0d0
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOO(m,n) = INTERMOO(m,n) -
!     & (2.0/3.0)*T1(m,e)*FOCK(e,n)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVOV(1:NOS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMOVOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOVOV(m,c,o,d) = INTERMOVOV(m,c,o,d) +
!      &  INTERMOO(m,n)*T2(n,o,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVOV(m,c,o,d)*T2(m,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T2+{menc}*T1{em}) * T2{aobn}) * F{od})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) -
!     & (1.0/6.0)*T2(m,n,e,c)*T1(m,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
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
!         INTERMVVOV(c,a,o,b) = INTERMVVOV(c,a,o,b) +
!      &  INTERMOV(n,c)*T2(o,n,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,a,o,b)*FOCK(o,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/12*( ( (T1+{me}*F{am}) * T2{bneo}) * T2+{ndoc})
c      Commutator = 3, PT Order = 4
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
!         INTERMVV(e,a) = INTERMVV(e,a) -
!     & (1.0/12.0)*T1(m,e)*FOCK(a,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,n,o) = INTERMVVOO(a,b,n,o) +
!      &  INTERMVV(e,a)*T2(n,o,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(a,b,n,o)*T2(n,o,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{me}*T2{ameo}) * T2+{ncod}) * F{bn})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,o) = INTERMVO(a,o) +
!      & (1.0/6.0)*T1(m,e)*T2(m,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVOVV(a,n,c,d) = INTERMVOVV(a,n,c,d) +
!      &  INTERMVO(a,o)*T2(n,o,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(a,n,c,d)*FOCK(b,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/12*( ( (F{od}*T1{eo}) * T2+{menc}) * T2{ambn})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(d,e) = INTERMVV(d,e) +
!      & (1.0/12.0)*FOCK(o,d)*T1(o,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOOV(NOS+1:NOS+NVS,1:NOS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVOOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOOV(d,m,n,c) = INTERMVOOV(d,m,n,c) +
!      &  INTERMVV(d,e)*T2(m,n,e,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOOV(d,m,n,c)*T2(m,n,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T2+{neoc}*F{en}) * T2{ambo}) * T1+{md})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(o,c) = INTERMOV(o,c) -
!     & (1.0/6.0)*T2(n,o,e,c)*FOCK(e,n)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
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
!         DO  o  = 1, NOS
!         INTERMVVOV(c,a,m,b) = INTERMVVOV(c,a,m,b) +
!      &  INTERMOV(o,c)*T2(m,o,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,a,m,b)*T1(m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( ( (F{oe}*T1{bo}) * T2{amen}) * T2+{mdnc})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/4.0)*FOCK(o,e)*T1(o,b)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,m,n) = INTERMVVOO(b,a,m,n) +
!      &  INTERMVV(e,b)*T2(m,n,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(b,a,m,n)*T2(m,n,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/12*( ( (T1+{me}*F{bm}) * T2{aneo}) * T2+{ndoc})
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
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(e,b) = INTERMVV(e,b) +
!      & (1.0/12.0)*T1(m,e)*FOCK(b,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,n,o) = INTERMVVOO(b,a,n,o) +
!      &  INTERMVV(e,b)*T2(n,o,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOO(b,a,n,o)*T2(n,o,d,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T2+{neod}*F{en}) * T2{ambo}) * T1+{mc})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(o,d) = INTERMOV(o,d) +
!      & (1.0/6.0)*T2(n,o,e,d)*FOCK(e,n)
!         ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVVOV(d,a,m,b) = INTERMVVOV(d,a,m,b) +
!      &  INTERMOV(o,d)*T2(m,o,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,a,m,b)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/4*( ( (T1+{md}*F{em}) * T2+{ncoe}) * T2{anbo})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVV(NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVV = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVV(d,e) = INTERMVV(d,e) -
!     & (1.0/4.0)*T1(m,d)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVO(d,n,c,o) = INTERMVOVO(d,n,c,o) +
!      &  INTERMVV(d,e)*T2(n,o,c,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVO(d,n,c,o)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/4*( ( (T1+{mc}*F{em}) * T2+{ndoe}) * T2{anbo})
c      Commutator = 3, PT Order = 4
!
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
!      & (1.0/4.0)*T1(m,c)*FOCK(e,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVO(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVOVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVO(c,n,d,o) = INTERMVOVO(c,n,d,o) +
!      &  INTERMVV(c,e)*T2(n,o,d,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVO(c,n,d,o)*T2(n,o,a,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T2+{mend}*T1{am}) * T2{bneo}) * F{oc})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,d,a) = INTERMVOVV(e,n,d,a) -
!     & (1.0/6.0)*T2(m,n,e,d)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(d,a,b,o) = INTERMVVVO(d,a,b,o) +
!      &  INTERMVOVV(e,n,d,a)*T2(n,o,b,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(d,a,b,o)*FOCK(o,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T2+{menc}*T1{am}) * T2{bneo}) * F{od})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
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
!      & (1.0/6.0)*T2(m,n,e,c)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(c,a,b,o)*FOCK(o,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T1+{md}*T2{bmeo}) * T2+{neoc}) * F{an})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(d,b,e,o) = INTERMVVVO(d,b,e,o) -
!     & (1.0/6.0)*T1(m,d)*T2(m,o,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(d,b,n,c) = INTERMVVOV(d,b,n,c) +
!      &  INTERMVVVO(d,b,e,o)*T2(n,o,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,b,n,c)*FOCK(a,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{mc}*T2{bmeo}) * T2+{neod}) * F{an})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
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
!         INTERMVVVO(c,b,e,o) = INTERMVVVO(c,b,e,o) +
!      & (1.0/6.0)*T1(m,c)*T2(m,o,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(c,b,n,d) = INTERMVVOV(c,b,n,d) +
!      &  INTERMVVVO(c,b,e,o)*T2(n,o,e,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,b,n,d)*FOCK(a,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T2+{menc}*T1{bm}) * T2{aneo}) * F{od})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
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
!     & (1.0/6.0)*T2(m,n,e,c)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(c,b,a,o)*FOCK(o,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T2+{mend}*T1{bm}) * T2{aneo}) * F{oc})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVOVV(e,n,d,b) = INTERMVOVV(e,n,d,b) +
!      & (1.0/6.0)*T2(m,n,e,d)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(d,b,a,o) = INTERMVVVO(d,b,a,o) +
!      &  INTERMVOVV(e,n,d,b)*T2(n,o,a,e)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  o  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVVO(d,b,a,o)*FOCK(o,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{md}*T2{ameo}) * T2+{neoc}) * F{bn})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
!         INTERMVVVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVVO(d,a,e,o) = INTERMVVVO(d,a,e,o) +
!      & (1.0/6.0)*T1(m,d)*T2(m,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))

!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(d,a,n,c) = INTERMVVOV(d,a,n,c) +
!      &  INTERMVVVO(d,a,e,o)*T2(n,o,e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,a,n,c)*FOCK(b,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T1+{mc}*T2{ameo}) * T2+{neod}) * F{bn})
c      Commutator = 3, PT Order = 4
!
!         IF(NONHF)THEN
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVVO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS))
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
!         INTERMVVVO(c,a,e,o) = INTERMVVVO(c,a,e,o) -
!     & (1.0/6.0)*T1(m,c)*T2(m,o,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOV(c,a,n,d) = INTERMVVOV(c,a,n,d) +
!      &  INTERMVVVO(c,a,e,o)*T2(n,o,e,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVVO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,a,n,d)*FOCK(b,n)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (F{be}*T1{en}) * T2+{mcnd}) * T1{am})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(b,n) = INTERMVO(b,n) +
!      & (1.0/6.0)*FOCK(b,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(b,m,c,d) = INTERMVOVV(b,m,c,d) +
!      &  INTERMVO(b,n)*T2(m,n,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(b,m,c,d)*T1(m,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (F{ae}*T1{en}) * T2+{mcnd}) * T1{bm})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVO(a,n) = INTERMVO(a,n) -
!     & (1.0/6.0)*FOCK(a,e)*T1(n,e)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVOVV(NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMVOVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVOVV(a,m,c,d) = INTERMVOVV(a,m,c,d) +
!      &  INTERMVO(a,n)*T2(m,n,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVOVV(a,m,c,d)*T1(m,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVOVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/3*( ( (F{ae}*T2{bmen}) * T1+{md}) * T1+{nc})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(a,b,m,n) = INTERMVVOO(a,b,m,n) +
!      & (1.0/3.0)*FOCK(a,e)*T2(m,n,b,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOV(a,b,n,d) = INTERMVVOV(a,b,n,d) +
!      &  INTERMVVOO(a,b,m,n)*T1(m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(a,b,n,d)*T1(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( ( (F{om}*T1{bo}) * T2+{mdnc}) * T1{an})
c      Commutator = 3, PT Order = 5
!
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
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOV(b,d,n,c) = INTERMVVOV(b,d,n,c) +
!      &  INTERMOV(m,b)*T2(m,n,d,c)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(b,d,n,c)*T1(n,a)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/3*( ( (F{be}*T2{amen}) * T1+{md}) * T1+{nc})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVVOO(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,1:NOS))
!         INTERMVVOO = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMVVOO(b,a,m,n) = INTERMVVOO(b,a,m,n) -
!     & (1.0/3.0)*FOCK(b,e)*T2(m,n,a,e)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOV(b,a,n,d) = INTERMVVOV(b,a,n,d) +
!      &  INTERMVVOO(b,a,m,n)*T1(m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOO)
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(b,a,n,d)*T1(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( ( (F{om}*T1{ao}) * T2+{mcnd}) * T1{bn})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMOV(m,a) = INTERMOV(m,a) -
!     & (1.0/2.0)*FOCK(o,m)*T1(o,a)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOV(a,c,n,d) = INTERMVVOV(a,c,n,d) +
!      &  INTERMOV(m,a)*T2(m,n,c,d)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(a,c,n,d)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/6*( ( (T1+{ne}*F{ed}) * T2{ambn}) * T1+{mc})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,d) = INTERMOV(n,d) +
!      & (1.0/6.0)*T1(n,e)*FOCK(e,d)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMVVOV(d,a,m,b) = INTERMVVOV(d,a,m,b) +
!      &  INTERMOV(n,d)*T2(m,n,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,a,m,b)*T1(m,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/6*( ( (T1+{ne}*F{ec}) * T2{ambn}) * T1+{md})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOV(1:NOS,NOS+1:NOS+NVS))
!         INTERMOV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOV(n,c) = INTERMOV(n,c) -
!     & (1.0/6.0)*T1(n,e)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
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
!         INTERMVVOV(c,a,m,b) = INTERMVVOV(c,a,m,b) +
!      &  INTERMOV(n,c)*T2(m,n,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,a,m,b)*T1(m,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/2*( ( (T1+{md}*F{om}) * T2{anbo}) * T1+{nc})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVO(d,o) = INTERMVO(d,o) -
!     & (1.0/2.0)*T1(m,d)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  b  = NOS+1, NOS+NVS
!         IF( b .GT.(NOS+NACTV).AND. b .LE.(NOS+NVAS))CYCLE
!         IF( b .GT.(NOS+NVAS+NACTV).AND. b .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVVOV(d,a,n,b) = INTERMVVOV(d,a,n,b) +
!      &  INTERMVO(d,o)*T2(n,o,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(d,a,n,b)*T1(n,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/2*( ( (T1+{mc}*F{om}) * T2{anbo}) * T1+{nd})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMVO(NOS+1:NOS+NVS,1:NOS))
!         INTERMVO = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  o  = 1, NOS
!         INTERMVO(c,o) = INTERMVO(c,o) +
!      & (1.0/2.0)*T1(m,c)*FOCK(o,m)
!         ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMVVOV(NOS+1:NOS+NVS,NOS+1:NOS+NVS,1:NOS,NOS+1:NOS+NVS))
!         INTERMVVOV = 0.0d0
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
!         INTERMVVOV(c,a,n,b) = INTERMVVOV(c,a,n,b) +
!      &  INTERMVO(c,o)*T2(n,o,a,b)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMVVOV(c,a,n,b)*T1(n,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMVVOV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     -1/3*( ( (T2+{menc}*F{ed}) * T1{am}) * T1{bn})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,c,d) = INTERMOOVV(m,n,c,d) -
!     & (1.0/3.0)*T2(m,n,e,c)*FOCK(e,d)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVVV(n,c,d,a) = INTERMOVVV(n,c,d,a) +
!      &  INTERMOOVV(m,n,c,d)*T1(m,a)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVVV(n,c,d,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF
!
!Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! !Zero! 
c     1/3*( ( (T2+{mend}*F{ec}) * T1{am}) * T1{bn})
c      Commutator = 3, PT Order = 5
!
!         IF(ANY((/6/)==DUCC_MODEL))THEN
!         ALLOCATE(INTERMOOVV(1:NOS,1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOOVV = 0.0d0
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         DO  e  = NOS+1, NOS+NVS
!         INTERMOOVV(m,n,d,c) = INTERMOOVV(m,n,d,c) +
!      & (1.0/3.0)*T2(m,n,e,d)*FOCK(e,c)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         ALLOCATE(INTERMOVVV(1:NOS,NOS+1:NOS+NVS,NOS+1:NOS+NVS,NOS+1:NOS+NVS))
!         INTERMOVVV = 0.0d0
!         DO  a  = NOS+1, NOS+NVS
!         IF( a .GT.(NOS+NACTV).AND. a .LE.(NOS+NVAS))CYCLE
!         IF( a .GT.(NOS+NVAS+NACTV).AND. a .LE.(NOS+NVS))CYCLE
!         DO  c  = NOS+1, NOS+NVS
!         IF( c .GT.(NOS+NACTV).AND. c .LE.(NOS+NVAS))CYCLE
!         IF( c .GT.(NOS+NVAS+NACTV).AND. c .LE.(NOS+NVS))CYCLE
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  m  = 1, NOS
!         DO  n  = 1, NOS
!         INTERMOVVV(n,d,c,a) = INTERMOVVV(n,d,c,a) +
!      &  INTERMOOVV(m,n,d,c)*T1(m,a)
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
!         DO  d  = NOS+1, NOS+NVS
!         IF( d .GT.(NOS+NACTV).AND. d .LE.(NOS+NVAS))CYCLE
!         IF( d .GT.(NOS+NVAS+NACTV).AND. d .LE.(NOS+NVS))CYCLE
!         DO  n  = 1, NOS
!         VT(a,b,c,d) = VT(a,b,c,d) +
!      &  INTERMOVVV(n,d,c,a)*T1(n,b)
!         ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
!         DEALLOCATE(INTERMOVVV)
!         ENDIF

      END
