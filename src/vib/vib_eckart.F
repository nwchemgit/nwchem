      SUBROUTINE vib_eckart( HESS, HESSP, HESST, COORD, VC ,
     ,     w1,w2)
* $Id$
      IMPLICIT none ! REAL*8 (A-H,O-Z)
c
#include "util.fh"
#include "stdio.fh"
c
      double precision zero, one
      PARAMETER( ZERO=0.D0, ONE=1.D0 )
c
      integer nat, nat3, nhess, nhesst
      COMMON /cvib_HESS/ NAT,NAT3,NHESS,NHESST
c
      double precision HESS(NAT3,NAT3)
      double precision HESSP(NAT3,NAT3)
      double precision HESST(NAT3*(NAT3+1)/2)
      double precision COORD(3,NAT)
      double precision VC(NAT3,6)
      double precision w1(nat3,*),w2(nat3,*)
c
      DOUBLE PRECISION UNIVEC(3), TEST(6,6), VNORM, temp, dotval, rnorm
      double precision test_norm
      integer i, j, k, l, m, n, mu, nu, indx, iatom, iaxis, itemp
C****
C**** construct translation unit vectors;  these are stored in the
C**** first three columns of array VC, the rotation vectors will
C**** be stored in the other 3 columns
C****
      CALL Dfill (6*NAT3,0.0d00,VC,1)
      VNORM = ONE/SQRT(dble(NAT))
      DO 10 I=1,3 ! unit translation vector for x-, y-, and z-directions
        DO 20 IATOM=1,NAT
           VC(3*(IATOM-1)+I,I) = VNORM
   20   CONTINUE
   10 CONTINUE
C****
C**** construct rotation unit vectors; store in last 3 columns of VC
C****
      DO 40 IAXIS=1,3         ! loop over rotation axes
        DO 50 J=1,3
          IF (J.EQ.IAXIS) THEN
            UNIVEC(J) = ONE
          ELSE
            UNIVEC(J) = ZERO
          END IF
   50   CONTINUE
        DO 60 IATOM=1,NAT       ! loop over atoms
          CALL vib_CROSS(COORD(1,IATOM),UNIVEC,
     &      VC(3*(IATOM-1)+1,IAXIS+3))
   60   CONTINUE
   40 CONTINUE
C****
C****   Schmidt orthogonalize the constraint vectors
C****
      DO 70 I=2,6     ! orthogonalize vector I to each of vectors J
        DO 80 J=1,I-1
           dotval=ddot(nat3,vc(1,J),1,vc(1,i),1)
          DO 100 K=1,NAT3
             VC(K,I) = VC(K,I)-dotval*VC(K,J)
  100     CONTINUE
   80   CONTINUE
        RNORM = ddot(nat3,vc(1,I),1,vc(1,I),1)
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
          call dscal(nat3,rnorm,vc(1,I),1)
        ELSE ! dependent constraint vector (linear molecule; set to zero)
          DO 125 K=1,NAT3
            VC(K,I) = ZERO
  125     CONTINUE
          ITEMP = I-3
          WRITE(luout,*)
     &        'Dependent rotation vector no.',ITEMP
          WRITE(luout,*)
     &        ' found in ECKART; assuming linear geometry'
        END IF
   70 CONTINUE
C****
C**** remove rotations and translations from the Hessian matrix
C****
!      DO 130 M=1,NAT3
!        DO 140 N=1,M
!          HESSP(M,N) = HESS(M,N)
!          DO 150 MU=1,6
!            DO 160 J=1,NAT3
!              HESSP(M,N) = HESSP(M,N) !- VC(M,MU)*VC(J,MU)*HESS(J,N)
!     &                               !- HESS(M,J)*VC(J,MU)*VC(N,MU)
!              DO 170 NU=1,6
!                DO 180 K=1,NAT3
!                  HESSP(M,N) = HESSP(M,N) +
!     &                VC(M,MU)*VC(J,MU)*HESS(J,K)*VC(K,NU)*VC(N,NU)
!  180           CONTINUE
!  170         CONTINUE
!  160       CONTINUE
!  150     CONTINUE
!          if(M.ne.N) HESSP(N,M) = HESSP(M,N)
!  140   CONTINUE
!  130 CONTINUE
      call dcopy(nat3*nat3,hess,1,hessp,1)

      call dgemm('n','t',nat3,nat3,6,1d0,vc,nat3,vc,nat3,
     .     0d0,w1,nat3)
      call dgemm('n','n',nat3,nat3,nat3,1d0,w1,nat3,hess,nat3,
     .     0d0,w2,nat3)
      call dgemm('n','n',nat3,nat3,nat3,1d0,w2,nat3,w1,nat3,
     .     1d0,hessp,nat3)
      call dgemm('n','n',nat3,nat3,nat3,-1d0,w1,nat3,hess,nat3,
     .     1d0,hessp,nat3)
      call dgemm('n','n',nat3,nat3,nat3,-1d0,hess,nat3,w1,nat3,
     .     1d0,hessp,nat3)
C****
C**** update triangular hessian to agree with HESS
C****
      DO 190 I=1,NAT3
        DO 200 J=1,I
          indx = i*(i-1)/2 + j
          HESST(indx) = HESSP(I,J)
  200   CONTINUE
  190 CONTINUE
C****
C**** place the projected Hessian in array HESS, and the difference
C**** between projected and unprojected in HESSP
C****
      DO 210 I=1,NAT3
        DO 220 J=1,I
          TEMP = HESS(I,J)
          HESS(I,J) = HESSP(I,J)
          HESSP(I,J) = HESS(I,J)-TEMP
          IF (DABS(HESSP(I,J)).LT.1.D-9) HESSP(I,J) = ZERO
          IF(J.NE.I) then
             HESSP(J,I)=HESSP(I,J)
             HESS(J,I)=HESS(I,J)
          endif
  220   CONTINUE
  210 CONTINUE
C****
C**** output projected Hessian
C****
      if (util_print('eckart',print_debug)) then
        WRITE(luout,*)
     &      'Hessian after projecting out external modes:'
        WRITE(luout,*)HESS
        WRITE(luout,*)
     &      'difference of Hessian with projected hessian:'
        WRITE(luout,*)HESSP
      endif
C****
C**** construct Hessian in trans-rot subspace (should be zeroes)
C****
      DO 230 MU=1,6
        DO 240 NU=1,6
          TEST(MU,NU) = ZERO
          DO 250 K=1,NAT3
            DO 260 L=1,NAT3
              TEST(MU,NU) = TEST(MU,NU)+VC(K,MU)*HESS(K,L)*VC(L,NU)
  260       CONTINUE
  250     CONTINUE
  240   CONTINUE
  230 CONTINUE
      if (util_print('eckart',print_medium)) then
        test_norm = ddot((6*6),test,1,test,1)
        write(luout,'(a,1pd10.4)')
     &      ' Projected Nuclear Hessian trans-rot subspace norm:',
     &      test_norm
        write(luout,*)
     &      '                        (should be close to zero!) '
      endif
      if (util_print('eckart',print_debug)) then
        WRITE(6,*)'Hessian projected into trans-rot subspace ',
     &      '(should be zeros):'
        WRITE(6,*)TEST
      endif
      RETURN
      END


*******************************************************************


      SUBROUTINE vib_eckart_trans( HESS, HESSP, HESST, COORD, VC ,
     ,     w1,w2)
      IMPLICIT none ! REAL*8 (A-H,O-Z)
c
#include "util.fh"
#include "stdio.fh"
c
      double precision zero, one
      PARAMETER( ZERO=0.D0, ONE=1.D0 )
c
      integer nat, nat3, nhess, nhesst
      COMMON /cvib_HESS/ NAT,NAT3,NHESS,NHESST
c
      double precision HESS(NAT3,NAT3)
      double precision HESSP(NAT3,NAT3)
      double precision HESST(NAT3*(NAT3+1)/2)
      double precision COORD(3,NAT)
      double precision VC(NAT3,3)
      double precision w1(nat3,*),w2(nat3,*)
c
      DOUBLE PRECISION UNIVEC(3), TEST(3,3), VNORM, temp, dotval, rnorm
      double precision test_norm
      integer i, j, k, l, m, n, mu, nu, indx, iatom, iaxis, itemp
C****
C**** construct translation unit vectors;  these are stored in the
C**** first three columns of array VC, the rotation vectors will
C**** be stored in the other 3 columns
C****
      CALL Dfill (3*NAT3,0.0d00,VC,1)
      VNORM = ONE/SQRT(dble(NAT))
      DO 10 I=1,3 ! unit translation vector for x-, y-, and z-directions
        DO 20 IATOM=1,NAT
           VC(3*(IATOM-1)+I,I) = VNORM
   20   CONTINUE
   10 CONTINUE

cC****
cC**** construct rotation unit vectors; store in last 3 columns of VC
cC****
c      DO 40 IAXIS=1,3         ! loop over rotation axes
c        DO 50 J=1,3
c          IF (J.EQ.IAXIS) THEN
c            UNIVEC(J) = ONE
c          ELSE
c            UNIVEC(J) = ZERO
c          END IF
c   50   CONTINUE
c        DO 60 IATOM=1,NAT       ! loop over atoms
c          CALL vib_CROSS(COORD(1,IATOM),UNIVEC,
c     &      VC(3*(IATOM-1)+1,IAXIS+3))
c   60   CONTINUE
c   40 CONTINUE
cC****
C****   Schmidt orthogonalize the constraint vectors
C****
      DO 70 I=2,3     ! orthogonalize vector I to each of vectors J
        DO 80 J=1,I-1
           dotval=ddot(nat3,vc(1,J),1,vc(1,i),1)
          DO 100 K=1,NAT3
             VC(K,I) = VC(K,I)-dotval*VC(K,J)
  100     CONTINUE
   80   CONTINUE
        RNORM = ddot(nat3,vc(1,I),1,vc(1,I),1)
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
          call dscal(nat3,rnorm,vc(1,I),1)
        ELSE ! dependent constraint vector (linear molecule; set to zero)
          DO 125 K=1,NAT3
            VC(K,I) = ZERO
  125     CONTINUE
          ITEMP = I-3
          WRITE(luout,*)
     &        'Dependent rotation vector no.',ITEMP
          WRITE(luout,*)
     &        ' found in ECKART; assuming linear geometry'
        END IF
   70 CONTINUE
C****
C**** remove rotations and translations from the Hessian matrix
C****
!      DO 130 M=1,NAT3
!        DO 140 N=1,M
!          HESSP(M,N) = HESS(M,N)
!          DO 150 MU=1,6
!            DO 160 J=1,NAT3
!              HESSP(M,N) = HESSP(M,N) !- VC(M,MU)*VC(J,MU)*HESS(J,N)
!     &                               !- HESS(M,J)*VC(J,MU)*VC(N,MU)
!              DO 170 NU=1,6
!                DO 180 K=1,NAT3
!                  HESSP(M,N) = HESSP(M,N) +
!     &                VC(M,MU)*VC(J,MU)*HESS(J,K)*VC(K,NU)*VC(N,NU)
!  180           CONTINUE
!  170         CONTINUE
!  160       CONTINUE
!  150     CONTINUE
!          if(M.ne.N) HESSP(N,M) = HESSP(M,N)
!  140   CONTINUE
!  130 CONTINUE
      call dcopy(nat3*nat3,hess,1,hessp,1)

      call dgemm('n','t',nat3,nat3,3,1d0,vc,nat3,vc,nat3,
     .     0d0,w1,nat3)
      call dgemm('n','n',nat3,nat3,nat3,1d0,w1,nat3,hess,nat3,
     .     0d0,w2,nat3)
      call dgemm('n','n',nat3,nat3,nat3,1d0,w2,nat3,w1,nat3,
     .     1d0,hessp,nat3)
      call dgemm('n','n',nat3,nat3,nat3,-1d0,w1,nat3,hess,nat3,
     .     1d0,hessp,nat3)
      call dgemm('n','n',nat3,nat3,nat3,-1d0,hess,nat3,w1,nat3,
     .     1d0,hessp,nat3)
C****
C**** update triangular hessian to agree with HESS
C****
      DO 190 I=1,NAT3
        DO 200 J=1,I
          indx = i*(i-1)/2 + j
          HESST(indx) = HESSP(I,J)
  200   CONTINUE
  190 CONTINUE
C****
C**** place the projected Hessian in array HESS, and the difference
C**** between projected and unprojected in HESSP
C****
      DO 210 I=1,NAT3
        DO 220 J=1,I
          TEMP = HESS(I,J)
          HESS(I,J) = HESSP(I,J)
          HESSP(I,J) = HESS(I,J)-TEMP
          IF (DABS(HESSP(I,J)).LT.1.D-9) HESSP(I,J) = ZERO
          IF(J.NE.I) then
             HESSP(J,I)=HESSP(I,J)
             HESS(J,I)=HESS(I,J)
          endif
  220   CONTINUE
  210 CONTINUE
C****
C**** output projected Hessian
C****
      if (util_print('eckart',print_debug)) then
        WRITE(luout,*)
     &      'Hessian after projecting out external modes:'
        WRITE(luout,*)HESS
        WRITE(luout,*)
     &      'difference of Hessian with projected hessian:'
        WRITE(luout,*)HESSP
      endif
C****
C**** construct Hessian in trans-rot subspace (should be zeroes)
C****
      DO 230 MU=1,3
        DO 240 NU=1,3
          TEST(MU,NU) = ZERO
          DO 250 K=1,NAT3
            DO 260 L=1,NAT3
              TEST(MU,NU) = TEST(MU,NU)+VC(K,MU)*HESS(K,L)*VC(L,NU)
  260       CONTINUE
  250     CONTINUE
  240   CONTINUE
  230 CONTINUE
      if (util_print('eckart',print_medium)) then
        test_norm = ddot((3*3),test,1,test,1)
        write(luout,'(a,1pd10.4)')
     &      ' Projected Nuclear Hessian trans-rot subspace norm:',
     &      test_norm
        write(luout,*)
     &      '                        (should be close to zero!) '
      endif
      if (util_print('eckart',print_debug)) then
        WRITE(luout,*)'Hessian projected into trans-rot subspace ',
     &      '(should be zeros):'
        WRITE(luout,*)TEST
      endif
      RETURN
      END
