      SUBROUTINE vib_eckart( HESS, HESSP, HESST, COORD, VC )
* $Id: vib_eckart.F,v 1.1 1995-07-29 23:19:29 d3e129 Exp $
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER( ZERO=0.D0, ONE=1.D0 )
      COMMON /HESS/ NATOM,NAT3,NHESS,NHESST
      DIMENSION HESS(NAT3,NAT3), COORD(3,NATOM), VC(NAT3,6),
     &  UNIVEC(3), HESST(NAT3*(NAT3+1)/2), TEST(6,6),
     &  HESSP(NAT3,NAT3)
C****
C**** construct translation unit vectors;  these are stored in the
C**** first three columns of array VC, the rotation vectors will
C**** be stored in the other 3 columns
C****
      CALL Dfill (6*NAT3,0.0d00,VC,1)
      VNORM = ONE/SQRT(dble(NATOM))
      DO 10 I=1,3 ! unit translation vector for x-, y-, and z-directions
        DO 20 IATOM=1,NATOM
          DO 30 IAXIS=1,3
            IF( IAXIS.NE.I ) THEN
              VC(3*(IATOM-1)+IAXIS,I) = ZERO
            ELSE
              VC(3*(IATOM-1)+IAXIS,I) = VNORM
            END IF
   30     CONTINUE
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
        DO 60 IATOM=1,NATOM       ! loop over atoms
          CALL CROSS(COORD(1,IATOM),UNIVEC,
     &      VC(3*(IATOM-1)+1,IAXIS+3))
   60   CONTINUE
   40 CONTINUE
C****
C****   Schmidt orthogonalize the constraint vectors
C****
      DO 70 I=2,6     ! orthogonalize vector I to each of vectors J
        DO 80 J=1,I-1
          DOT = ZERO
          DO 90 K=1,NAT3
            DOT = DOT+VC(K,J)*VC(K,I)
   90     CONTINUE
          DO 100 K=1,NAT3
            VC(K,I) = VC(K,I)-DOT*VC(K,J)
  100     CONTINUE
   80   CONTINUE
        RNORM = ZERO  ! normalize vector I
        DO 110 K=1,NAT3
          RNORM = RNORM+VC(K,I)*VC(K,I)
  110   CONTINUE
        IF ( RNORM.GT.1.D-8 ) THEN
          RNORM = ONE/SQRT(RNORM)
          DO 120 K=1,NAT3
            VC(K,I) = VC(K,I)*RNORM
  120     CONTINUE
        ELSE ! dependent constraint vector (linear molecule; set to zero)
          DO 125 K=1,NAT3
            VC(K,I) = ZERO
  125     CONTINUE
          ITEMP = I-3
          WRITE(6,*)'Dependent rotation vector no.',ITEMP
          WRITE(6,*)' found in ECKART; assuming linear geometry'
        END IF
   70 CONTINUE
C****
C**** remove rotations and translations from the Hessian matrix
C****
      DO 130 M=1,NAT3
        DO 140 N=1,NAT3
          HESSP(M,N) = HESS(M,N)
          DO 150 MU=1,6
            DO 160 J=1,NAT3
              HESSP(M,N) = HESSP(M,N) - VC(M,MU)*VC(J,MU)*HESS(J,N)
     &                               - HESS(M,J)*VC(J,MU)*VC(N,MU)
              DO 170 NU=1,6
                DO 180 K=1,NAT3
                  HESSP(M,N) = HESSP(M,N) +
     &                VC(M,MU)*VC(J,MU)*HESS(J,K)*VC(K,NU)*VC(N,NU)
  180           CONTINUE
  170         CONTINUE
  160       CONTINUE
  150     CONTINUE
  140   CONTINUE
  130 CONTINUE
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
        DO 220 J=1,NAT3
          TEMP = HESS(I,J)
          HESS(I,J) = HESSP(I,J)
          HESSP(I,J) = HESS(I,J)-TEMP
          IF (DABS(HESSP(I,J)).LT.1.D-9) HESSP(I,J) = ZERO
  220   CONTINUE
  210 CONTINUE
C****
C**** output projected Hessian
C****
*      WRITE(6,*)'Hessian after projecting out external modes:'
*      WRITE(6,*)HESS
C      WRITE(6,*)'difference of Hessian with projected hessian:'
C      WRITE(6,*)HESSP
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
C      WRITE(6,*)'Hessian projected into trans-rot subspace:'
C      WRITE(6,*)TEST
      RETURN
      END
