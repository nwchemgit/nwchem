*
* $Id$
*
      SUBROUTINE EIGEN(NEMAX,NE,HML,EIG,WORK,IERR)
*     =======================================
*     DIAGONALIZATION OF HAMILTONIAN MATRIX 
*     =======================================
      implicit double precision(a-h, o-z)
      DIMENSION EIG(*),HML(NEMAX,*),WORK(*)

      CALL DCOPY(NEMAX,0.0d0,0,EIG,1)
      IF(NE.EQ.0) RETURN
      IF(NE.EQ.1) THEN
        EIG(1)=HML(1,1)
        HML(1,1)=1.0d0
        RETURN
      ENDIF
      CALL TRED3(HML,NE,NEMAX,EIG,WORK)
      CALL TQLI(EIG,WORK,NE,NEMAX,HML,IERR)

      CALL EIGSRT(EIG,HML,NE,NEMAX)

 1000 FORMAT(I8,E20.10)
      RETURN
      END

*     =======================================
*     DIAGONALIZATION OF HAMILTONIAN MATRIX 
*     =======================================
      SUBROUTINE EIGEN_UNSORT(NEMAX,NE,HML,EIG,WORK,IERR)
      implicit double precision(a-h, o-z)
      DIMENSION EIG(*),HML(NEMAX,*),WORK(*)

      CALL DCOPY(NEMAX,0.0d0,0,EIG,1)
      IF(NE.EQ.0) RETURN
      IF(NE.EQ.1) THEN
        EIG(1)=HML(1,1)
        HML(1,1)=1.0d0
        RETURN
      ENDIF
      CALL TRED3(HML,NE,NEMAX,EIG,WORK)
      CALL TQLI(EIG,WORK,NE,NEMAX,HML,IERR)

      RETURN
      END




      SUBROUTINE EIGSRT(D,V,N,NP)
*     -----------------------------------------------------------------
*     SORT EIGENVALUES INTO DESCENDING ORDER AND REARRANGE EIGENVECTORS
*     -----------------------------------------------------------------
      implicit double precision(a-h, o-z)
      double precision D(NP),V(NP,NP)
!$OMP MASTER
      DO 13 I=1,N-1
         K=I
         P=D(I)
         DO 11 J=I+1,N
            IF(D(J).GE.P) THEN
               K=J
               P=D(J)
            ENDIF
   11    CONTINUE
         IF(K.NE.I) THEN
            D(K)=D(I)
            D(I)=P
            DO 12 J=1,N
               P=V(J,I)
               V(J,I)=V(J,K)
               V(J,K)=P
   12       CONTINUE
         ENDIF
   13 END DO
!$OMP END MASTER
!$OMP BARRIER
      RETURN
      END


      SUBROUTINE DMADD(NX,N,X,Y,Z)
*     ==============
*     MATRIX ADDTION
*     ==============
      implicit double precision(a-h, o-z)
      DIMENSION X(*),Y(*),Z(*)
      DO 1 I=1,NX*N
        Z(I)=X(I)+Y(I)
    1 CONTINUE
      RETURN
      END

      SUBROUTINE DMMUL(NX,N,X,Y,Z)
*     =====================
*     MATRIX MULTIPLICATION
*     =====================
      implicit double precision(a-h, o-z)
      DIMENSION X(NX,*),Y(NX,*),Z(NX,*)
      DO 1 J=1,N
      DO 1 I=1,N
        Z(I,J)=DDOT(N,X(I,1),NX,Y(1,J),1)
    1 CONTINUE
      RETURN
      END

      SUBROUTINE DMSUB(NX,N,X,Y,Z)
*     ===================
*     MUTRIX SUBSTRACTION
*     ===================
      implicit double precision(a-h, o-z)
      DIMENSION X(*),Y(*),Z(*)
      DO 1 I=1,NX*N
        Z(I)=X(I)-Y(I)
    1 CONTINUE
      RETURN
      END


      SUBROUTINE TQLI(D,E,N,NP,Z,IERR)
*     ==========================================
*     EIGENVALUES AND EIGENVECTORS
*     FOR A REAL SYMMETRIC TRIDIAGONAL MATRIX 
*     ==========================================
      implicit double precision(a-h, o-z)
      double precision D(NP),E(NP),Z(NP,NP)
      IF(N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
   11   CONTINUE
        E(N)=0
        DO 15 L=1,N
          ITER=0
    1     DO 12 M=L,N-1
            DD=DABS(D(M))+DABS(D(M+1))
            IF(DABS(E(M))+DD.EQ.DD) GO TO 2
   12     CONTINUE
          M=N
    2     IF(M.NE.L) THEN
            IF(ITER.EQ.300) THEN
              IERR=1
              RETURN
            ENDIF
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.0d0*E(L))
            R=DSQRT(G**2+1.0d0)
            G=D(M)-D(L)+E(L)/(G+DSIGN(R,G))
            S=1.0d0
            C=1.0d0
            P=0.0d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(DABS(F).GE.DABS(G)) THEN
                C=G/F
                R=DSQRT(C**2+1.0d0)
                E(I+1)=F*R
                S=1.0d0/R
                C=C*S
              ELSE
                S=F/G
                R=DSQRT(S**2+1.0d0)
                E(I+1)=G*R
                C=1.0d0/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.0d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
   13         CONTINUE
   14       CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.0d0
            GO TO 1
          ENDIF
   15   CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE TRED3(A,N,NP,D,E)
*     ==================================================
*     HOUSEHOLDER REDUCTION OF A REAL SYMMETRIC MATRIX 
*     ==================================================
      implicit double precision(a-h, o-z)
      double precision A(NP,NP),D(NP),E(NP)
      IF(N.GT.1) THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.0d0
          SCALE=0.0d0
          IF(L.GT.1) THEN
            DO 11 K=1,L
              SCALE=SCALE+DABS(A(I,K))
   11       CONTINUE
            IF(SCALE.EQ.0.) THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
   12         CONTINUE
              F=A(I,L)
              G=-DSIGN(DSQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.0d0
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.0d0
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
   13           CONTINUE
                IF(L.GT.J) THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
   14             CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
   15         CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
   16           CONTINUE
   17         CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
   18   CONTINUE
      ENDIF
      D(1)=0.0d0
      E(1)=0.0d0
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.0d0) THEN
          DO 21 J=1,L
            G=0.0d0
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
   19       CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
   20       CONTINUE
   21     CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.0d0
        IF(L.GE.1) THEN
          DO 22 J=1,L
            A(I,J)=0.0d0
            A(J,I)=0.0d0
   22     CONTINUE
        ENDIF
   23 CONTINUE
      RETURN
      END
