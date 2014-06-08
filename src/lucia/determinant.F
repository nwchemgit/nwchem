      SUBROUTINE DTRM(A,N,D,INDX)
C
C Subroutine for evaluating the determinant of a matrix using 
C the partial-pivoting Gaussian elimination scheme.
C
      DIMENSION A(N,N),INDX(N)
C
      CALL ELGS(A,N,INDX)
C
      D    = 1.0
      DO     100 I = 1, N
         D = D*A(INDX(I),I)
  100 CONTINUE
C
      MSGN = 1
      DO     200 I = 1, N
        DO   150 WHILE (I.NE.INDX(I))
          MSGN = -MSGN
          J = INDX(I)
          INDX(I) = INDX(J)
          INDX(J) = J
  150   END DO
  200 CONTINUE
      D = MSGN*D
C
      RETURN
      END
C
      SUBROUTINE ELGS(A,N,INDX)
C
C Subroutine to perform the partial-pivoting Gaussian elimination.
C A(N,N) is the original matrix in the input and transformed
C matrix plus the pivoting element ratios below the diagonal in
C the output.  INDX(N) records the pivoting order.

C
      DIMENSION A(N,N),INDX(N),C(N)
C
C Initialize the index
C
      DO     50    I = 1, N
        INDX(I) = I
   50 CONTINUE
C
C Find the rescaling factors, one from each row
C
        DO     100   I = 1, N
          C1= 0.0
          DO    90   J = 1, N
            C1 = AMAX1(C1,ABS(A(I,J)))
   90     CONTINUE
          C(I) = C1
  100   CONTINUE
C
C Search the pivoting (largest) element from each column
C
      DO     200   J = 1, N-1
        PI1 = 0.0
        DO   150   I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
  150   CONTINUE
C
C Interchange the rows via INDX(N) to record pivoting order
C
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO   170   I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)
C
C Record pivoting ratios below the diagonal
C
          A(INDX(I),J) = PJ
C
C Modify other elements accordingly
C
          DO 160   K = J+1, N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
  160     CONTINUE
  170   CONTINUE
  200 CONTINUE
C
      RETURN
      END
c $Id$
