      SUBROUTINE vib_CEBDAB(C,NC,A,NA,B)
C* MODULE ZMATRX *DECK CEBDAB
* $Id: vib_cebdab.F,v 1.1 1995-07-29 23:19:22 d3e129 Exp $
C
C       +
C  C = B A B
C
C  C square matrix target (not same memory location as A)
C  B transformation matrix B(order of A,order of C)
C  A square matrix to be transformed
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(NC,NC)
      DIMENSION A(NA,NA)
      DIMENSION B(NA,NC)
C
      DO 00100 I = 1,NC
         DO 00100 J = 1,NC
            C(I,J) = 0.0D00
            DO 00200 L = 1,NA
               DO 00200 K = 1,NA
                  C(I,J) = C(I,J) + B(L,I)*B(K,J)*A(L,K)
00200       CONTINUE
00100 CONTINUE
C
      RETURN
      END
