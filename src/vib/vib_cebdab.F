      SUBROUTINE vib_CEBDAB(C,NC,A,NA,B)
C* MODULE ZMATRX *DECK CEBDAB
* $Id$
C
C       +
C  C = B A B
C
C  C square matrix target (not same memory location as A)
C  B transformation matrix B(order of A,order of C)
C  A square matrix to be transformed
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
      integer nc, na
      DOUBLE PRECISION C(NC,NC)
      DOUBLE PRECISION A(NA,NA)
      DOUBLE PRECISION B(NA,NC)
c
      integer i, j, l, k
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
