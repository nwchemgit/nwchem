      SUBROUTINE vib_TRANSP(A,NDIM,N)
C*MODULE VIB  *DECK TRANSP
* $Id: vib_transp.F,v 1.2 1995-10-21 09:16:28 d3e129 Exp $
      IMPLICIT NONE
      integer ndim, n
      double precision A(NDIM,NDIM)
c
      double precision temp
      integer i, j
C
C  TRANSPOSE MATRIX A, STORE RESULT IN A.
C
      DO 10 I=1,N
        DO 10 J=1,I
          TEMP=A(I,J)
          A(I,J)=A(J,I)
   10     A(J,I)=TEMP
      RETURN
      END
