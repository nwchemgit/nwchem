      SUBROUTINE vib_DTRNGL(S,A,N,MS)
* $Id$
C
C COMPRESSION OF SQUARE MATRIX TO LOWER TRIANGULAR PACKING
C INPUT:
C     S     SQUARE PACKED ARRAY
C     N     ORDER OF S
C     MS    ROW DIMENSION OF S IN CALLING ROUTINE
C OUTPUT:
C     A    TRIANGULAR PACKED ARRAY
C
      implicit none
      integer N, MS
      double precision S(*),A(*)
c
      integer ij,ijs, i, j, nj
c
      IJ=0
      DO 11 I=1,N
        NJ=0
        DO 10 J=1,I
          IJ=IJ+1
          IJS=I+NJ
          A(IJ)=S(IJS)
          NJ=NJ+MS
   10   CONTINUE
   11 CONTINUE
      RETURN
      END
