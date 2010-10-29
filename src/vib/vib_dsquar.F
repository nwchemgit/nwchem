      SUBROUTINE vib_DSQUAR(A,S,N,MS)
* $Id$
C
C EXPANSION OF LOWER TRIANGULAR PACKED ARRAY TO FULL SQUARE
C INPUT:
C     A   LOWER TRIANGULAR PACKED ARRAY
C     N   ROW ORDER OF A
C     MS  ROW DIMENSION OF S IN CALLING ROUTINE
C OUTPUT:
C     S  SQUARE PACKED ARRAY CORRESPONDING TO A
C
      implicit none
      integer n,ms
      double precision S(*),A(*),VAL
c
      integer ij,ni,nj,ijs,jis,i,j
c
      IJ=0
      NI=0
      DO 11 I=1,N
        NJ=0
        DO 10 J=1,I
          IJ=IJ+1
          VAL=A(IJ)
          IJS=I+NJ
          JIS=J+NI
          S(IJS)=VAL
          S(JIS)=VAL
          NJ=NJ+MS
   10   CONTINUE
        NI=NI+MS
   11 CONTINUE
      RETURN
      END
