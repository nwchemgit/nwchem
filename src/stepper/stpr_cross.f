      SUBROUTINE stpr_cross(V1,V2,V3)
c $Id: stpr_cross.f,v 1.1 1994-06-24 20:42:18 d3e129 Exp $
C
C     form cross product v1 x v2.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      DIMENSION V1(3), V2(3), V3(3)
      V3(1) = V1(2)*V2(3)-V1(3)*V2(2)
      V3(2) = V1(3)*V2(1)-V1(1)*V2(3)
      V3(3) = V1(1)*V2(2)-V1(2)*V2(1)
      RETURN
      END
