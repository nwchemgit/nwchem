      SUBROUTINE vib_CROSS(V1,V2,V3)
* $Id$
C****
C**** form cross product v1 x v2
C****
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
      double precision V1(3), V2(3), V3(3)
c
      V3(1) = V1(2)*V2(3)-V1(3)*V2(2)
      V3(2) = V1(3)*V2(1)-V1(1)*V2(3)
      V3(3) = V1(1)*V2(2)-V1(2)*V2(1)
c
      RETURN
      END
