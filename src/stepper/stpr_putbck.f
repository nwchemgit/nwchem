      SUBROUTINE stpr_putbck ( COORD, C )
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / CFACE / IWCTR,NATOM,ICALC
      DIMENSION COORD(3,NATOM), C(3,NATOM)
C
C     Get coordinates.
C
      DO 2 I = 1, 3
        DO 1 N = 1, NATOM
          C(I,N) = COORD(I,N)
    1   CONTINUE
    2 CONTINUE
      RETURN
      END
