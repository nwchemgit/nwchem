      SUBROUTINE stpr_mktri (H, HTRI)
c $Id$
C
C     opposite of BLOWUP; contract square symmetric matrix
C     to lower triangular form.
C
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION H(NAT3,NAT3), HTRI(NAT3TR)
      IJ = 0
      DO 1 I=1,3*NATOM
        DO 2 J=1,I
          IJ = IJ+1
          HTRI(IJ) = H(J,I)
    2   CONTINUE
    1 CONTINUE
      RETURN
      END
