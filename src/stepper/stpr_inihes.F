      SUBROUTINE stpr_inihes ( HESIAN )
c $Id: stpr_inihes.F,v 1.1 1994-06-24 20:42:26 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
#include "chow.fh"
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION HESIAN(NAT3TR)
      call dfill (NAT3TR,0.0d00, HESIAN,1)
      IJ = 0
      DO 2 I = 1,NAT3
        DO 1 J = 1,I
          IJ = IJ + 1
          IF (I.EQ.J) THEN
            IF (MWC) THEN
              HESIAN(IJ) = 4.0D-04
            ELSE
              HESIAN(IJ) = ONE
            ENDIF
          ELSE
            HESIAN(IJ) = ZERO
          ENDIF
    1   CONTINUE
    2 CONTINUE
      RETURN
      END
