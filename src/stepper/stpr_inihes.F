      SUBROUTINE stpr_inihes ( HESIAN )
c $Id: stpr_inihes.F,v 1.2 1994-10-26 23:11:30 d3e129 Exp $
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
*debug:      write(6,*)'hesian stpr_inihes', nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hesian(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hesian (',i,') =',hesian(i)
*debug:      enddo
      RETURN
      END
