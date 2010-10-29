      SUBROUTINE YGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
c $Id$
*
      implicit none
#include "y64.fh"
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
c
      INTGR4            IHI4, ILO4, INFO4, LDA4, N4      
c
      LDA4=lda
      N4=n
c
      call DGEBAL( JOB, N4, A, LDA4, ILO4, IHI4, SCALE,INFO4)
      IHI=ihi4
      ILO=ilo4
c
      info4=info
      return
      end
